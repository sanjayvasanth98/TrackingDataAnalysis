function plot_proximity_activation_analysis(allProximityActivation, outDir, plotOpts, opts)
%PLOT_PROXIMITY_ACTIVATION_ANALYSIS  Publication plots for neighbor activation.

if nargin < 3 || isempty(plotOpts)
    plotOpts = struct();
end
if ~isfield(plotOpts, 'themes') || isempty(plotOpts.themes)
    plotOpts.themes = "normal";
end
if nargin < 4 || isempty(opts)
    opts = struct();
end
opts = apply_plot_defaults(opts);

if ~isfolder(outDir)
    mkdir(outDir);
end

[pairTable, eventTable, binnedStats] = proximity_activation_to_tables(allProximityActivation);
if isempty(pairTable) || height(pairTable) == 0
    warning('plot_proximity_activation_analysis: no neighbor pairs. Skipping plots.');
    return;
end

plot_standoff_response(pairTable, binnedStats, outDir, plotOpts, opts);
plot_secondary_probability(binnedStats, outDir, plotOpts, opts);
plot_radius_rate_coupling(pairTable, binnedStats, outDir, plotOpts, opts);
plot_nearest_survival_cdf(eventTable, outDir, plotOpts, opts);
fprintf('Saved proximity activation plots to: %s\n', outDir);
end


% =========================================================================
function opts = apply_plot_defaults(opts)
opts = default_field(opts, 'maxGamma', 20);
opts = default_field(opts, 'extendedMaxGamma', 40);
opts = default_field(opts, 'gammaBins', [0 2 5 10 20 40]);
opts = default_field(opts, 'responseYLim', []);
opts = default_field(opts, 'rdotRatioYLim', []);
opts = default_field(opts, 'minPairsForTrend', 3);
end


% =========================================================================
function s = default_field(s, name, value)
if ~isfield(s, name) || isempty(s.(name))
    s.(name) = value;
end
end


% =========================================================================
function plot_standoff_response(pairTable, ~, outDir, plotOpts, opts)
ReVals = unique(double(pairTable.Re));
ReVals = ReVals(isfinite(ReVals));
if isempty(ReVals), return; end

for theme = reshape(plotOpts.themes, 1, [])
    fontName = resolve_plot_font_name();
    [nRows, nCols] = subplot_grid(numel(ReVals));
    f = figure('Color', 'w', 'Position', [80 80 620*nCols 560*nRows]);

    for ri = 1:numel(ReVals)
        ax = subplot(nRows, nCols, ri, 'Parent', f);
        hold(ax, 'on');
        Rei = ReVals(ri);
        reMask = double(pairTable.Re) == Rei;
        caseKeys = case_key_table(pairTable(reMask, :));
        cmap = scientific_line_colormap(height(caseKeys));
        lgd = gobjects(0, 1);
        lgdTxt = strings(0, 1);

        randomVals = double(pairTable.randomMaxAbsDeltaR_over_R0(reMask));
        randomVals = randomVals(isfinite(randomVals) & randomVals >= 0);
        if numel(randomVals) >= 3
            medRandom = finite_prctile(randomVals, 50);
            p90Random = finite_prctile(randomVals, 90);
            patch(ax, [0 opts.extendedMaxGamma opts.extendedMaxGamma 0], ...
                [medRandom medRandom p90Random p90Random], [0.82 0.82 0.82], ...
                'FaceAlpha', 0.34, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            plot(ax, [0 opts.extendedMaxGamma], [p90Random p90Random], ':', ...
                'Color', [0.35 0.35 0.35], 'LineWidth', 1.3, 'HandleVisibility', 'off');
        end

        for ci = 1:height(caseKeys)
            cName = string(caseKeys.Case(ci));
            kD = double(caseKeys.kD(ci));
            [xBin, yMed, yQ25, yQ75, yP90] = deformation_bin_summary(pairTable, reMask, cName, kD, opts);
            if isempty(xBin)
                continue;
            end
            col = cmap(ci, :);
            yNeg = max(yMed - yQ25, 0);
            yPos = max(yQ75 - yMed, 0);
            h = errorbar(ax, xBin, yMed, yNeg, yPos, 'o-', ...
                'Color', col, ...
                'MarkerFaceColor', lighten_color(col, 0.12), ...
                'MarkerEdgeColor', [0.05 0.05 0.05], ...
                'MarkerSize', 6.5, ...
                'CapSize', 8, ...
                'LineWidth', 2.2);
            lgd(end+1,1) = h; %#ok<AGROW>
            lgdTxt(end+1,1) = sprintf('k/d = %.4g', kD); %#ok<AGROW>
            plot(ax, xBin, yP90, '--', 'Color', lighten_color(col, 0.25), ...
                'LineWidth', 1.7, 'HandleVisibility', 'off');
        end

        xlim(ax, [0 opts.extendedMaxGamma]);
        if ~isempty(opts.responseYLim) && numel(opts.responseYLim) >= 2
            ylim(ax, opts.responseYLim);
        else
            yAuto = finite_prctile(double(pairTable.maxAbsDeltaR_over_R0(reMask)), 98);
            if isfinite(yAuto) && yAuto > 0
                ylim(ax, [0 1.2*yAuto]);
            end
        end
        draw_gamma_reference_lines(ax, opts);
        xlabel(ax, 'Projected standoff $\gamma=d/R_{\max}$', 'Interpreter', 'latex');
        ylabel(ax, '$\max|\Delta R|/R_0$', 'Interpreter', 'latex');
        title(ax, sprintf('Binned neighbor deformation (median, IQR), Re = %g', Rei), 'FontName', fontName);
        set(ax, 'FontName', fontName);
        grid(ax, 'off');
        box(ax, 'on');
        apply_plot_theme(ax, char(theme));
        if ~isempty(lgd)
            leg = legend(ax, lgd, cellstr(lgdTxt), 'Location', 'southoutside', ...
                'NumColumns', min(2, numel(lgdTxt)), 'Box', 'off');
            style_legend_for_theme(leg, char(theme));
        end
    end

    save_fig_dual_safe(f, fullfile(outDir, "Standoff_Deformation_Map_" + theme), plotOpts);
    close_if_needed(f, plotOpts);
end
end


% =========================================================================
function [xBin, yMed, yQ25, yQ75, yP90, nBin] = deformation_bin_summary(pairTable, reMask, caseName, kD, opts)
xBin = nan(0, 1);
yMed = nan(0, 1);
yQ25 = nan(0, 1);
yQ75 = nan(0, 1);
yP90 = nan(0, 1);
nBin = nan(0, 1);

caseMask = reMask & string(pairTable.Case) == string(caseName) & ...
    abs(double(pairTable.kD) - kD) <= eps(max(1, abs(kD)));
gamma = double(pairTable.gammaForPlot);
y = double(pairTable.maxAbsDeltaR_over_R0);

for bi = 1:(numel(opts.gammaBins) - 1)
    lo = opts.gammaBins(bi);
    hi = opts.gammaBins(bi + 1);
    inBin = caseMask & gamma >= lo & gamma < hi;
    if bi == (numel(opts.gammaBins) - 1)
        inBin = caseMask & gamma >= lo & gamma <= hi;
    end

    vals = y(inBin);
    vals = vals(isfinite(vals) & vals >= 0);
    if numel(vals) < opts.minPairsForTrend
        continue;
    end

    xBin(end+1, 1) = (lo + hi) / 2; %#ok<AGROW>
    yMed(end+1, 1) = finite_median(vals); %#ok<AGROW>
    yQ25(end+1, 1) = finite_prctile(vals, 25); %#ok<AGROW>
    yQ75(end+1, 1) = finite_prctile(vals, 75); %#ok<AGROW>
    yP90(end+1, 1) = finite_prctile(vals, 90); %#ok<AGROW>
    nBin(end+1, 1) = numel(vals); %#ok<AGROW>
end
end


% =========================================================================
function plot_secondary_probability(binnedStats, outDir, plotOpts, opts)
if isempty(binnedStats) || height(binnedStats) == 0
    return;
end
ReVals = unique(double(binnedStats.Re));
ReVals = ReVals(isfinite(ReVals));
for theme = reshape(plotOpts.themes, 1, [])
    fontName = resolve_plot_font_name();
    [nRows, nCols] = subplot_grid(numel(ReVals));
    f = figure('Color', 'w', 'Position', [80 80 620*nCols 560*nRows]);

    for ri = 1:numel(ReVals)
        ax = subplot(nRows, nCols, ri, 'Parent', f);
        hold(ax, 'on');
        Rei = ReVals(ri);
        reMask = double(binnedStats.Re) == Rei;
        caseKeys = case_key_table(binnedStats(reMask, :));
        cmap = scientific_line_colormap(height(caseKeys));
        lgd = gobjects(0, 1);
        lgdTxt = strings(0, 1);

        for ci = 1:height(caseKeys)
            cName = string(caseKeys.Case(ci));
            kD = double(caseKeys.kD(ci));
            T = case_binned_rows(binnedStats, cName, Rei, kD);
            if isempty(T) || height(T) == 0
                continue;
            end
            [~, ord] = sort(double(T.gammaBinCenter));
            T = T(ord, :);
            x = double(T.gammaBinCenter);
            p = 100 * double(T.secondaryActivationProbability);
            ciLow = 100 * double(T.secondaryActivationCI_low);
            ciHigh = 100 * double(T.secondaryActivationCI_high);
            eLow = max(0, p - ciLow);
            eHigh = max(0, ciHigh - p);
            col = cmap(ci, :);
            h = errorbar(ax, x, p, eLow, eHigh, 'o-', ...
                'Color', col, 'MarkerFaceColor', col, 'MarkerEdgeColor', [0.05 0.05 0.05], ...
                'LineWidth', 1.6, 'MarkerSize', 8, 'CapSize', 7);
            lgd(end+1,1) = h; %#ok<AGROW>
            lgdTxt(end+1,1) = sprintf('%s, k/d=%.4g', cName, kD); %#ok<AGROW>
            for j = 1:numel(x)
                text(ax, x(j), p(j), sprintf(' n=%d', T.nPairs(j)), ...
                    'FontName', fontName, 'FontSize', 8, ...
                    'Color', theme_text_color(char(theme)), ...
                    'VerticalAlignment', 'bottom');
            end
        end

        xlim(ax, [0 opts.extendedMaxGamma]);
        ylim(ax, [0 100]);
        draw_gamma_reference_lines(ax, opts);
        xlabel(ax, 'Projected standoff bin center, $\gamma$', 'Interpreter', 'latex');
        ylabel(ax, 'Secondary activation probability (\%)', 'Interpreter', 'latex');
        title(ax, sprintf('Viable-neighbor cavitation probability, Re = %g', Rei), 'FontName', fontName);
        set(ax, 'FontName', fontName);
        grid(ax, 'off');
        box(ax, 'on');
        apply_plot_theme(ax, char(theme));
        if ~isempty(lgd)
            leg = legend(ax, lgd, cellstr(lgdTxt), 'Location', 'southoutside', ...
                'NumColumns', min(2, numel(lgdTxt)), 'Box', 'off');
            style_legend_for_theme(leg, char(theme));
        end
    end
    save_fig_dual_safe(f, fullfile(outDir, "Secondary_Cavitation_Probability_vs_Standoff_" + theme), plotOpts);
    close_if_needed(f, plotOpts);
end
end


% =========================================================================
function plot_radius_rate_coupling(pairTable, ~, outDir, plotOpts, opts)
ReVals = unique(double(pairTable.Re));
ReVals = ReVals(isfinite(ReVals));
if isempty(ReVals), return; end

for theme = reshape(plotOpts.themes, 1, [])
    fontName = resolve_plot_font_name();
    [nRows, nCols] = subplot_grid(numel(ReVals));
    f = figure('Color', 'w', 'Position', [80 80 620*nCols 560*nRows]);

    for ri = 1:numel(ReVals)
        ax = subplot(nRows, nCols, ri, 'Parent', f);
        hold(ax, 'on');
        Rei = ReVals(ri);
        reMask = double(pairTable.Re) == Rei;
        caseKeys = case_key_table(pairTable(reMask, :));
        cmap = scientific_line_colormap(height(caseKeys));
        lgd = gobjects(0, 1);
        lgdTxt = strings(0, 1);

        for ci = 1:height(caseKeys)
            cName = string(caseKeys.Case(ci));
            kD = double(caseKeys.kD(ci));
            [xBin, yMed, yQ25, yQ75] = radius_rate_bin_summary(pairTable, reMask, cName, kD, opts);
            if isempty(xBin)
                continue;
            end
            col = cmap(ci, :);
            yNeg = max(yMed - yQ25, 0);
            yPos = max(yQ75 - yMed, 0);
            h = errorbar(ax, xBin, yMed, yNeg, yPos, 'o-', ...
                'Color', col, ...
                'MarkerFaceColor', lighten_color(col, 0.12), ...
                'MarkerEdgeColor', [0.05 0.05 0.05], ...
                'MarkerSize', 6.5, ...
                'CapSize', 8, ...
                'LineWidth', 2.2);
            lgd(end+1,1) = h; %#ok<AGROW>
            lgdTxt(end+1,1) = sprintf('k/d = %.4g', kD); %#ok<AGROW>
        end

        plot(ax, [0 opts.extendedMaxGamma], [1 1], ':', 'Color', [0.2 0.2 0.2], ...
            'LineWidth', 1.3, 'HandleVisibility', 'off');
        xlim(ax, [0 opts.extendedMaxGamma]);
        set(ax, 'YScale', 'log', 'FontName', fontName);
        if ~isempty(opts.rdotRatioYLim) && numel(opts.rdotRatioYLim) >= 2
            ylim(ax, opts.rdotRatioYLim);
        end
        draw_gamma_reference_lines(ax, opts);
        xlabel(ax, 'Projected standoff $\gamma=d/R_{\max}$', 'Interpreter', 'latex');
        ylabel(ax, '$\mathrm{RMS}(\dot R_n)/\mathrm{RMS}(\dot R_p)$', 'Interpreter', 'latex');
        title(ax, sprintf('Binned radius-rate response (median, IQR), Re = %g', Rei), 'FontName', fontName);
        grid(ax, 'off');
        box(ax, 'on');
        apply_plot_theme(ax, char(theme));
        if ~isempty(lgd)
            leg = legend(ax, lgd, cellstr(lgdTxt), 'Location', 'southoutside', ...
                'NumColumns', min(2, numel(lgdTxt)), 'Box', 'off');
            style_legend_for_theme(leg, char(theme));
        end
    end
    save_fig_dual_safe(f, fullfile(outDir, "Radius_Rate_Coupling_vs_Standoff_" + theme), plotOpts);
    close_if_needed(f, plotOpts);
end
end


% =========================================================================
function [xBin, yMed, yQ25, yQ75, nBin] = radius_rate_bin_summary(pairTable, reMask, caseName, kD, opts)
xBin = nan(0, 1);
yMed = nan(0, 1);
yQ25 = nan(0, 1);
yQ75 = nan(0, 1);
nBin = nan(0, 1);

caseMask = reMask & string(pairTable.Case) == string(caseName) & ...
    abs(double(pairTable.kD) - kD) <= eps(max(1, abs(kD)));
gamma = double(pairTable.gammaForPlot);
y = double(pairTable.rdotRmsRatio);

for bi = 1:(numel(opts.gammaBins) - 1)
    lo = opts.gammaBins(bi);
    hi = opts.gammaBins(bi + 1);
    inBin = caseMask & gamma >= lo & gamma < hi;
    if bi == (numel(opts.gammaBins) - 1)
        inBin = caseMask & gamma >= lo & gamma <= hi;
    end

    vals = y(inBin);
    vals = vals(isfinite(vals) & vals > 0);
    if numel(vals) < opts.minPairsForTrend
        continue;
    end

    xBin(end+1, 1) = (lo + hi) / 2; %#ok<AGROW>
    yMed(end+1, 1) = finite_median(vals); %#ok<AGROW>
    yQ25(end+1, 1) = finite_prctile(vals, 25); %#ok<AGROW>
    yQ75(end+1, 1) = finite_prctile(vals, 75); %#ok<AGROW>
    nBin(end+1, 1) = numel(vals); %#ok<AGROW>
end
end


% =========================================================================
function plot_nearest_survival_cdf(eventTable, outDir, plotOpts, opts)
if isempty(eventTable) || height(eventTable) == 0
    return;
end
ReVals = unique(double(eventTable.Re));
ReVals = ReVals(isfinite(ReVals));
for theme = reshape(plotOpts.themes, 1, [])
    fontName = resolve_plot_font_name();
    [nRows, nCols] = subplot_grid(numel(ReVals));
    f = figure('Color', 'w', 'Position', [80 80 620*nCols 560*nRows]);
    for ri = 1:numel(ReVals)
        ax = subplot(nRows, nCols, ri, 'Parent', f);
        hold(ax, 'on');
        Rei = ReVals(ri);
        reMask = double(eventTable.Re) == Rei;
        caseKeys = case_key_table(eventTable(reMask, :));
        cmap = scientific_line_colormap(height(caseKeys));
        colorLgd = gobjects(0, 1);
        colorLgdTxt = strings(0, 1);
        for ci = 1:height(caseKeys)
            cName = string(caseKeys.Case(ci));
            kD = double(caseKeys.kD(ci));
            mask = reMask & string(eventTable.Case) == cName & double(eventTable.kD) == kD;
            stableGamma = double(eventTable.nearestStableGamma(mask));
            stableGamma = stableGamma(isfinite(stableGamma) & stableGamma >= 0 & stableGamma <= opts.extendedMaxGamma);
            colorHandle = gobjects(0);
            if ~isempty(stableGamma)
                [x, y] = empirical_cdf(stableGamma);
                h = plot(ax, x, y, '-', 'Color', cmap(ci,:), 'LineWidth', 2.4);
                colorHandle = h;
            end
            activatedGamma = double(eventTable.nearestActivatedGamma(mask));
            activatedGamma = activatedGamma(isfinite(activatedGamma) & activatedGamma >= 0 & activatedGamma <= opts.extendedMaxGamma);
            if ~isempty(activatedGamma)
                [x, y] = empirical_cdf(activatedGamma);
                hAct = plot(ax, x, y, '--', 'Color', lighten_color(cmap(ci,:), 0.25), ...
                    'LineWidth', 1.8, 'HandleVisibility', 'off');
                if isempty(colorHandle)
                    colorHandle = hAct;
                end
            end
            if ~isempty(colorHandle) && isgraphics(colorHandle)
                colorLgd(end+1,1) = colorHandle; %#ok<AGROW>
                colorLgdTxt(end+1,1) = sprintf('k/d = %.4g', kD); %#ok<AGROW>
            end
        end
        hStableStyle = plot(ax, nan, nan, '-', 'Color', [0 0 0], 'LineWidth', 2.4);
        hActivatedStyle = plot(ax, nan, nan, '--', 'Color', [0 0 0], 'LineWidth', 1.8);
        xlim(ax, [0 opts.extendedMaxGamma]);
        ylim(ax, [0 1]);
        draw_gamma_reference_lines(ax, opts);
        xlabel(ax, 'Nearest viable neighbor standoff $\gamma$', 'Interpreter', 'latex');
        ylabel(ax, 'Empirical CDF', 'Interpreter', 'latex');
        title(ax, sprintf('CDF of nearest viable neighbor, Re = %g', Rei), 'FontName', fontName);
        set(ax, 'FontName', fontName);
        grid(ax, 'off');
        box(ax, 'on');
        apply_plot_theme(ax, char(theme));
        lgd = [colorLgd; hStableStyle; hActivatedStyle];
        lgdTxt = [colorLgdTxt; "non-activated"; "secondary activated"];
        if ~isempty(lgd)
            leg = legend(ax, lgd, cellstr(lgdTxt), 'Location', 'southoutside', ...
                'NumColumns', min(3, numel(lgdTxt)), 'Box', 'off');
            style_legend_for_theme(leg, char(theme));
        end
    end
    save_fig_dual_safe(f, fullfile(outDir, "Nearest_Viable_Survival_CDF_" + theme), plotOpts);
    close_if_needed(f, plotOpts);
end
end


% =========================================================================
function draw_gamma_reference_lines(ax, opts)
yLim = ylim(ax);
refs = [5 10 20];
for i = 1:numel(refs)
    if refs(i) > 0 && refs(i) < opts.extendedMaxGamma
        plot(ax, [refs(i) refs(i)], yLim, ':', ...
            'Color', [0.25 0.25 0.25], 'LineWidth', 1.1, 'HandleVisibility', 'off');
    end
end
ylim(ax, yLim);
end


% =========================================================================
function T = case_key_table(Tin)
if isempty(Tin) || height(Tin) == 0
    T = table();
    return;
end
caseNames = string(Tin.Case);
kD = double(Tin.kD);
[~, ia] = unique([caseNames, string(kD)], 'rows', 'stable');
T = table(caseNames(ia), kD(ia), 'VariableNames', {'Case','kD'});
[~, ord] = sort(T.kD);
T = T(ord, :);
end


% =========================================================================
function T = case_binned_rows(binnedStats, caseName, Re, kD)
T = table();
if isempty(binnedStats) || height(binnedStats) == 0
    return;
end
mask = string(binnedStats.Case) == string(caseName) & ...
    double(binnedStats.Re) == Re & abs(double(binnedStats.kD) - kD) <= eps(max(1, abs(kD)));
T = binnedStats(mask, :);
if height(T) > 0
    [~, ord] = sort(double(T.gammaBinCenter));
    T = T(ord, :);
end
end


% =========================================================================
function [nRows, nCols] = subplot_grid(n)
if n <= 1
    nRows = 1; nCols = 1;
elseif n == 2
    nRows = 1; nCols = 2;
else
    nCols = ceil(sqrt(n));
    nRows = ceil(n / nCols);
end
end


% =========================================================================
function [x, y] = empirical_cdf(vals)
vals = sort(vals(isfinite(vals(:))));
n = numel(vals);
if n == 0
    x = nan(0,1);
    y = nan(0,1);
else
    x = vals(:);
    y = (1:n).' ./ n;
end
end


% =========================================================================
function q = finite_prctile(x, p)
x = sort(x(isfinite(x(:))));
n = numel(x);
if n == 0
    q = NaN;
    return;
end
idx = 1 + (n - 1) * p / 100;
i0 = floor(idx);
i1 = ceil(idx);
if i0 == i1
    q = x(i0);
else
    q = x(i0) + (idx - i0) * (x(i1) - x(i0));
end
end


% =========================================================================
function m = finite_median(x)
x = x(isfinite(x(:)));
if isempty(x)
    m = NaN;
else
    m = median(x);
end
end


% =========================================================================
function cmap = scientific_line_colormap(n)
base = [ ...
    0.05 0.28 0.63
    0.83 0.22 0.13
    0.09 0.55 0.32
    0.67 0.35 0.09
    0.34 0.26 0.55
    0.00 0.55 0.62
    0.55 0.12 0.34
    0.20 0.20 0.20];
if n <= size(base, 1)
    cmap = base(1:n, :);
else
    cmap = lines(n);
end
end


% =========================================================================
function c = lighten_color(c, frac)
c = c + frac * (1 - c);
c = max(0, min(1, c));
end


% =========================================================================
function c = theme_text_color(theme)
if strcmpi(theme, 'poster')
    c = [1 1 1];
else
    c = [0 0 0];
end
end


% =========================================================================
function style_legend_for_theme(leg, theme)
if isempty(leg) || ~isgraphics(leg)
    return;
end
if strcmpi(theme, 'poster')
    leg.TextColor = [1 1 1];
    leg.Color = [0 0 0];
else
    leg.TextColor = [0 0 0];
    leg.Color = 'none';
end
end


% =========================================================================
function close_if_needed(f, plotOpts)
if ~isfield(plotOpts, 'keepFiguresOpen') || ~plotOpts.keepFiguresOpen
    close(f);
end
end
