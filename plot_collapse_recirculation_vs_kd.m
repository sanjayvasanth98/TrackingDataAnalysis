function plot_collapse_recirculation_vs_kd(allCollapseRecirculation, outDir, plotOpts)
%PLOT_COLLAPSE_RECIRCULATION_VS_KD  Collapse-recirculation metrics vs k/d.
%
%   Saves:
%     CollapseGeneratedPerCollapse_vs_kD_<theme>
%     PctActivationFromCollapse_vs_kD_<theme>

if nargin < 3 || ~isfield(plotOpts, 'themes') || isempty(plotOpts.themes)
    plotOpts.themes = "normal";
end

T = collapse_recirculation_to_table(allCollapseRecirculation);
if isempty(T) || height(T) == 0
    warning('plot_collapse_recirculation_vs_kd: no cases. Skipping.');
    return;
end

if ~isfolder(outDir)
    mkdir(outDir);
end

plot_generated_per_collapse(T, outDir, plotOpts);
plot_pct_activation_from_collapse(T, outDir, plotOpts);
fprintf('Saved collapse recirculation vs k/d plots to: %s\n', outDir);
end


% =========================================================================
function plot_generated_per_collapse(T, outDir, plotOpts)
kD = T.kD(:);
Re = T.Re(:);
nGenerated = double(T.nCollapseGeneratedMicrobubbles(:));
nCollapse = double(T.nCollapseEvents(:));

y = nGenerated ./ nCollapse;
[ciLow, ciHigh] = poisson_rate_ci(nGenerated, nCollapse, 0.05);
valid = isfinite(kD) & isfinite(Re) & isfinite(y) & nCollapse > 0 & nGenerated >= 0;

plot_metric_with_errorbars(kD(valid), Re(valid), y(valid), ciLow(valid), ciHigh(valid), outDir, plotOpts, ...
    'CollapseGeneratedPerCollapse_vs_kD', ...
    'Collapse-generated microbubbles per collapse', ...
    '$N_{\mathrm{generated}} / N_{\mathrm{collapse}}$');
end


% =========================================================================
function plot_pct_activation_from_collapse(T, outDir, plotOpts)
kD = T.kD(:);
Re = T.Re(:);
nFromCollapse = double(T.nCollapseGeneratedActivationEvents(:));
nActivation = double(T.nTotalActivationEvents(:));

y = 100 * nFromCollapse ./ nActivation;
[pLow, pHigh] = wilson_proportion_ci(nFromCollapse, nActivation, 0.05);
ciLow = 100 * pLow;
ciHigh = 100 * pHigh;
valid = isfinite(kD) & isfinite(Re) & isfinite(y) & ...
    nActivation > 0 & nFromCollapse >= 0 & nFromCollapse <= nActivation;

plot_metric_with_errorbars(kD(valid), Re(valid), y(valid), ciLow(valid), ciHigh(valid), outDir, plotOpts, ...
    'PctActivationFromCollapse_vs_kD', ...
    'Activation events attributed to collapse-generated microbubbles', ...
    'Collapse-attributed activation events (\%)');
end


% =========================================================================
function plot_metric_with_errorbars(kD, Re, y, ciLow, ciHigh, outDir, plotOpts, outStem, titleStr, yLabelStr)
valid = isfinite(kD) & isfinite(Re) & isfinite(y);
kD = kD(valid);
Re = Re(valid);
y = y(valid);
ciLow = ciLow(valid);
ciHigh = ciHigh(valid);

if isempty(kD)
    warning('plot_collapse_recirculation_vs_kd:%s: no valid rows. Skipping.', outStem);
    return;
end

ciLow(~isfinite(ciLow)) = y(~isfinite(ciLow));
ciHigh(~isfinite(ciHigh)) = y(~isfinite(ciHigh));
ciLow = min(ciLow, y);
ciHigh = max(ciHigh, y);
errLow = max(0, y - ciLow);
errHigh = max(0, ciHigh - y);

ReVals = unique(Re(:));
markerStyles = {'o', 's', '^', 'd', 'v', '>', '<', 'p', 'h'};

for theme = reshape(plotOpts.themes, 1, [])
    fontName = resolve_plot_font_name();
    f = figure('Color', 'w', 'Position', [100 100 1000 700]);
    ax = axes(f);
    hold(ax, 'on');

    cmap = lines(max(numel(ReVals), 1));
    lgd = gobjects(0, 1);
    lgdTxt = strings(0, 1);

    for r = 1:numel(ReVals)
        Rei = ReVals(r);
        mask = Re == Rei;
        if ~any(mask)
            continue;
        end

        x_r = kD(mask);
        y_r = y(mask);
        eLow_r = errLow(mask);
        eHigh_r = errHigh(mask);
        [x_r, si] = sort(x_r);
        y_r = y_r(si);
        eLow_r = eLow_r(si);
        eHigh_r = eHigh_r(si);

        col = cmap(r, :);
        markerStyle = markerStyles{mod(r - 1, numel(markerStyles)) + 1};

        if numel(x_r) >= 2
            plot(ax, x_r, y_r, '-', ...
                'Color', col, ...
                'LineWidth', 0.75, ...
                'HandleVisibility', 'off');
        end

        hErr = errorbar(ax, x_r, y_r, eLow_r, eHigh_r, markerStyle, ...
            'LineStyle', 'none', ...
            'LineWidth', 0.9, ...
            'MarkerSize', 12, ...
            'CapSize', 8, ...
            'Color', col, ...
            'MarkerFaceColor', col, ...
            'MarkerEdgeColor', [0 0 0]);
        lgd(end+1,1) = hErr; %#ok<AGROW>
        lgdTxt(end+1,1) = sprintf('Re=%g', Rei); %#ok<AGROW>
    end

    xlabel(ax, '$k/d$', 'Interpreter', 'latex');
    ylabel(ax, yLabelStr, 'Interpreter', 'latex');
    title(ax, titleStr, 'FontName', fontName, 'FontSize', 12);
    set(ax, 'YScale', 'linear', 'FontName', fontName);
    grid(ax, 'off');
    box(ax, 'on');

    yTop = max(y + errHigh, [], 'omitnan');
    if isfinite(yTop) && yTop > 0
        ylim(ax, [0, 1.15 * yTop]);
    end

    if ~isempty(lgd)
        leg = legend(ax, lgd, cellstr(lgdTxt), ...
            'Location', 'southoutside', 'NumColumns', 2, 'Box', 'off');
    else
        leg = [];
    end

    apply_plot_theme(ax, char(theme));
    style_legend_for_theme(leg, char(theme));

    outBase = fullfile(outDir, sprintf('%s_%s', outStem, char(theme)));
    save_fig_dual_safe(f, outBase, plotOpts);
    if ~isfield(plotOpts, 'keepFiguresOpen') || ~plotOpts.keepFiguresOpen
        close(f);
    end
end
end


% =========================================================================
function [rateLow, rateHigh] = poisson_rate_ci(counts, exposure, alpha)
counts = counts(:);
exposure = exposure(:);
rateLow = nan(size(counts));
rateHigh = nan(size(counts));

for i = 1:numel(counts)
    k = counts(i);
    n = exposure(i);
    if ~(isfinite(k) && isfinite(n) && k >= 0 && n > 0)
        continue;
    end

    try
        if k == 0
            countLow = 0;
        else
            countLow = gammaincinv(alpha / 2, k);
        end
        countHigh = gammaincinv(1 - alpha / 2, k + 1);
    catch
        % Fallback if gammaincinv is unavailable: normal approximation.
        countLow = max(0, k - 1.96 * sqrt(max(k, 1)));
        countHigh = k + 1.96 * sqrt(max(k, 1));
    end

    rateLow(i) = countLow / n;
    rateHigh(i) = countHigh / n;
end
end


% =========================================================================
function [pLow, pHigh] = wilson_proportion_ci(successes, totals, alpha)
successes = successes(:);
totals = totals(:);
pLow = nan(size(successes));
pHigh = nan(size(successes));
z = 1.95996398454005; % 95% normal quantile
if nargin >= 3 && isfinite(alpha) && alpha ~= 0.05
    z = sqrt(2) * erfinv(1 - alpha);
end

for i = 1:numel(successes)
    k = successes(i);
    n = totals(i);
    if ~(isfinite(k) && isfinite(n) && n > 0 && k >= 0 && k <= n)
        continue;
    end

    phat = k / n;
    denom = 1 + z^2 / n;
    center = (phat + z^2 / (2 * n)) / denom;
    halfWidth = z * sqrt((phat * (1 - phat) / n) + (z^2 / (4 * n^2))) / denom;
    pLow(i) = max(0, center - halfWidth);
    pHigh(i) = min(1, center + halfWidth);
end
end


% =========================================================================
function style_legend_for_theme(leg, theme)
if isempty(leg) || ~isgraphics(leg)
    return;
end

if strcmp(theme, 'poster')
    leg.TextColor = [1 1 1];
    leg.Color = [0 0 0];
else
    leg.TextColor = [0 0 0];
    leg.Color = 'none';
end
end
