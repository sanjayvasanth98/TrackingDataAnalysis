function plot_lagrangian_acceleration_analysis(allLagAccel, outDir, plotOpts, lagAccelOpts)
%PLOT_LAGRANGIAN_ACCELERATION_ANALYSIS  Publication plots for acceleration proxy analysis.

if nargin < 3 || isempty(plotOpts)
    plotOpts = struct();
end
if ~isfield(plotOpts, 'themes') || isempty(plotOpts.themes)
    plotOpts.themes = "normal";
end
if nargin < 4 || isempty(lagAccelOpts)
    lagAccelOpts = struct();
end
lagAccelOpts = apply_plot_defaults(lagAccelOpts);

if isempty(allLagAccel)
    warning('plot_lagrangian_acceleration_analysis: no cases. Skipping.');
    return;
end

if ~isfolder(outDir)
    mkdir(outDir);
end

if lagAccelOpts.makeAllFramePdfPlot
    plot_accel_pdf_by_re(allLagAccel, outDir, plotOpts, lagAccelOpts, "allframe");
end
if lagAccelOpts.makeTriggerWindowPdfPlot
    plot_accel_pdf_by_re(allLagAccel, outDir, plotOpts, lagAccelOpts, "trigger");
end
if isfield(lagAccelOpts, 'makeTriggerSurvivalPlot') && lagAccelOpts.makeTriggerSurvivalPlot
    plot_trigger_survival_by_re(allLagAccel, outDir, plotOpts, lagAccelOpts);
end
if isfield(lagAccelOpts, 'makeTriggerThresholdSummaryPlot') && lagAccelOpts.makeTriggerThresholdSummaryPlot
    plot_trigger_survival_threshold_summary(allLagAccel, outDir, plotOpts, lagAccelOpts);
end
if lagAccelOpts.makePeakTriggerGrowthPlot
    plot_peak_trigger_vs_growth(allLagAccel, outDir, plotOpts, lagAccelOpts);
end
if lagAccelOpts.makeHeatmapPlots
    plot_accel_heatmaps_by_re(allLagAccel, outDir, plotOpts, lagAccelOpts);
end
if lagAccelOpts.makeSanityCheckPlots
    plot_accel_sanity_checks(allLagAccel, fullfile(outDir, "SanityChecks"), plotOpts, lagAccelOpts);
end
fprintf('Saved Lagrangian acceleration plots to: %s\n', outDir);
end


% =========================================================================
function opts = apply_plot_defaults(opts)
opts = default_field(opts, 'throatHeight_mm', 10);
opts = default_field(opts, 'xLimNorm', [0 0.5]);
opts = default_field(opts, 'yLimNorm', [0 0.12]);
opts = default_field(opts, 'pdfGridN', 260);
opts = default_field(opts, 'pdfPercentileRange', [0.5 99.5]);
opts = default_field(opts, 'allFramePdfLegendLocation', 'southoutside');
opts = default_field(opts, 'triggerPdfLegendLocation', 'southoutside');
opts = default_field(opts, 'makeTriggerPdfZoomPlot', false);
opts = default_field(opts, 'triggerPdfZoomXLim', [1e-1 1]);
opts = default_field(opts, 'makeTriggerSurvivalPlot', false);
opts = default_field(opts, 'triggerSurvivalXLim', [1e-1 1]);
opts = default_field(opts, 'triggerSurvivalLegendLocation', 'southwest');
opts = default_field(opts, 'showTriggerSurvivalLegend', true);
opts = default_field(opts, 'makeTriggerThresholdSummaryPlot', false);
opts = default_field(opts, 'triggerSurvivalThresholds', [0.3 0.5 0.75]);
opts = default_field(opts, 'peakGrowthLegendLocation', 'northwest');
opts = default_field(opts, 'peakGrowthLegendFontSize', 9);
opts = default_field(opts, 'makeAllFramePdfPlot', true);
opts = default_field(opts, 'makeTriggerWindowPdfPlot', true);
opts = default_field(opts, 'makePeakTriggerGrowthPlot', true);
opts = default_field(opts, 'makeHeatmapPlots', true);
opts = default_field(opts, 'makeSanityCheckPlots', true);
opts = default_field(opts, 'heatmapGridSize', [20 20]);
opts = default_field(opts, 'heatmapStats', ["median", "p90"]);
opts = default_field(opts, 'heatmapMinSamplesPerBin', 10);
opts = default_field(opts, 'heatmapColormap', "cbrewer2:YlGnBu");
opts = default_field(opts, 'heatmapColormapByStat', default_heatmap_colormap_by_stat());
opts = default_field(opts, 'heatmapPreserveSpatialAspect', true);
opts = default_field(opts, 'heatmapSquareBins', true);
opts = default_field(opts, 'heatmapShowActivationContours', true);
opts = default_field(opts, 'heatmapShowAccelerationRidge', false);
opts = default_field(opts, 'heatmapAccelerationRidgeStats', ["median", "mean", "p90"]);
opts = default_field(opts, 'heatmapAccelerationRidgePercentile', 0);
opts = default_field(opts, 'heatmapAccelerationRidgeSmoothWindow', 9);
opts = default_field(opts, 'heatmapColorPercentileRange', []);
opts = default_field(opts, 'heatmapActivationMarkerSize', 46);
opts = default_field(opts, 'heatmapLegendLocation', 'southoutside');
opts = default_field(opts, 'heatmapLegendFontSize', []);
opts = default_field(opts, 'activationOverlayPerCase', 75);
opts = default_field(opts, 'activationOverlayEdgeMarginNorm', 0);
opts = default_field(opts, 'randomSeed', 42);
opts = default_field(opts, 'minSpearmanN', 10);
opts = default_field(opts, 'maxSanityTracks', 250);
opts = default_field(opts, 'maxSanityTracksToPlot', opts.maxSanityTracks);
if ~(isfinite(opts.maxSanityTracksToPlot) && opts.maxSanityTracksToPlot >= 0)
    opts.maxSanityTracksToPlot = 250;
end
opts.maxSanityTracksToPlot = max(0, round(opts.maxSanityTracksToPlot));
end


% =========================================================================
function s = default_field(s, name, value)
if ~isfield(s, name) || isempty(s.(name))
    s.(name) = value;
end
end


% =========================================================================
function cmapByStat = default_heatmap_colormap_by_stat()
cmapByStat = struct();
cmapByStat.mean = "cbrewer2:GnBu";
cmapByStat.median = "cbrewer2:PuBu";
cmapByStat.p90 = "cbrewer2:Blues";
end


% =========================================================================
function plot_accel_pdf_by_re(allLagAccel, outDir, plotOpts, opts, mode)
ReVals = unique([allLagAccel.Re]);
ReVals = ReVals(isfinite(ReVals));
if isempty(ReVals), return; end

for theme = reshape(plotOpts.themes, 1, [])
    for ri = 1:numel(ReVals)
        Rei = ReVals(ri);
        caseIdx = find([allLagAccel.Re] == Rei);
        [xGrid, hasGrid] = common_pdf_grid(allLagAccel(caseIdx), mode, opts);
        if ~hasGrid
            continue;
        end

        fontName = resolve_plot_font_name();
        f = figure('Color', 'w', 'Position', [100 100 1120 760]);
        ax = axes(f);
        hold(ax, 'on');

        cmap = scientific_line_colormap(numel(caseIdx));
        lgd = gobjects(0, 1);
        lgdTxt = strings(0, 1);

        for j = 1:numel(caseIdx)
            d = allLagAccel(caseIdx(j));
            if mode == "trigger"
                yAct = positive_values(d.activatedTriggerAstar);
                yNon = positive_values(d.nonActivatedRandomAstar);
            else
                yAct = positive_values(d.activatedAllAstar);
                yNon = positive_values(d.nonActivatedAllAstar);
            end

            col = cmap(j, :);
            if numel(yAct) >= 5
                pdfAct = log_kde_pdf(yAct, xGrid);
                h = plot(ax, xGrid, pdfAct, '-', 'Color', col, 'LineWidth', 2.4);
                lgd(end+1,1) = h; %#ok<AGROW>
                lgdTxt(end+1,1) = sprintf('k/d=%.4g activated', d.kD); %#ok<AGROW>
            end
            if numel(yNon) >= 5
                pdfNon = log_kde_pdf(yNon, xGrid);
                h = plot(ax, xGrid, pdfNon, '--', 'Color', col, 'LineWidth', 2.1);
                lgd(end+1,1) = h; %#ok<AGROW>
                lgdTxt(end+1,1) = sprintf('k/d=%.4g non-activated', d.kD); %#ok<AGROW>
            end
        end

        set(ax, 'XScale', 'log', 'YScale', 'log', 'FontName', fontName);
        xlabel(ax, '$|a^*| = |\mathbf{a}| d_{\mathrm{mean}} / U_{\mathrm{ref}}^2$', 'Interpreter', 'latex');
        ylabel(ax, 'Probability density', 'Interpreter', 'latex');
        if mode == "trigger"
            title(ax, sprintf('Trigger-window acceleration PDFs, Re = %g', Rei), 'FontName', fontName);
            outName = sprintf('LagrangianAccel_trigger_PDF_Re_%g_%s', Rei, char(theme));
            legendLocation = char(opts.triggerPdfLegendLocation);
        else
            title(ax, sprintf('All-frame acceleration PDFs, Re = %g', Rei), 'FontName', fontName);
            outName = sprintf('LagrangianAccel_allframe_PDF_Re_%g_%s', Rei, char(theme));
            legendLocation = char(opts.allFramePdfLegendLocation);
        end
        grid(ax, 'off');
        box(ax, 'on');

        if ~isempty(lgd)
            leg = legend(ax, lgd, cellstr(lgdTxt), ...
                'Location', legendLocation, 'NumColumns', pdf_legend_num_columns(legendLocation, numel(lgdTxt)), 'Box', 'off');
        else
            leg = [];
        end
        apply_plot_theme(ax, char(theme));
        style_legend_for_theme(leg, char(theme));
        save_fig_dual_safe(f, fullfile(outDir, outName), plotOpts);
        if mode == "trigger"
            save_trigger_pdf_zoom_if_requested(f, ax, outDir, Rei, char(theme), plotOpts, opts);
        end
        close_if_needed(f, plotOpts);
    end
end
end


% =========================================================================
function save_trigger_pdf_zoom_if_requested(f, ax, outDir, Rei, theme, plotOpts, opts)
if ~isfield(opts, 'makeTriggerPdfZoomPlot') || ~opts.makeTriggerPdfZoomPlot
    return;
end
if ~isfield(opts, 'triggerPdfZoomXLim') || numel(opts.triggerPdfZoomXLim) < 2
    return;
end

zoomXLim = sort(double(opts.triggerPdfZoomXLim(1:2)));
if ~(all(isfinite(zoomXLim)) && zoomXLim(1) > 0 && zoomXLim(2) > zoomXLim(1))
    return;
end

origXLim = xlim(ax);
origYLim = ylim(ax);
origTitle = ax.Title.String;
origTitleInterp = ax.Title.Interpreter;
origTitleFont = ax.Title.FontName;
leg = legend(ax);
hasLegend = ~isempty(leg) && isgraphics(leg, 'legend');
if hasLegend
    origLegendVisible = leg.Visible;
    leg.Visible = 'off';
end

xlim(ax, zoomXLim);
set_pdf_zoom_ylim_from_visible_lines(ax, zoomXLim);
title(ax, sprintf('Trigger-window acceleration PDFs, Re = %g, %.2g \\leq |a^*| \\leq %.2g', ...
    Rei, zoomXLim(1), zoomXLim(2)), ...
    'Interpreter', 'tex', 'FontName', origTitleFont);

outName = sprintf('LagrangianAccel_trigger_PDF_zoom_0p1_to_1_Re_%g_%s', Rei, theme);
save_fig_dual_safe(f, fullfile(outDir, outName), plotOpts);

xlim(ax, origXLim);
ylim(ax, origYLim);
title(ax, origTitle, 'Interpreter', origTitleInterp, 'FontName', origTitleFont);
if hasLegend
    leg.Visible = origLegendVisible;
end
end


% =========================================================================
function set_pdf_zoom_ylim_from_visible_lines(ax, xLimUse)
lineHandles = findobj(ax, 'Type', 'Line');
yVals = nan(0, 1);
for i = 1:numel(lineHandles)
    x = get(lineHandles(i), 'XData');
    y = get(lineHandles(i), 'YData');
    x = x(:);
    y = y(:);
    n = min(numel(x), numel(y));
    x = x(1:n);
    y = y(1:n);
    inWindow = isfinite(x) & isfinite(y) & y > 0 & x >= xLimUse(1) & x <= xLimUse(2);
    yVals = [yVals; y(inWindow)]; %#ok<AGROW>
end

if isempty(yVals)
    return;
end
yLo = min(yVals);
yHi = max(yVals);
if ~(isfinite(yLo) && isfinite(yHi) && yHi > yLo)
    return;
end
ylim(ax, [max(yLo * 0.75, realmin), yHi * 1.35]);
end


% =========================================================================
function nCols = pdf_legend_num_columns(location, nItems)
location = lower(char(string(location)));
if contains(location, 'outside')
    nCols = 2;
else
    nCols = 1;
end
nCols = max(1, min(nCols, max(1, nItems)));
end


% =========================================================================
function plot_trigger_survival_by_re(allLagAccel, outDir, plotOpts, opts)
ReVals = unique([allLagAccel.Re]);
ReVals = ReVals(isfinite(ReVals));
if isempty(ReVals), return; end

xLimUse = sort(double(opts.triggerSurvivalXLim(1:2)));
if ~(all(isfinite(xLimUse)) && xLimUse(1) > 0 && xLimUse(2) > xLimUse(1))
    xLimUse = [1e-1 1];
end
xGrid = logspace(log10(xLimUse(1)), log10(xLimUse(2)), opts.pdfGridN).';

for theme = reshape(plotOpts.themes, 1, [])
    for ri = 1:numel(ReVals)
        Rei = ReVals(ri);
        caseIdx = find([allLagAccel.Re] == Rei);
        if isempty(caseIdx), continue; end

        fontName = resolve_plot_font_name();
        f = figure('Color', 'w', 'Position', [100 100 1120 760]);
        ax = axes(f);
        hold(ax, 'on');

        cmap = scientific_line_colormap(numel(caseIdx));
        lgd = gobjects(0, 1);
        lgdTxt = strings(0, 1);

        for j = 1:numel(caseIdx)
            d = allLagAccel(caseIdx(j));
            yAct = positive_values(d.activatedTriggerAstar);
            yNon = positive_values(d.nonActivatedRandomAstar);
            col = cmap(j, :);

            if numel(yAct) >= 5
                survAct = empirical_survival(yAct, xGrid);
                h = plot(ax, xGrid, survAct, '-', 'Color', col, 'LineWidth', 2.4);
                lgd(end+1,1) = h; %#ok<AGROW>
                lgdTxt(end+1,1) = sprintf('k/d=%.4g activated', d.kD); %#ok<AGROW>
            end
            if numel(yNon) >= 5
                survNon = empirical_survival(yNon, xGrid);
                h = plot(ax, xGrid, survNon, '--', 'Color', col, 'LineWidth', 2.1);
                lgd(end+1,1) = h; %#ok<AGROW>
                lgdTxt(end+1,1) = sprintf('k/d=%.4g non-activated', d.kD); %#ok<AGROW>
            end
        end

        set(ax, 'XScale', 'log', 'YScale', 'linear', 'FontName', fontName);
        xlim(ax, xLimUse);
        ylim(ax, [0 1]);
        xlabel(ax, 'Threshold $A$ in $|a^*| \geq A$', 'Interpreter', 'latex');
        ylabel(ax, '$P(|a^*| \geq A)$', 'Interpreter', 'latex');
        title(ax, sprintf('Trigger-window acceleration survival, Re = %g', Rei), 'FontName', fontName);
        grid(ax, 'off');
        box(ax, 'on');

        if opts.showTriggerSurvivalLegend && ~isempty(lgd)
            legendLocation = char(opts.triggerSurvivalLegendLocation);
            leg = legend(ax, lgd, cellstr(lgdTxt), ...
                'Location', legendLocation, ...
                'NumColumns', pdf_legend_num_columns(legendLocation, numel(lgdTxt)), ...
                'Box', 'off');
        else
            leg = [];
        end

        apply_plot_theme(ax, char(theme));
        style_legend_for_theme(leg, char(theme));
        outName = sprintf('LagrangianAccel_trigger_survival_Re_%g_%s', Rei, char(theme));
        save_fig_dual_safe(f, fullfile(outDir, outName), plotOpts);
        close_if_needed(f, plotOpts);
    end
end
end


% =========================================================================
function plot_trigger_survival_threshold_summary(allLagAccel, outDir, plotOpts, opts)
ReVals = unique([allLagAccel.Re]);
ReVals = ReVals(isfinite(ReVals));
if isempty(ReVals), return; end

thresholds = positive_values(opts.triggerSurvivalThresholds);
thresholds = unique(thresholds(:).', 'stable');
if isempty(thresholds)
    thresholds = [0.3 0.5 0.75];
end
nT = numel(thresholds);

for theme = reshape(plotOpts.themes, 1, [])
    for ri = 1:numel(ReVals)
        Rei = ReVals(ri);
        caseIdx = find([allLagAccel.Re] == Rei);
        if isempty(caseIdx), continue; end

        [~, sortIdx] = sort([allLagAccel(caseIdx).kD]);
        caseIdx = caseIdx(sortIdx);
        kDVals = [allLagAccel(caseIdx).kD];
        validKD = isfinite(kDVals) & kDVals > 0;
        if ~any(validKD), continue; end
        caseIdx = caseIdx(validKD);
        kDVals = kDVals(validKD);

        actP = nan(numel(caseIdx), nT);
        nonP = nan(numel(caseIdx), nT);
        for j = 1:numel(caseIdx)
            d = allLagAccel(caseIdx(j));
            yAct = positive_values(d.activatedTriggerAstar);
            yNon = positive_values(d.nonActivatedRandomAstar);
            for ti = 1:nT
                if numel(yAct) >= 5
                    actP(j, ti) = mean(yAct >= thresholds(ti));
                end
                if numel(yNon) >= 5
                    nonP(j, ti) = mean(yNon >= thresholds(ti));
                end
            end
        end

        fontName = resolve_plot_font_name();
        f = figure('Color', 'w', 'Position', [100 100 max(960, 320*nT) 440]);
        tl = tiledlayout(f, 1, nT, 'TileSpacing', 'compact', 'Padding', 'compact');
        actColor = [0.05 0.28 0.63];
        nonColor = [0.35 0.35 0.35];
        markerSize = 4.8;
        xPlot = 1:numel(kDVals);
        yUpper = max([actP(:); nonP(:)], [], 'omitnan');
        if ~(isfinite(yUpper) && yUpper > 0)
            yUpper = 1;
        end
        yUpper = min(1, max(0.15, ceil((yUpper + 0.04) * 10) / 10));

        for ti = 1:nT
            ax = nexttile(tl);
            hold(ax, 'on');
            hAct = plot(ax, xPlot, actP(:, ti), '-o', ...
                'Color', actColor, ...
                'MarkerFaceColor', actColor, ...
                'MarkerEdgeColor', [0.05 0.05 0.05], ...
                'LineWidth', 1.8, ...
                'MarkerSize', markerSize);
            hNon = plot(ax, xPlot, nonP(:, ti), '--s', ...
                'Color', nonColor, ...
                'MarkerFaceColor', 'w', ...
                'MarkerEdgeColor', nonColor, ...
                'LineWidth', 1.6, ...
                'MarkerSize', markerSize);

            set(ax, 'FontName', fontName);
            xlim(ax, [0.6 numel(kDVals)+0.4]);
            ylim(ax, [0 yUpper]);
            xticks(ax, xPlot);
            xticklabels(ax, compose('%.3g', kDVals));
            xtickangle(ax, 0);
            xlabel(ax, 'k/d', 'Interpreter', 'tex', 'FontSize', 12);
            if ti == 1
                ylabel(ax, 'P(|a*| >= A)', 'Interpreter', 'tex', 'FontSize', 12);
            end
            title(ax, sprintf('A = %.2g', thresholds(ti)), 'FontName', fontName, 'FontSize', 14);
            grid(ax, 'off');
            box(ax, 'on');
            apply_plot_theme(ax, char(theme));
            set(ax, 'FontSize', 10, 'LineWidth', 1.0);
            ax.XLabel.FontSize = 12;
            if isgraphics(ax.YLabel), ax.YLabel.FontSize = 12; end
            if isgraphics(ax.Title), ax.Title.FontSize = 14; end

            if ti == 1
                leg = legend(ax, [hAct hNon], {'activated', 'non-activated'}, ...
                    'Location', 'northoutside', ...
                    'Orientation', 'horizontal', ...
                    'Box', 'off');
                style_legend_for_theme(leg, char(theme));
                leg.FontSize = 10;
            end
        end

        title(tl, sprintf('Trigger-window survival probability, Re = %g', Rei), ...
            'FontName', fontName, ...
            'FontSize', 16);
        outName = sprintf('LagrangianAccel_trigger_survival_threshold_summary_Re_%g_%s', Rei, char(theme));
        save_fig_dual_safe(f, fullfile(outDir, outName), plotOpts);
        close_if_needed(f, plotOpts);
    end
end
end


% =========================================================================
function p = empirical_survival(values, xGrid)
values = positive_values(values);
xGrid = xGrid(:);
p = nan(size(xGrid));
if isempty(values)
    return;
end
for i = 1:numel(xGrid)
    p(i) = mean(values >= xGrid(i));
end
end


% =========================================================================
function [xGrid, ok] = common_pdf_grid(caseData, mode, opts)
vals = nan(0,1);
for i = 1:numel(caseData)
    if mode == "trigger"
        vals = [vals; positive_values(caseData(i).activatedTriggerAstar); ...
            positive_values(caseData(i).nonActivatedRandomAstar)]; %#ok<AGROW>
    else
        vals = [vals; positive_values(caseData(i).activatedAllAstar); ...
            positive_values(caseData(i).nonActivatedAllAstar)]; %#ok<AGROW>
    end
end
vals = positive_values(vals);
ok = numel(vals) >= 5;
if ~ok
    xGrid = nan(0,1);
    return;
end
pRange = opts.pdfPercentileRange;
xMin = finite_prctile(vals, pRange(1));
xMax = finite_prctile(vals, pRange(2));
if ~(isfinite(xMin) && xMin > 0), xMin = min(vals); end
if ~(isfinite(xMax) && xMax > xMin), xMax = max(vals); end
if ~(isfinite(xMax) && xMax > xMin)
    xMax = xMin * 10;
end
xMin = max(xMin, realmin);
xGrid = logspace(log10(xMin), log10(xMax), opts.pdfGridN).';
end


% =========================================================================
function plot_peak_trigger_vs_growth(allLagAccel, outDir, plotOpts, opts)
ReVals = unique([allLagAccel.Re]);
ReVals = ReVals(isfinite(ReVals));
if isempty(ReVals), return; end

for theme = reshape(plotOpts.themes, 1, [])
    fontName = resolve_plot_font_name();
    nRe = numel(ReVals);
    f = figure('Color', 'w', 'Position', [80 100 max(920, 520*nRe) 720]);

    for ri = 1:nRe
        ax = subplot(1, nRe, ri, 'Parent', f);
        hold(ax, 'on');
        Rei = ReVals(ri);
        caseIdx = find([allLagAccel.Re] == Rei);
        cmap = scientific_line_colormap(numel(caseIdx));
        insideLegendText = strings(0, 1);
        insideLegendColors = zeros(0, 3);
        lgd = gobjects(0, 1);
        lgdTxt = strings(0, 1);

        for j = 1:numel(caseIdx)
            d = allLagAccel(caseIdx(j));
            x = d.peakTriggerAstar(:);
            y = d.growthRatio(:);
            n = min(numel(x), numel(y));
            x = x(1:n);
            y = y(1:n);
            valid = isfinite(x) & x > 0 & isfinite(y) & y > 0;
            x = x(valid);
            y = y(valid);
            if isempty(x)
                continue;
            end

            col = cmap(j, :);
            h = scatter(ax, x, y, 52, ...
                'Marker', 'o', ...
                'MarkerFaceColor', col, ...
                'MarkerEdgeColor', [0.05 0.05 0.05], ...
                'LineWidth', 0.75, ...
                'MarkerFaceAlpha', 0.72);
            lgd(end+1,1) = h; %#ok<AGROW>
            lgdTxt(end+1,1) = sprintf('k/d=%.4g', d.kD); %#ok<AGROW>

            if numel(x) >= opts.minSpearmanN
                rho = spearman_rho(x, y);
                insideLegendText(end+1,1) = sprintf('k/d=%.4g: \\rho=%.2f (n=%d)', d.kD, rho, numel(x)); %#ok<AGROW>
            else
                insideLegendText(end+1,1) = sprintf('k/d=%.4g: n=%d', d.kD, numel(x)); %#ok<AGROW>
            end
            insideLegendColors(end+1,:) = col; %#ok<AGROW>
        end

        set(ax, 'XScale', 'log', 'FontName', fontName);
        xlabel(ax, 'Peak trigger-window $|a^*|$', 'Interpreter', 'latex');
        ylabel(ax, '$d_{\max}/d_{\mathrm{birth}}$', 'Interpreter', 'latex');
        title(ax, sprintf('Re = %g', Rei), 'FontName', fontName);
        grid(ax, 'off');
        box(ax, 'on');
        apply_plot_theme(ax, char(theme));
        draw_peak_growth_inside_legend(ax, insideLegendText, insideLegendColors, opts, fontName, char(theme));
    end
    save_fig_dual_safe(f, fullfile(outDir, "LagrangianAccel_peakTrigger_vs_growth_" + theme), plotOpts);
    close_if_needed(f, plotOpts);
end
end


% =========================================================================
function draw_peak_growth_inside_legend(ax, labels, colors, opts, fontName, theme)
if isempty(labels)
    return;
end

xLim = xlim(ax);
yLim = ylim(ax);
if numel(xLim) < 2 || numel(yLim) < 2
    return;
end
xLim = double(xLim);
yLim = double(yLim);
xSpanLog = diff(log10(xLim));
ySpan = diff(yLim);
if ~(all(isfinite(xLim)) && all(isfinite(yLim)) && xLim(1) > 0 && xSpanLog > 0 && ySpan > 0)
    return;
end

loc = 'northwest';
if isfield(opts, 'peakGrowthLegendLocation') && ~isempty(opts.peakGrowthLegendLocation)
    loc = lower(char(string(opts.peakGrowthLegendLocation)));
end

fontSize = 9;
if isfield(opts, 'peakGrowthLegendFontSize') && ~isempty(opts.peakGrowthLegendFontSize)
    fs = double(opts.peakGrowthLegendFontSize);
    if isfinite(fs) && fs > 0
        fontSize = fs;
    end
end

isWest = any(strcmp(loc, {'northwest', 'southwest', 'west'}));
isSouth = any(strcmp(loc, {'southwest', 'southeast', 'south'}));
if isWest
    xMarker = 10^(log10(xLim(1)) + 0.060 * xSpanLog);
    xText = 10^(log10(xLim(1)) + 0.105 * xSpanLog);
    hAlign = 'left';
else
    xMarker = 10^(log10(xLim(2)) - 0.22 * xSpanLog);
    xText = 10^(log10(xLim(2)) - 0.055 * xSpanLog);
    hAlign = 'right';
end

yStep = 0.040 * ySpan;
if isSouth
    yStart = yLim(1) + 0.10 * ySpan;
    yForItem = @(i) yStart + (i - 1) * yStep;
else
    yStart = yLim(2) - 0.10 * ySpan;
    yForItem = @(i) yStart - (i - 1) * yStep;
end

for i = 1:numel(labels)
    y = yForItem(i);
    if y < yLim(1) || y > yLim(2)
        break;
    end
    scatter(ax, xMarker, y, 44, ...
        'Marker', 'o', ...
        'MarkerFaceColor', colors(i,:), ...
        'MarkerEdgeColor', [0.05 0.05 0.05], ...
        'LineWidth', 0.75, ...
        'MarkerFaceAlpha', 0.72, ...
        'HandleVisibility', 'off');
    text(ax, xText, y, char(labels(i)), ...
        'HorizontalAlignment', hAlign, ...
        'VerticalAlignment', 'middle', ...
        'Interpreter', 'tex', ...
        'FontName', fontName, ...
        'FontSize', fontSize, ...
        'Color', theme_text_color(theme), ...
        'BackgroundColor', theme_background_color(theme), ...
        'Margin', 2, ...
        'Clipping', 'on', ...
        'HandleVisibility', 'off');
end
end


% =========================================================================
function plot_accel_heatmaps_by_re(allLagAccel, outDir, plotOpts, opts)
ReVals = unique([allLagAccel.Re]);
ReVals = ReVals(isfinite(ReVals));
if isempty(ReVals), return; end

stats = string(opts.heatmapStats);
for theme = reshape(plotOpts.themes, 1, [])
    for ri = 1:numel(ReVals)
        Rei = ReVals(ri);
        caseIdx = find([allLagAccel.Re] == Rei);
        if isempty(caseIdx), continue; end

        x = nan(0,1); y = nan(0,1); a = nan(0,1);
        for j = 1:numel(caseIdx)
            d = allLagAccel(caseIdx(j));
            x = [x; d.sampleXNorm(:)]; %#ok<AGROW>
            y = [y; d.sampleYNorm(:)]; %#ok<AGROW>
            a = [a; d.sampleAstar(:)]; %#ok<AGROW>
        end
        valid = isfinite(x) & isfinite(y) & isfinite(a) & a >= 0 & ...
            x >= opts.xLimNorm(1) & x <= opts.xLimNorm(2) & ...
            y >= opts.yLimNorm(1) & y <= opts.yLimNorm(2);
        x = x(valid); y = y(valid); a = a(valid);
        if isempty(a)
            continue;
        end

        for si = 1:numel(stats)
            statName = lower(char(stats(si)));
            [Z, xCenters, yCenters] = binned_stat_2d(x, y, a, opts, statName);
            plot_one_heatmap(allLagAccel(caseIdx), Z, xCenters, yCenters, Rei, statName, outDir, plotOpts, opts, char(theme));
        end
    end
end
end


% =========================================================================
function plot_one_heatmap(caseData, Z, xCenters, yCenters, Rei, statName, outDir, plotOpts, opts, theme)
fontName = resolve_plot_font_name();
if ~strcmpi(fontName, 'Times New Roman')
    fontName = 'Times New Roman';
end
figSize = heatmap_figure_size(plotOpts);
f = figure('Color', 'w', 'Position', [120 120 figSize]);
set(f, 'DefaultAxesFontName', fontName, 'DefaultTextFontName', fontName);
ax = axes(f, 'Position', heatmap_axes_position());
hold(ax, 'on');

hImg = imagesc(ax, xCenters, yCenters, Z);
set(hImg, 'AlphaData', isfinite(Z));
set(ax, 'YDir', 'normal', 'Color', 'w');
axis(ax, 'tight');
colormap(ax, lagrangian_heat_colormap(256, opts, statName));
apply_heatmap_color_limits(ax, Z, opts);
cb = colorbar(ax);
cb.Label.String = sprintf('%s |a^*| per bin', upper(statName));
cb.Label.Interpreter = 'tex';
cb.FontName = fontName;

hasROI = isfield(plotOpts, 'roiData') && isstruct(plotOpts.roiData);
if hasROI
    yExtent_mm = heatmap_y_extent_mm(caseData, opts);
    draw_lagrangian_wall_patch(ax, plotOpts.roiData, yExtent_mm, opts.throatHeight_mm, ...
        opts.xLimNorm, opts.yLimNorm, true);
end

cmap = scientific_line_colormap(numel(caseData));
markers = {'o', 'd', 'p', 's', '^', 'v', 'h', '>', '<', '*'};
lgd = gobjects(0, 1);
lgdTxt = strings(0, 1);
lgdColors = zeros(0, 3);
lgdMarkers = {};
rng(opts.randomSeed + round(Rei), 'twister');

drawContours = isfield(opts, 'heatmapShowActivationContours') && opts.heatmapShowActivationContours;
edgeMargin = activation_overlay_edge_margin(opts);
contourOverlays = struct('Xg', {}, 'Yg', {}, 'D', {}, 'color', {});
for j = 1:numel(caseData)
    d = caseData(j);
    xy = d.activationXYNorm;
    xMinOverlay = opts.xLimNorm(1) + edgeMargin(1);
    xMaxOverlay = opts.xLimNorm(2) - edgeMargin(1);
    yMinOverlay = opts.yLimNorm(1) + edgeMargin(2);
    yMaxOverlay = opts.yLimNorm(2) - edgeMargin(2);
    valid = all(isfinite(xy), 2) & xy(:,1) >= xMinOverlay & xy(:,1) <= xMaxOverlay & ...
        xy(:,2) >= yMinOverlay & xy(:,2) <= yMaxOverlay;
    xy = xy(valid, :);
    if isempty(xy)
        continue;
    end
    if size(xy, 1) > opts.activationOverlayPerCase
        idx = sort(randperm(size(xy, 1), opts.activationOverlayPerCase));
        xyPlot = xy(idx, :);
    else
        xyPlot = xy;
    end
    marker = markers{mod(j - 1, numel(markers)) + 1};
    h = scatter(ax, xyPlot(:,1), xyPlot(:,2), opts.heatmapActivationMarkerSize, ...
        'Marker', marker, ...
        'MarkerFaceColor', cmap(j,:), ...
        'MarkerEdgeColor', [0.03 0.03 0.03], ...
        'LineWidth', 0.8, ...
        'MarkerFaceAlpha', 0.82);
    lgd(end+1,1) = h; %#ok<AGROW>
    lgdTxt(end+1,1) = sprintf('k/d=%.4g', d.kD); %#ok<AGROW>
    lgdColors(end+1,:) = cmap(j,:); %#ok<AGROW>
    lgdMarkers{end+1,1} = marker; %#ok<AGROW>

    if drawContours && size(xy, 1) >= 10
        [Xg, Yg] = meshgrid(linspace(opts.xLimNorm(1), opts.xLimNorm(2), 90), ...
            linspace(opts.yLimNorm(1), opts.yLimNorm(2), 90));
        D = kde2d_simple(xy(:,1), xy(:,2), Xg, Yg);
        contourOverlays(end+1).Xg = Xg; %#ok<AGROW>
        contourOverlays(end).Yg = Yg;
        contourOverlays(end).D = D;
        contourOverlays(end).color = cmap(j,:);
    end
end

if drawContours
    for ci = 1:numel(contourOverlays)
        contour(ax, contourOverlays(ci).Xg, contourOverlays(ci).Yg, contourOverlays(ci).D, 1, ...
            'LineColor', contourOverlays(ci).color, ...
            'LineWidth', 2.4, ...
            'HandleVisibility', 'off');
    end
end

draw_acceleration_ridge_overlay(ax, Z, xCenters, yCenters, statName, opts, fontName);

xlim(ax, opts.xLimNorm);
ylim(ax, opts.yLimNorm);
if opts.heatmapPreserveSpatialAspect
    pbaspect(ax, [diff(opts.xLimNorm) diff(opts.yLimNorm) 1]);
end
xlabel(ax, 'x/H', 'Interpreter', 'tex', 'FontName', fontName);
ylabel(ax, 'y/H', 'Interpreter', 'tex', 'FontName', fontName);
title(ax, sprintf('Lagrangian acceleration heatmap (%s), Re = %g', upper(statName), Rei), ...
    'Interpreter', 'tex', 'FontName', fontName);
set(ax, 'FontName', fontName, 'Layer', 'top');
grid(ax, 'off');
box(ax, 'on');

useInlineLegend = heatmap_uses_inline_legend(opts);
if ~isempty(lgd) && ~useInlineLegend
    leg = legend(ax, lgd, cellstr(lgdTxt), 'Location', char(opts.heatmapLegendLocation), ...
        'NumColumns', heatmap_legend_num_columns(opts, numel(lgdTxt)), 'Box', 'off');
else
    leg = [];
end

apply_plot_theme(ax, theme);
style_heatmap_white_theme(f, ax, cb);
style_legend_for_theme(leg, theme);
style_legend_for_heatmap_white(leg);
style_heatmap_legend_text(leg, opts, fontName);
draw_heatmap_inline_legend(ax, lgdTxt, lgdColors, lgdMarkers, opts, fontName, theme);
save_fig_dual_safe(f, fullfile(outDir, sprintf('LagrangianAccel_heatmap_%s_Re_%g_%s', statName, Rei, theme)), plotOpts);
close_if_needed(f, plotOpts);
end


% =========================================================================
function tf = heatmap_uses_inline_legend(opts)
tf = false;
if ~isfield(opts, 'heatmapLegendLocation') || isempty(opts.heatmapLegendLocation)
    return;
end
loc = lower(char(string(opts.heatmapLegendLocation)));
tf = any(strcmp(loc, {'northeast', 'northwest', 'upperright', 'upperleft'}));
end


% =========================================================================
function nCols = heatmap_legend_num_columns(opts, nItems)
nCols = min(4, nItems);
if isfield(opts, 'heatmapLegendLocation') && ~isempty(opts.heatmapLegendLocation)
    loc = lower(char(string(opts.heatmapLegendLocation)));
    if ~contains(loc, 'outside')
        nCols = 1;
    end
end
nCols = max(1, nCols);
end


% =========================================================================
function style_heatmap_legend_text(leg, opts, fontName)
if isempty(leg) || ~isgraphics(leg)
    return;
end
leg.FontName = fontName;
if isfield(opts, 'heatmapLegendFontSize') && ~isempty(opts.heatmapLegendFontSize)
    fs = double(opts.heatmapLegendFontSize);
    if isfinite(fs) && fs > 0
        leg.FontSize = fs;
    end
end
end


% =========================================================================
function draw_acceleration_ridge_overlay(ax, Z, xCenters, yCenters, statName, opts, fontName)
if ~isfield(opts, 'heatmapShowAccelerationRidge') || ~opts.heatmapShowAccelerationRidge
    return;
end
if isempty(Z) || isempty(xCenters) || isempty(yCenters)
    return;
end

ridgeStats = string(opts.heatmapAccelerationRidgeStats);
if ~any(strcmpi(string(statName), ridgeStats))
    return;
end

zFinite = Z(isfinite(Z));
if numel(zFinite) < 3
    return;
end

ridgePct = double(opts.heatmapAccelerationRidgePercentile);
if ~(isfinite(ridgePct) && ridgePct >= 0 && ridgePct <= 100)
    ridgePct = 0;
end
ridgeThreshold = finite_prctile(zFinite, ridgePct);

xLine = xCenters(:);
yLine = nan(size(xLine));
for ix = 1:numel(xCenters)
    col = Z(:, ix);
    valid = isfinite(col);
    if ~any(valid)
        continue;
    end
    validRows = find(valid);
    [maxVal, localIdx] = max(col(valid));
    if isfinite(maxVal) && (ridgePct <= 0 || maxVal >= ridgeThreshold)
        yLine(ix) = yCenters(validRows(localIdx));
    end
end

[xLine, yLine] = interpolate_ridge_full_span(xLine, yLine);
validLine = isfinite(xLine) & isfinite(yLine);
if nnz(validLine) < 3
    return;
end

smoothWindow = max(1, round(double(opts.heatmapAccelerationRidgeSmoothWindow)));
if smoothWindow > 1
    yLine(validLine) = moving_nan_mean(yLine(validLine), smoothWindow);
end

validLine = isfinite(xLine) & isfinite(yLine);
if nnz(validLine) < 3
    return;
end

xDense = linspace(min(xLine(validLine)), max(xLine(validLine)), max(160, 8 * nnz(validLine)));
yDense = interp1(xLine(validLine), yLine(validLine), xDense, 'pchip');
validDense = isfinite(xDense) & isfinite(yDense);
if nnz(validDense) < 3
    return;
end

plot(ax, xDense(validDense), yDense(validDense), '--', ...
    'Color', [0.03 0.03 0.03], ...
    'LineWidth', 2.0, ...
    'HandleVisibility', 'off');
end


% =========================================================================
function [xFull, yFull] = interpolate_ridge_full_span(xLine, yLine)
xFull = xLine(:);
yFull = yLine(:);
valid = isfinite(xFull) & isfinite(yFull);
if nnz(valid) < 3
    return;
end
firstIdx = find(valid, 1, 'first');
lastIdx = find(valid, 1, 'last');
spanIdx = firstIdx:lastIdx;
spanValid = valid(spanIdx);
if nnz(spanValid) < 3
    return;
end
yFull(spanIdx) = interp1(xFull(spanIdx(spanValid)), yFull(spanIdx(spanValid)), ...
    xFull(spanIdx), 'pchip');
end


% =========================================================================
function ySmooth = moving_nan_mean(y, windowSize)
y = y(:);
ySmooth = y;
halfWindow = floor(windowSize / 2);
for i = 1:numel(y)
    lo = max(1, i - halfWindow);
    hi = min(numel(y), i + halfWindow);
    vals = y(lo:hi);
    vals = vals(isfinite(vals));
    if ~isempty(vals)
        ySmooth(i) = mean(vals);
    end
end
end


% =========================================================================
function draw_heatmap_inline_legend(ax, labels, colors, markers, opts, fontName, theme)
if ~heatmap_uses_inline_legend(opts) || isempty(labels)
    return;
end

xLim = xlim(ax);
yLim = ylim(ax);
xSpan = diff(xLim);
ySpan = diff(yLim);
if ~(isfinite(xSpan) && isfinite(ySpan) && xSpan > 0 && ySpan > 0)
    return;
end

fontSize = 8;
if isfield(opts, 'heatmapLegendFontSize') && ~isempty(opts.heatmapLegendFontSize)
    fs = double(opts.heatmapLegendFontSize);
    if isfinite(fs) && fs > 0
        fontSize = fs;
    end
end

textColor = theme_text_color(theme);
loc = lower(char(string(opts.heatmapLegendLocation)));
if any(strcmp(loc, {'northwest', 'upperleft'}))
    xMarker = xLim(1) + 0.035 * xSpan;
    xText = xLim(1) + 0.070 * xSpan;
    textAlign = 'left';
else
    xText = xLim(2) - 0.035 * xSpan;
    xMarker = xLim(2) - 0.18 * xSpan;
    textAlign = 'right';
end
yTop = yLim(2) - 0.10 * ySpan;
yStep = 0.055 * ySpan;

for i = 1:numel(labels)
    y = yTop - (i - 1) * yStep;
    if y < yLim(1)
        break;
    end
    scatter(ax, xMarker, y, opts.heatmapActivationMarkerSize, ...
        'Marker', markers{i}, ...
        'MarkerFaceColor', colors(i,:), ...
        'MarkerEdgeColor', [0.03 0.03 0.03], ...
        'LineWidth', 0.8, ...
        'HandleVisibility', 'off');
    text(ax, xText, y, char(labels(i)), ...
        'HorizontalAlignment', textAlign, ...
        'VerticalAlignment', 'middle', ...
        'Interpreter', 'tex', ...
        'FontName', fontName, ...
        'FontSize', fontSize, ...
        'Color', textColor, ...
        'Clipping', 'on', ...
        'HandleVisibility', 'off');
end
end


% =========================================================================
function apply_heatmap_color_limits(ax, Z, opts)
if ~isfield(opts, 'heatmapColorPercentileRange') || isempty(opts.heatmapColorPercentileRange)
    return;
end

pRange = double(opts.heatmapColorPercentileRange(:).');
if numel(pRange) < 2
    return;
end
pRange = sort(pRange(1:2));
pRange(1) = max(0, min(100, pRange(1)));
pRange(2) = max(0, min(100, pRange(2)));
if pRange(2) <= pRange(1)
    return;
end

vals = Z(isfinite(Z));
if numel(vals) < 2
    return;
end

cMin = finite_prctile(vals, pRange(1));
cMax = finite_prctile(vals, pRange(2));
if ~(isfinite(cMin) && isfinite(cMax) && cMax > cMin)
    cMin = min(vals);
    cMax = max(vals);
end
if isfinite(cMin) && isfinite(cMax) && cMax > cMin
    caxis(ax, [cMin cMax]);
end
end


% =========================================================================
function edgeMargin = activation_overlay_edge_margin(opts)
edgeMargin = [0 0];
if ~isfield(opts, 'activationOverlayEdgeMarginNorm') || isempty(opts.activationOverlayEdgeMarginNorm)
    return;
end

rawMargin = double(opts.activationOverlayEdgeMarginNorm(:).');
if isempty(rawMargin)
    return;
elseif numel(rawMargin) == 1
    edgeMargin = [rawMargin rawMargin];
else
    edgeMargin = rawMargin(1:2);
end

edgeMargin(~isfinite(edgeMargin)) = 0;
edgeMargin = max(edgeMargin, 0);
xHalfSpan = max(diff(opts.xLimNorm) / 2, 0);
yHalfSpan = max(diff(opts.yLimNorm) / 2, 0);
edgeMargin = min(edgeMargin, [xHalfSpan yHalfSpan]);
end


% =========================================================================
function [Z, xCenters, yCenters] = binned_stat_2d(x, y, a, opts, statName)
[nx, ny] = heatmap_bin_counts(opts);
xEdges = linspace(opts.xLimNorm(1), opts.xLimNorm(2), nx + 1);
yEdges = linspace(opts.yLimNorm(1), opts.yLimNorm(2), ny + 1);
xCenters = (xEdges(1:end-1) + xEdges(2:end)) / 2;
yCenters = (yEdges(1:end-1) + yEdges(2:end)) / 2;
Z = nan(ny, nx);

for ix = 1:nx
    inX = x >= xEdges(ix) & x < xEdges(ix + 1);
    if ix == nx
        inX = x >= xEdges(ix) & x <= xEdges(ix + 1);
    end
    for iy = 1:ny
        inY = y >= yEdges(iy) & y < yEdges(iy + 1);
        if iy == ny
            inY = y >= yEdges(iy) & y <= yEdges(iy + 1);
        end
        vals = a(inX & inY);
        vals = vals(isfinite(vals));
        if numel(vals) < opts.heatmapMinSamplesPerBin
            continue;
        end
        switch lower(char(statName))
            case 'mean'
                Z(iy, ix) = mean(vals);
            case 'p90'
                Z(iy, ix) = finite_prctile(vals, 90);
            otherwise
                Z(iy, ix) = median(vals);
        end
    end
end
end


% =========================================================================
function [nx, ny] = heatmap_bin_counts(opts)
gridSize = opts.heatmapGridSize;
if numel(gridSize) < 2
    gridSize = [gridSize gridSize];
end

nx = max(1, round(gridSize(1)));
ny = max(1, round(gridSize(2)));

if ~isfield(opts, 'heatmapSquareBins') || ~opts.heatmapSquareBins
    return;
end

xSpan = diff(opts.xLimNorm);
ySpan = diff(opts.yLimNorm);
if ~(isfinite(xSpan) && isfinite(ySpan) && xSpan > 0 && ySpan > 0)
    return;
end

% Keep the requested vertical resolution and choose the horizontal count so
% each heatmap cell has the same size in x/H and y/H units.
targetBinSize = ySpan / ny;
nx = max(1, round(xSpan / targetBinSize));
end


% =========================================================================
function plot_accel_sanity_checks(allLagAccel, sanityDir, plotOpts, opts)
if ~isfolder(sanityDir)
    mkdir(sanityDir);
end
for i = 1:numel(allLagAccel)
    d = allLagAccel(i);
    if isfield(d, 'opts') && isstruct(d.opts) && isfield(d.opts, 'makeSanityPlots') && ~d.opts.makeSanityPlots
        continue;
    end
    plot_raw_vs_smoothed_accel_hist(d, sanityDir, plotOpts);
    plot_raw_vs_smoothed_tracks(d, sanityDir, plotOpts, opts);
    plot_stationary_noise_floor(d, sanityDir, plotOpts);
end
end


% =========================================================================
function plot_raw_vs_smoothed_accel_hist(d, sanityDir, plotOpts)
rawVals = positive_values(d.sampleRawAstar);
smoothVals = positive_values(d.sampleAstar);
if numel(rawVals) < 5 || numel(smoothVals) < 5
    return;
end
allVals = [rawVals; smoothVals];
xMin = max(finite_prctile(allVals, 0.5), realmin);
xMax = finite_prctile(allVals, 99.5);
if ~(isfinite(xMax) && xMax > xMin), return; end
edges = logspace(log10(xMin), log10(xMax), 70);
xCenters = sqrt(edges(1:end-1) .* edges(2:end));
rawCounts = histcounts(rawVals, edges, 'Normalization', 'pdf');
smoothCounts = histcounts(smoothVals, edges, 'Normalization', 'pdf');

for theme = reshape(plotOpts.themes, 1, [])
    fontName = resolve_plot_font_name();
    f = figure('Color', 'w', 'Position', [100 100 980 680]);
    ax = axes(f);
    hold(ax, 'on');
    plot(ax, xCenters, rawCounts, '-', 'Color', [0.74 0.18 0.12], 'LineWidth', 2.1);
    plot(ax, xCenters, smoothCounts, '-', 'Color', [0.05 0.33 0.62], 'LineWidth', 2.4);
    set(ax, 'XScale', 'log', 'YScale', 'log', 'FontName', fontName);
    xlabel(ax, '$|a^*|$', 'Interpreter', 'latex');
    ylabel(ax, 'PDF', 'Interpreter', 'latex');
    title(ax, sprintf('Raw vs smoothed acceleration, %s', char(d.caseName)), 'FontName', fontName);
    txt = sprintf('Raw p90 %.3g, p99 %.3g\nSmoothed p90 %.3g, p99 %.3g', ...
        finite_prctile(rawVals, 90), finite_prctile(rawVals, 99), ...
        finite_prctile(smoothVals, 90), finite_prctile(smoothVals, 99));
    text(ax, 0.03, 0.96, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'FontName', fontName, 'Color', theme_text_color(char(theme)), ...
        'BackgroundColor', theme_background_color(char(theme)), 'Margin', 5);
    leg = legend(ax, {'Raw derivative', 'Smoothed derivative'}, 'Location', 'southoutside', 'NumColumns', 2, 'Box', 'off');
    grid(ax, 'off'); box(ax, 'on');
    apply_plot_theme(ax, char(theme));
    style_legend_for_theme(leg, char(theme));
    save_fig_dual_safe(f, fullfile(sanityDir, sprintf('Sanity_raw_vs_smoothed_accel_%s_Re_%g_kD_%g_%s', ...
        sanitize_case_token(d.caseName), d.Re, d.kD, char(theme))), plotOpts);
    close_if_needed(f, plotOpts);
end
end


% =========================================================================
function plot_raw_vs_smoothed_tracks(d, sanityDir, plotOpts, opts)
if ~isfield(d, 'sanityTracks') || isempty(d.sanityTracks)
    return;
end
maxPlotTracks = max(0, round(opts.maxSanityTracksToPlot));
nPlotTracks = min(numel(d.sanityTracks), maxPlotTracks);
if nPlotTracks < 1
    return;
end
sanityTracks = d.sanityTracks(1:nPlotTracks);
for theme = reshape(plotOpts.themes, 1, [])
    fontName = resolve_plot_font_name();
    figSize = heatmap_figure_size(plotOpts);
    f = figure('Color', 'w', 'Position', [120 120 figSize]);
    ax = axes(f, 'Position', heatmap_axes_position());
    hold(ax, 'on');
    cmap = scientific_line_colormap(nPlotTracks);

    if isfield(plotOpts, 'roiData') && isstruct(plotOpts.roiData)
        yExtent_mm = heatmap_y_extent_mm(d, opts);
        draw_lagrangian_wall_patch(ax, plotOpts.roiData, yExtent_mm, opts.throatHeight_mm, ...
            opts.xLimNorm, opts.yLimNorm, true);
    end

    hRawLegend = plot(ax, nan, nan, 'o-', ...
        'Color', [0.56 0.56 0.56], ...
        'MarkerSize', 4, ...
        'LineWidth', 0.9);
    hSmoothLegend = plot(ax, nan, nan, '-', ...
        'Color', [0.05 0.05 0.05], ...
        'LineWidth', 2.5);
    hActLegend = scatter(ax, nan, nan, 80, ...
        'Marker', '*', ...
        'MarkerFaceColor', [0.05 0.05 0.05], ...
        'MarkerEdgeColor', [0 0 0]);

    for i = 1:nPlotTracks
        tr = sanityTracks(i);
        col = cmap(i,:);
        plot(ax, tr.xRawNorm, tr.yRawNorm, 'o-', ...
            'Color', lighten_color(col, 0.55), ...
            'MarkerSize', 3.5, ...
            'LineWidth', 0.8, ...
            'HandleVisibility', 'off');
        plot(ax, tr.xSmoothNorm, tr.ySmoothNorm, '-', ...
            'Color', col, ...
            'LineWidth', 2.4, ...
            'HandleVisibility', 'off');
        if all(isfinite(tr.activationXYNorm))
            scatter(ax, tr.activationXYNorm(1), tr.activationXYNorm(2), 70, ...
                'Marker', '*', 'MarkerFaceColor', col, 'MarkerEdgeColor', [0 0 0], ...
                'HandleVisibility', 'off');
        end
    end
    xlim(ax, opts.xLimNorm);
    ylim(ax, opts.yLimNorm);
    xlabel(ax, '$x/H$', 'Interpreter', 'latex');
    ylabel(ax, '$y/H$', 'Interpreter', 'latex');
    title(ax, sprintf('Raw vs smoothed trajectories, %s', char(d.caseName)), 'FontName', fontName);
    grid(ax, 'off'); box(ax, 'on');
    if opts.heatmapPreserveSpatialAspect
        pbaspect(ax, [diff(opts.xLimNorm) diff(opts.yLimNorm) 1]);
    end
    leg = legend(ax, [hRawLegend hSmoothLegend hActLegend], ...
        {'Raw position: pale circles', 'Smoothed position: thick line', 'Activation location: star'}, ...
        'Location', 'southoutside', ...
        'NumColumns', 3, ...
        'Box', 'off');
    apply_plot_theme(ax, char(theme));
    set(ax, 'Layer', 'top');
    style_legend_for_theme(leg, char(theme));
    save_fig_dual_safe(f, fullfile(sanityDir, sprintf('Sanity_raw_vs_smoothed_tracks_%s_Re_%g_kD_%g_%s', ...
        sanitize_case_token(d.caseName), d.Re, d.kD, char(theme))), plotOpts);
    close_if_needed(f, plotOpts);
end
end


% =========================================================================
function plot_stationary_noise_floor(d, sanityDir, plotOpts)
stationaryVals = positive_values(d.stationaryAstar);
allVals = positive_values(d.sampleAstar);
if numel(stationaryVals) < 5 || numel(allVals) < 5
    return;
end
plotVals = [stationaryVals; allVals];
xMin = max(finite_prctile(plotVals, 0.5), realmin);
xMax = finite_prctile(plotVals, 99.5);
if ~(isfinite(xMax) && xMax > xMin), return; end
edges = logspace(log10(xMin), log10(xMax), 70);
xCenters = sqrt(edges(1:end-1) .* edges(2:end));
stationaryCounts = histcounts(stationaryVals, edges, 'Normalization', 'pdf');
allCounts = histcounts(allVals, edges, 'Normalization', 'pdf');
p95 = finite_prctile(stationaryVals, 95);

for theme = reshape(plotOpts.themes, 1, [])
    fontName = resolve_plot_font_name();
    f = figure('Color', 'w', 'Position', [100 100 980 680]);
    ax = axes(f);
    hold(ax, 'on');
    plot(ax, xCenters, allCounts, '-', 'Color', [0.15 0.40 0.72], 'LineWidth', 2.1);
    plot(ax, xCenters, stationaryCounts, '-', 'Color', [0.85 0.40 0.07], 'LineWidth', 2.4);
    set(ax, 'XScale', 'log', 'YScale', 'log', 'FontName', fontName);
    if isfinite(p95)
        yNow = ylim(ax);
        plot(ax, [p95 p95], yNow, ':', 'Color', [0.05 0.05 0.05], ...
            'LineWidth', 1.6, 'HandleVisibility', 'off');
        ylim(ax, yNow);
    end
    xlabel(ax, '$|a^*|$', 'Interpreter', 'latex');
    ylabel(ax, 'PDF', 'Interpreter', 'latex');
    title(ax, sprintf('Stationary-track acceleration noise floor, %s', char(d.caseName)), 'FontName', fontName);
    txt = sprintf('Stationary tracks: %d\nStationary p95 |a*| = %.3g', d.nStationaryTracks, p95);
    text(ax, 0.03, 0.96, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'FontName', fontName, 'Color', theme_text_color(char(theme)), ...
        'BackgroundColor', theme_background_color(char(theme)), 'Margin', 5);
    leg = legend(ax, {'All selected tracks', 'Stationary candidates'}, ...
        'Location', 'southoutside', 'NumColumns', 2, 'Box', 'off');
    grid(ax, 'off'); box(ax, 'on');
    apply_plot_theme(ax, char(theme));
    style_legend_for_theme(leg, char(theme));
    save_fig_dual_safe(f, fullfile(sanityDir, sprintf('Sanity_stationary_noise_floor_%s_Re_%g_kD_%g_%s', ...
        sanitize_case_token(d.caseName), d.Re, d.kD, char(theme))), plotOpts);
    close_if_needed(f, plotOpts);
end
end


% =========================================================================
function vals = positive_values(vals)
vals = vals(:);
vals = vals(isfinite(vals) & vals > 0);
end


% =========================================================================
function pdf = log_kde_pdf(values, xGrid)
values = positive_values(values);
xGrid = xGrid(:);
if isempty(values)
    pdf = nan(size(xGrid));
    return;
end
z = log10(values);
zGrid = log10(xGrid);
hz = silverman_bw(z);
densityZ = zeros(size(zGrid));
for i = 1:numel(z)
    dz = (zGrid - z(i)) / hz;
    densityZ = densityZ + exp(-0.5 * dz.^2);
end
densityZ = densityZ ./ (numel(z) * sqrt(2*pi) * hz);
pdf = densityZ ./ (xGrid * log(10));
end


% =========================================================================
function density = kde2d_simple(x, y, Xg, Yg)
x = x(:); y = y(:);
valid = isfinite(x) & isfinite(y);
x = x(valid); y = y(valid);
n = numel(x);
if n == 0
    density = nan(size(Xg));
    return;
end
hx = silverman_bw(x);
hy = silverman_bw(y);
density = zeros(size(Xg));
for i = 1:n
    dx = (Xg - x(i)) / hx;
    dy = (Yg - y(i)) / hy;
    density = density + exp(-0.5 * (dx.^2 + dy.^2));
end
density = density / (n * 2 * pi * hx * hy);
end


% =========================================================================
function h = silverman_bw(x)
x = x(isfinite(x));
n = numel(x);
if n < 2
    h = 1;
    return;
end
s = std(x, 0);
iqrVal = finite_prctile(x, 75) - finite_prctile(x, 25);
scale = min(s, iqrVal / 1.34);
if ~(isfinite(scale) && scale > 0)
    scale = max(s, iqrVal / 1.34);
end
if ~(isfinite(scale) && scale > 0)
    scale = 1;
end
h = max(0.9 * scale * n^(-1/5), eps);
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
function rho = spearman_rho(x, y)
x = x(:); y = y(:);
valid = isfinite(x) & isfinite(y);
x = x(valid); y = y(valid);
if numel(x) < 2
    rho = NaN;
    return;
end
rx = tied_ranks(x);
ry = tied_ranks(y);
rx = rx - mean(rx);
ry = ry - mean(ry);
den = sqrt(sum(rx.^2) * sum(ry.^2));
if den <= 0
    rho = NaN;
else
    rho = sum(rx .* ry) / den;
end
end


% =========================================================================
function r = tied_ranks(x)
[xs, ord] = sort(x(:));
rSorted = nan(size(xs));
i = 1;
while i <= numel(xs)
    j = i;
    while j < numel(xs) && xs(j+1) == xs(i)
        j = j + 1;
    end
    rSorted(i:j) = (i + j) / 2;
    i = j + 1;
end
r = nan(size(x(:)));
r(ord) = rSorted;
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
function cmap = scientific_heat_colormap(n)
anchors = [ ...
    0.015 0.045 0.130
    0.020 0.190 0.340
    0.000 0.420 0.470
    0.420 0.650 0.360
    0.920 0.710 0.260
    0.870 0.330 0.120
    0.560 0.060 0.080];
x = linspace(0, 1, size(anchors, 1));
xi = linspace(0, 1, n);
cmap = interp1(x, anchors, xi, 'linear');
cmap = max(0, min(1, cmap));
end


% =========================================================================
function cmap = lagrangian_heat_colormap(n, opts, statName)
mapName = "sky";
if isfield(opts, 'heatmapColormap') && ~isempty(opts.heatmapColormap)
    mapName = string(opts.heatmapColormap);
end
if nargin >= 3 && isfield(opts, 'heatmapColormapByStat') && isstruct(opts.heatmapColormapByStat)
    statField = matlab.lang.makeValidName(lower(char(string(statName))));
    if isfield(opts.heatmapColormapByStat, statField) && ~isempty(opts.heatmapColormapByStat.(statField))
        mapName = string(opts.heatmapColormapByStat.(statField));
    end
end

mapKey = lower(char(mapName));
if startsWith(mapKey, 'cbrewer2:') || startsWith(mapKey, 'brewer2:')
    parts = split(string(mapName), ':');
    cmap = colorbrewer2_seq_colormap(char(parts(end)), n);
    return;
end

switch mapKey
    case 'abyss'
        cmap = abyss_colormap_compat(n);
        return;
    case 'sky'
        cmap = sky_colormap_compat(n);
        return;
    case 'ylorrd'
        cmap = colorbrewer2_seq_colormap('YlOrRd', n);
        return;
    case 'ylgnbu'
        cmap = colorbrewer2_seq_colormap('YlGnBu', n);
        return;
    case 'orrd'
        cmap = colorbrewer2_seq_colormap('OrRd', n);
        return;
    case 'gnbu'
        cmap = colorbrewer2_seq_colormap('GnBu', n);
        return;
    case 'pubu'
        cmap = colorbrewer2_seq_colormap('PuBu', n);
        return;
    case 'blues'
        cmap = colorbrewer2_seq_colormap('Blues', n);
        return;
end

try
    cmap = feval(char(mapName), n);
catch
    cmap = scientific_heat_colormap(n);
end
end


% =========================================================================
function cmap = colorbrewer2_seq_colormap(mapName, n)
if nargin < 2 || isempty(n)
    n = 256;
end

ensure_cbrewer2_on_path();
try
    cmap = cbrewer2('seq', char(mapName), n, 'pchip', 'rgb');
    cmap = normalize_colormap_rgb(cmap);
    return;
catch
end

try
    cmap = cbrewer('seq', char(mapName), n, 'pchip');
    cmap = normalize_colormap_rgb(cmap);
    return;
catch
end

switch lower(char(mapName))
    case 'gnbu'
        anchors = [ ...
            247 252 240
            224 243 219
            204 235 197
            168 221 181
            123 204 196
             78 179 211
             43 140 190
              8 104 172
              8  64 129] / 255;
    case 'pubu'
        anchors = [ ...
            255 247 251
            236 231 242
            208 209 230
            166 189 219
            116 169 207
             54 144 192
              5 112 176
              4  90 141
              2  56  88] / 255;
    case 'blues'
        anchors = [ ...
            247 251 255
            222 235 247
            198 219 239
            158 202 225
            107 174 214
             66 146 198
             33 113 181
              8  81 156
              8  48 107] / 255;
    case 'ylgnbu'
        anchors = [ ...
            255 255 217
            237 248 177
            199 233 180
            127 205 187
             65 182 196
             29 145 192
             34  94 168
             37  52 148
              8  29  88] / 255;
    case 'orrd'
        anchors = [ ...
            255 247 236
            254 232 200
            253 212 158
            253 187 132
            252 141  89
            239 101  72
            215  48  31
            179   0   0
            127   0   0] / 255;
    otherwise % YlOrRd
        anchors = [ ...
            255 255 204
            255 237 160
            254 217 118
            254 178  76
            253 141  60
            252  78  42
            227  26  28
            189   0  38
            128   0  38] / 255;
end

x = linspace(0, 1, size(anchors, 1));
xi = linspace(0, 1, n);
cmap = interp1(x, anchors, xi, 'pchip');
cmap = normalize_colormap_rgb(cmap);
end


% =========================================================================
function ensure_cbrewer2_on_path()
if exist('cbrewer2', 'file') == 2
    return;
end

thisDir = fileparts(mfilename('fullpath'));
cbrewer2Dir = fullfile(thisDir, 'external', 'cbrewer2', 'cbrewer2');
if isfolder(cbrewer2Dir)
    addpath(cbrewer2Dir);
end
end


% =========================================================================
function cmap = normalize_colormap_rgb(cmap)
cmap = double(cmap);
if isempty(cmap) || size(cmap, 2) ~= 3
    cmap = scientific_heat_colormap(256);
    return;
end
if max(cmap(:), [], 'omitnan') > 1
    cmap = cmap / 255;
end
cmap = max(0, min(1, cmap));
end


% =========================================================================
function cmap = sky_colormap_compat(n)
% Use MATLAB's sky map when available; otherwise use a readable blue ramp.
if nargin < 1 || isempty(n)
    n = 256;
end

try
    cmap = sky(n);
    return;
catch
end

anchors = [ ...
    0.925 0.975 1.000
    0.760 0.910 1.000
    0.500 0.760 0.960
    0.230 0.560 0.850
    0.075 0.330 0.680
    0.020 0.120 0.360];
x = linspace(0, 1, size(anchors, 1));
xi = linspace(0, 1, n);
cmap = interp1(x, anchors, xi, 'pchip');
cmap = max(0, min(1, cmap));
end


% =========================================================================
function cmap = abyss_colormap_compat(n)
% Use MATLAB's abyss map when available; otherwise use a close dark-blue ramp.
if nargin < 1 || isempty(n)
    n = 256;
end

try
    cmap = abyss(n);
    return;
catch
end

anchors = [ ...
    0.000 0.000 0.020
    0.000 0.020 0.080
    0.010 0.070 0.180
    0.020 0.150 0.320
    0.030 0.260 0.500
    0.050 0.420 0.660
    0.170 0.610 0.780
    0.500 0.820 0.900];
x = linspace(0, 1, size(anchors, 1));
xi = linspace(0, 1, n);
cmap = interp1(x, anchors, xi, 'pchip');
cmap = max(0, min(1, cmap));
end


% =========================================================================
function figSize = heatmap_figure_size(plotOpts)
figSize = [1280 512];
if isfield(plotOpts, 'lagrangianHeatmapImageSize_px') && numel(plotOpts.lagrangianHeatmapImageSize_px) >= 2
    figSize = double(plotOpts.lagrangianHeatmapImageSize_px(1:2));
elseif isfield(plotOpts, 'inceptionImageSize_px') && numel(plotOpts.inceptionImageSize_px) >= 2
    figSize = double([plotOpts.inceptionImageSize_px(1), 512]);
end
figSize = reshape(figSize, 1, []);
if any(~isfinite(figSize)) || any(figSize <= 0)
    figSize = [1280 512];
end
end


% =========================================================================
function pos = heatmap_axes_position()
% Match the wide inception-location canvas while leaving room for colorbar/legend.
pos = [0.08 0.24 0.78 0.58];
end


% =========================================================================
function yExtent_mm = heatmap_y_extent_mm(caseData, opts)
pixelSizeVals = nan(0, 1);
for i = 1:numel(caseData)
    if isfield(caseData(i), 'pixelSize_mm') && isfinite(caseData(i).pixelSize_mm) && caseData(i).pixelSize_mm > 0
        pixelSizeVals(end+1, 1) = caseData(i).pixelSize_mm; %#ok<AGROW>
    end
end

if ~isempty(pixelSizeVals) && isfield(opts, 'imageSize_px') && numel(opts.imageSize_px) >= 2
    yExtent_mm = double(opts.imageSize_px(2)) * median(pixelSizeVals);
else
    yExtent_mm = opts.yLimNorm(2) * opts.throatHeight_mm;
end
end


% =========================================================================
function draw_lagrangian_wall_patch(ax, roiData, yExtent_mm, throatHeight_mm, xLimNorm, yLimNorm, doNormalize)
if nargin < 7
    doNormalize = true;
end
if ~isfield(roiData, 'wallMask') || ~isfield(roiData, 'maskPixelSize')
    return;
end

wallMask = roiData.wallMask;
ps = roiData.maskPixelSize;
wallCols = find(any(wallMask, 1));
if isempty(wallCols) || ~(isfinite(ps) && ps > 0) || ~(isfinite(throatHeight_mm) && throatHeight_mm > 0)
    return;
end

wallTopRow = zeros(size(wallCols));
for i = 1:numel(wallCols)
    wallTopRow(i) = find(wallMask(:, wallCols(i)), 1, 'first');
end

wallX_mm = wallCols(:) * ps;
wallYSurface_mm = yExtent_mm - wallTopRow(:) * ps;
if doNormalize
    wallX = wallX_mm / throatHeight_mm;
    wallYSurface = wallYSurface_mm / throatHeight_mm;
else
    wallX = wallX_mm;
    wallYSurface = wallYSurface_mm;
end

inRange = wallX >= xLimNorm(1) & wallX <= xLimNorm(2);
wallX = wallX(inRange);
wallYSurface = wallYSurface(inRange);
if isempty(wallX)
    return;
end

patchX = [wallX; flipud(wallX)];
patchY = [wallYSurface; repmat(yLimNorm(1), numel(wallX), 1)];
patch(ax, patchX, patchY, [0.80 0.80 0.80], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 1.0, ...
    'HandleVisibility', 'off');
end


% =========================================================================
function style_heatmap_white_theme(f, ax, cb)
set(f, 'Color', 'w', 'InvertHardcopy', 'off');
set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'Layer', 'top');
if isgraphics(ax.Title), ax.Title.Color = [0 0 0]; end
if isgraphics(ax.XLabel), ax.XLabel.Color = [0 0 0]; end
if isgraphics(ax.YLabel), ax.YLabel.Color = [0 0 0]; end
if ~isempty(cb) && isgraphics(cb)
    cb.Color = [0 0 0];
    if isprop(cb, 'Label') && isgraphics(cb.Label)
        cb.Label.Color = [0 0 0];
    end
end
end


% =========================================================================
function style_legend_for_heatmap_white(leg)
if isempty(leg) || ~isgraphics(leg)
    return;
end
leg.TextColor = [0 0 0];
leg.Color = 'none';
end


% =========================================================================
function c = lighten_color(c, frac)
c = c + frac * (1 - c);
c = max(0, min(1, c));
end


% =========================================================================
function bg = theme_background_color(theme)
if strcmpi(theme, 'poster')
    bg = [0 0 0];
else
    bg = [1 1 1];
end
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
function style_colorbar_for_theme(cb, theme)
if isempty(cb) || ~isgraphics(cb)
    return;
end
if strcmpi(theme, 'poster')
    cb.Color = [1 1 1];
    if isprop(cb, 'Label') && isgraphics(cb.Label)
        cb.Label.Color = [1 1 1];
    end
else
    cb.Color = [0 0 0];
    if isprop(cb, 'Label') && isgraphics(cb.Label)
        cb.Label.Color = [0 0 0];
    end
end
end


% =========================================================================
function close_if_needed(f, plotOpts)
if ~isfield(plotOpts, 'keepFiguresOpen') || ~plotOpts.keepFiguresOpen
    close(f);
end
end


% =========================================================================
function token = sanitize_case_token(caseName)
token = char(string(caseName));
token = regexprep(token, '[^A-Za-z0-9]+', '_');
token = regexprep(token, '^_+', '');
token = regexprep(token, '_+$', '');
if isempty(token)
    token = 'case';
end
end
