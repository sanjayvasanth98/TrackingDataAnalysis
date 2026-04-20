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
opts = default_field(opts, 'makeAllFramePdfPlot', true);
opts = default_field(opts, 'makeTriggerWindowPdfPlot', true);
opts = default_field(opts, 'makePeakTriggerGrowthPlot', true);
opts = default_field(opts, 'makeHeatmapPlots', true);
opts = default_field(opts, 'makeSanityCheckPlots', true);
opts = default_field(opts, 'heatmapGridSize', [20 20]);
opts = default_field(opts, 'heatmapStats', ["median", "p90"]);
opts = default_field(opts, 'heatmapMinSamplesPerBin', 10);
opts = default_field(opts, 'heatmapColormap', "sky");
opts = default_field(opts, 'heatmapPreserveSpatialAspect', true);
opts = default_field(opts, 'heatmapSquareBins', true);
opts = default_field(opts, 'activationOverlayPerCase', 75);
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
        else
            title(ax, sprintf('All-frame acceleration PDFs, Re = %g', Rei), 'FontName', fontName);
            outName = sprintf('LagrangianAccel_allframe_PDF_Re_%g_%s', Rei, char(theme));
        end
        grid(ax, 'off');
        box(ax, 'on');

        if ~isempty(lgd)
            leg = legend(ax, lgd, cellstr(lgdTxt), ...
                'Location', 'southoutside', 'NumColumns', 2, 'Box', 'off');
        else
            leg = [];
        end
        apply_plot_theme(ax, char(theme));
        style_legend_for_theme(leg, char(theme));
        save_fig_dual_safe(f, fullfile(outDir, outName), plotOpts);
        close_if_needed(f, plotOpts);
    end
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
        rhoLines = strings(0, 1);
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
                rhoLines(end+1,1) = sprintf('k/d=%.4g: \\rho=%.2f (n=%d)', d.kD, rho, numel(x)); %#ok<AGROW>
            else
                rhoLines(end+1,1) = sprintf('k/d=%.4g: n=%d', d.kD, numel(x)); %#ok<AGROW>
            end
        end

        set(ax, 'XScale', 'log', 'FontName', fontName);
        xlabel(ax, 'Peak trigger-window $|a^*|$', 'Interpreter', 'latex');
        ylabel(ax, '$d_{\max}/d_{\mathrm{birth}}$', 'Interpreter', 'latex');
        title(ax, sprintf('Re = %g', Rei), 'FontName', fontName);
        grid(ax, 'off');
        box(ax, 'on');
        if ~isempty(rhoLines)
            text(ax, 0.03, 0.97, strjoin(cellstr(rhoLines), sprintf('\n')), ...
                'Units', 'normalized', 'VerticalAlignment', 'top', ...
                'FontName', fontName, 'FontSize', 9, ...
                'Color', theme_text_color(char(theme)), ...
                'BackgroundColor', theme_background_color(char(theme)), ...
                'Margin', 4);
        end
        apply_plot_theme(ax, char(theme));

        if ~isempty(lgd)
            leg = legend(ax, lgd, cellstr(lgdTxt), ...
                'Location', 'southoutside', 'NumColumns', min(2, numel(lgdTxt)), 'Box', 'off');
        else
            leg = [];
        end
        style_legend_for_theme(leg, char(theme));
    end
    save_fig_dual_safe(f, fullfile(outDir, "LagrangianAccel_peakTrigger_vs_growth_" + theme), plotOpts);
    close_if_needed(f, plotOpts);
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
figSize = heatmap_figure_size(plotOpts);
f = figure('Color', 'w', 'Position', [120 120 figSize]);
ax = axes(f, 'Position', heatmap_axes_position());
hold(ax, 'on');

hImg = imagesc(ax, xCenters, yCenters, Z);
set(hImg, 'AlphaData', isfinite(Z));
set(ax, 'YDir', 'normal', 'Color', 'w');
axis(ax, 'tight');
colormap(ax, lagrangian_heat_colormap(256, opts));
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
rng(opts.randomSeed + round(Rei), 'twister');

contourOverlays = struct('Xg', {}, 'Yg', {}, 'D', {}, 'color', {});
for j = 1:numel(caseData)
    d = caseData(j);
    xy = d.activationXYNorm;
    valid = all(isfinite(xy), 2) & xy(:,1) >= opts.xLimNorm(1) & xy(:,1) <= opts.xLimNorm(2) & ...
        xy(:,2) >= opts.yLimNorm(1) & xy(:,2) <= opts.yLimNorm(2);
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
    h = scatter(ax, xyPlot(:,1), xyPlot(:,2), 46, ...
        'Marker', marker, ...
        'MarkerFaceColor', cmap(j,:), ...
        'MarkerEdgeColor', [0.03 0.03 0.03], ...
        'LineWidth', 0.8, ...
        'MarkerFaceAlpha', 0.82);
    lgd(end+1,1) = h; %#ok<AGROW>
    lgdTxt(end+1,1) = sprintf('k/d=%.4g', d.kD); %#ok<AGROW>

    if size(xy, 1) >= 10
        [Xg, Yg] = meshgrid(linspace(opts.xLimNorm(1), opts.xLimNorm(2), 90), ...
            linspace(opts.yLimNorm(1), opts.yLimNorm(2), 90));
        D = kde2d_simple(xy(:,1), xy(:,2), Xg, Yg);
        contourOverlays(end+1).Xg = Xg; %#ok<AGROW>
        contourOverlays(end).Yg = Yg;
        contourOverlays(end).D = D;
        contourOverlays(end).color = cmap(j,:);
    end
end

for ci = 1:numel(contourOverlays)
    contour(ax, contourOverlays(ci).Xg, contourOverlays(ci).Yg, contourOverlays(ci).D, 1, ...
        'LineColor', contourOverlays(ci).color, ...
        'LineWidth', 2.4, ...
        'HandleVisibility', 'off');
end

xlim(ax, opts.xLimNorm);
ylim(ax, opts.yLimNorm);
if opts.heatmapPreserveSpatialAspect
    pbaspect(ax, [diff(opts.xLimNorm) diff(opts.yLimNorm) 1]);
end
xlabel(ax, '$x/H$', 'Interpreter', 'latex');
ylabel(ax, '$y/H$', 'Interpreter', 'latex');
title(ax, sprintf('Lagrangian acceleration heatmap (%s), Re = %g', upper(statName), Rei), 'FontName', fontName);
set(ax, 'FontName', fontName, 'Layer', 'top');
grid(ax, 'off');
box(ax, 'on');

if ~isempty(lgd)
    leg = legend(ax, lgd, cellstr(lgdTxt), 'Location', 'southoutside', ...
        'NumColumns', min(4, numel(lgdTxt)), 'Box', 'off');
else
    leg = [];
end

apply_plot_theme(ax, theme);
style_heatmap_white_theme(f, ax, cb);
style_legend_for_theme(leg, theme);
style_legend_for_heatmap_white(leg);
save_fig_dual_safe(f, fullfile(outDir, sprintf('LagrangianAccel_heatmap_%s_Re_%g_%s', statName, Rei, theme)), plotOpts);
close_if_needed(f, plotOpts);
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
        if strcmpi(statName, 'p90')
            Z(iy, ix) = finite_prctile(vals, 90);
        else
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
function cmap = lagrangian_heat_colormap(n, opts)
mapName = "sky";
if isfield(opts, 'heatmapColormap') && ~isempty(opts.heatmapColormap)
    mapName = string(opts.heatmapColormap);
end

switch lower(char(mapName))
    case 'abyss'
        cmap = abyss_colormap_compat(n);
        return;
    case 'sky'
        cmap = sky_colormap_compat(n);
        return;
end

try
    cmap = feval(char(mapName), n);
catch
    cmap = scientific_heat_colormap(n);
end
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
