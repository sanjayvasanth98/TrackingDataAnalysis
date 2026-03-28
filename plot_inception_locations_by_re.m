function plot_inception_locations_by_re(allLoc, outDir, plotOpts)

if nargin < 3 || ~isfield(plotOpts, 'themes') || isempty(plotOpts.themes)
    plotOpts.themes = "normal";
end

nCases = numel(allLoc.caseName);
if nCases == 0
    warning('No cases in allLoc. Skipping inception location plot.');
    return;
end

ReVals = unique(allLoc.Re(:));
if isempty(ReVals)
    warning('No Reynolds numbers found. Skipping inception location plot.');
    return;
end

% Shared calibration assumption: all cases use the same pixel size.
pixelSizeVals = allLoc.pixelSize(isfinite(allLoc.pixelSize) & allLoc.pixelSize > 0);
if isempty(pixelSizeVals)
    warning('No valid pixel sizes found. Skipping inception location plot.');
    return;
end

yExtent_mm = plotOpts.inceptionImageSize_px(2) * median(pixelSizeVals);

% Throat height used for normalization
throatHeight_mm = 10;

xLim = [0 4.8];
if isfield(plotOpts, 'inceptionXLim_mm') && numel(plotOpts.inceptionXLim_mm) >= 2
    xLim = double(plotOpts.inceptionXLim_mm(1:2));
end
yLim = [0 1.2];
if isfield(plotOpts, 'inceptionYLim_mm') && numel(plotOpts.inceptionYLim_mm) >= 2
    yLim = double(plotOpts.inceptionYLim_mm(1:2));
end

% Normalized limits (x/H and y/H)
xLimNorm = xLim / throatHeight_mm;
xLimNorm(2) = min(xLimNorm(2), 0.5);  % restrict x/H to 0.5
yLimNorm = yLim / throatHeight_mm;

hasPoints = false;
for i = 1:nCases
    if ~isempty(allLoc.inception2x_xy{i})
        hasPoints = true;
        break;
    end
end
if ~hasPoints
    warning('No activation points found on left-moving tracks.');
    return;
end

% Marginal histogram settings
nBins = 40;
histAlpha = 0.35;
histFrac = 0.18; % fraction of figure for marginal panels

for theme = reshape(plotOpts.themes, 1, [])
    for r = 1:numel(ReVals)
        Rei = ReVals(r);
        idxRe = find(allLoc.Re == Rei);
        nReCases = numel(idxRe);
        cmap = lines(max(nReCases, 1));

        % --- Plot 1: dimensional (mm) ---
        plot_one_inception(allLoc, idxRe, nReCases, cmap, yExtent_mm, ...
            throatHeight_mm, xLim, yLim, false, nBins, histAlpha, histFrac, ...
            '$x\;(\mathrm{mm})$', '$y\;(\mathrm{mm})$', ...
            sprintf('Activation locations, Re=%g', Rei), ...
            char(theme), ...
            fullfile(outDir, sprintf('Inception2x_locations_Re_%g_%s', Rei, theme)), ...
            plotOpts);

        % --- Plot 2: normalized by throat height (x/H, y/H) ---
        plot_one_inception(allLoc, idxRe, nReCases, cmap, yExtent_mm, ...
            throatHeight_mm, xLimNorm, yLimNorm, true, nBins, histAlpha, histFrac, ...
            '$x/H$', '$y/H$', ...
            sprintf('Activation locations, Re=%g', Rei), ...
            char(theme), ...
            fullfile(outDir, sprintf('Inception2x_locations_normalized_Re_%g_%s', Rei, theme)), ...
            plotOpts);
    end
end

end

%% ---- main plotting helper ----
function plot_one_inception(allLoc, idxRe, nReCases, cmap, yExtent_mm, ...
    throatHeight_mm, xLim, yLim, doNormalize, nBins, histAlpha, histFrac, ...
    xLabel, yLabel, titleStr, theme, outBase, plotOpts)

f = figure('Color', 'w', 'Position', [120 120 plotOpts.inceptionImageSize_px], 'Visible', 'on');

% Layout: main scatter + top marginal + right marginal
gap = 0.02;
mainL = 0.10; mainB = 0.10;
mainW = 1 - mainL - histFrac - 3*gap;
mainH = 1 - mainB - histFrac - 3*gap;
topL = mainL; topB = mainB + mainH + gap; topW = mainW; topH = histFrac;
rightL = mainL + mainW + gap; rightB = mainB; rightW = histFrac; rightH = mainH;

axMain  = axes(f, 'Position', [mainL  mainB  mainW  mainH]);
axTop   = axes(f, 'Position', [topL   topB   topW   topH]);
axRight = axes(f, 'Position', [rightL rightB rightW rightH]);

hold(axMain, 'on');
hold(axTop, 'on');
hold(axRight, 'on');

lgd = gobjects(0,1);
lgdTxt = strings(0,1);

xEdges = linspace(xLim(1), xLim(2), nBins+1);
yEdges = linspace(yLim(1), yLim(2), nBins+1);

for j = 1:nReCases
    ci = idxRe(j);
    xy = allLoc.inception2x_xy{ci};
    if isempty(xy), continue; end

    if doNormalize
        xPts = xy(:,1) / throatHeight_mm;
        yPts = (yExtent_mm - xy(:,2)) / throatHeight_mm;
    else
        xPts = xy(:,1);
        yPts = yExtent_mm - xy(:,2);
    end

    h = scatter(axMain, xPts, yPts, 12, ...
        'Marker', 'o', ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceColor', cmap(j,:), ...
        'MarkerFaceAlpha', 0.35);

    lgd(end+1,1) = h; %#ok<AGROW>
    lgdTxt(end+1,1) = sprintf('k/d=%.4g', allLoc.kD(ci)); %#ok<AGROW>

    % Marginal histograms
    nxCounts = histcounts(xPts, xEdges);
    nyCounts = histcounts(yPts, yEdges);
    xCenters = (xEdges(1:end-1) + xEdges(2:end)) / 2;
    yCenters = (yEdges(1:end-1) + yEdges(2:end)) / 2;

    bar(axTop, xCenters, nxCounts, 1, ...
        'FaceColor', cmap(j,:), 'FaceAlpha', histAlpha, 'EdgeColor', cmap(j,:), 'LineWidth', 0.5);
    barh(axRight, yCenters, nyCounts, 1, ...
        'FaceColor', cmap(j,:), 'FaceAlpha', histAlpha, 'EdgeColor', cmap(j,:), 'LineWidth', 0.5);
end

% Main axes styling
xlim(axMain, xLim);
ylim(axMain, yLim);
xlabel(axMain, xLabel, 'Interpreter', 'latex');
ylabel(axMain, yLabel, 'Interpreter', 'latex');
title(axMain, titleStr, 'FontName', 'Times New Roman', 'FontSize', 12);
set(axMain, 'XLimMode', 'manual', 'YLimMode', 'manual');
grid(axMain, 'on'); box(axMain, 'on');
set(axMain, 'GridColor', [0.75 0.75 0.75], 'GridAlpha', 0.4);

% Top marginal styling
xlim(axTop, xLim);
set(axTop, 'XTickLabel', [], 'YTickLabel', [], 'XLimMode', 'manual');
set(axTop, 'GridColor', [0.75 0.75 0.75], 'GridAlpha', 0.4);
box(axTop, 'on');

% Right marginal styling
ylim(axRight, yLim);
set(axRight, 'XTickLabel', [], 'YTickLabel', [], 'YLimMode', 'manual');
set(axRight, 'GridColor', [0.75 0.75 0.75], 'GridAlpha', 0.4);
box(axRight, 'on');

% Legend inside main axes
if ~isempty(lgd)
    leg = legend(axMain, lgd, cellstr(lgdTxt), ...
        'Location', 'northwest', 'Box', 'off', 'FontSize', 8);
else
    leg = [];
end

apply_plot_theme(axMain, theme);
style_legend_for_theme(leg, theme);

save_fig_dual_safe(f, outBase, plotOpts);
close(f);
end

%% ---- theme helpers ----
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
