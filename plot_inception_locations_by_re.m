function plot_inception_locations_by_re(allLoc, outDir, plotOpts)

if nargin < 3 || ~isfield(plotOpts, 'themes') || isempty(plotOpts.themes)
    plotOpts.themes = "normal";
end

locationField = 'inception2x_xy';
if isfield(plotOpts, 'inceptionLocationField') && ~isempty(plotOpts.inceptionLocationField)
    locationField = char(string(plotOpts.inceptionLocationField));
end

outputStem = 'Inception2x_locations';
if isfield(plotOpts, 'inceptionLocationOutputStem') && ~isempty(plotOpts.inceptionLocationOutputStem)
    outputStem = char(string(plotOpts.inceptionLocationOutputStem));
end

warningLabel = 'activation points on left-moving tracks';
if isfield(plotOpts, 'inceptionLocationWarningLabel') && ~isempty(plotOpts.inceptionLocationWarningLabel)
    warningLabel = char(string(plotOpts.inceptionLocationWarningLabel));
end

plotDimensional = true;
if isfield(plotOpts, 'inceptionLocationPlotDimensional')
    plotDimensional = logical(plotOpts.inceptionLocationPlotDimensional);
end

plotNormalized = true;
if isfield(plotOpts, 'inceptionLocationPlotNormalized')
    plotNormalized = logical(plotOpts.inceptionLocationPlotNormalized);
end

if ~plotDimensional && ~plotNormalized
    warning('Both dimensional and normalized inception location plots are disabled. Skipping.');
    return;
end

nCases = numel(allLoc.caseName);
if nCases == 0
    warning('No cases in allLoc. Skipping inception location plot.');
    return;
end

if ~isfield(allLoc, locationField)
    warning('allLoc.%s was not found. Skipping inception location plot.', locationField);
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
    if ~isempty(allLoc.(locationField){i})
        hasPoints = true;
        break;
    end
end
if ~hasPoints
    warning('No %s found.', warningLabel);
    return;
end

% Marginal PDF settings
nPdfGridX = 200;
nPdfGridY = 200;
pdfLineWidth = 1.8;
marginalFrac = 0.12; % fraction of figure for marginal panels
plotVariants = resolve_inception_plot_variants(plotOpts);

for theme = reshape(plotOpts.themes, 1, [])
    for r = 1:numel(ReVals)
        Rei = ReVals(r);
        idxRe = find(allLoc.Re == Rei);
        nReCases = numel(idxRe);
        cmap = inception_colormap(nReCases);
        for vi = 1:numel(plotVariants)
            variant = plotVariants(vi);

            % --- Plot 1: dimensional (mm) ---
            if plotDimensional
                plot_one_inception(allLoc, idxRe, nReCases, cmap, yExtent_mm, ...
                    throatHeight_mm, xLim, yLim, false, nPdfGridX, nPdfGridY, pdfLineWidth, marginalFrac, ...
                    '$x\;(\mathrm{mm})$', '$y\;(\mathrm{mm})$', ...
                    char(theme), ...
                    fullfile(outDir, sprintf('%s%s_Re_%g_%s', outputStem, variant.fileSuffix, Rei, theme)), ...
                    plotOpts, variant, locationField);
            end

            % --- Plot 2: normalized by throat height (x/H, y/H) ---
            if plotNormalized
                plot_one_inception(allLoc, idxRe, nReCases, cmap, yExtent_mm, ...
                    throatHeight_mm, xLimNorm, yLimNorm, true, nPdfGridX, nPdfGridY, pdfLineWidth, marginalFrac, ...
                    '$x/H$', '$y/H$', ...
                    char(theme), ...
                    fullfile(outDir, sprintf('%s_normalized%s_Re_%g_%s', outputStem, variant.fileSuffix, Rei, theme)), ...
                    plotOpts, variant, locationField);
            end
        end
    end
end

end

%% ---- main plotting helper ----
function plot_one_inception(allLoc, idxRe, nReCases, cmap, yExtent_mm, ...
    throatHeight_mm, xLim, yLim, doNormalize, nPdfGridX, nPdfGridY, pdfLineWidth, marginalFrac, ...
    xLabel, yLabel, theme, outBase, plotOpts, variant, locationField)

if nargin < 20 || isempty(locationField)
    locationField = 'inception2x_xy';
end
locationField = char(string(locationField));
variant = complete_inception_plot_variant(variant);

fontName = resolve_plot_font_name();
f = figure('Color', 'w', 'Position', [120 120 plotOpts.inceptionImageSize_px]);

% Layout: main scatter + top marginal + right marginal
% Use separate gaps because figure aspect ratio is very wide (1280x320),
% so the same normalized gap is much larger horizontally than vertically.
figW = plotOpts.inceptionImageSize_px(1);
figH = plotOpts.inceptionImageSize_px(2);
gapPx = 10;  % uniform physical gap in pixels
gapX = gapPx / figW;
gapY = gapPx / figH;
padR = 0.01; padT = 0.01;
mainL = 0.08; mainB = 0.18;
mainW = 1 - mainL - marginalFrac - gapX - padR;
mainH = 1 - mainB - marginalFrac - gapY - padT;
topL = mainL; topB = mainB + mainH + gapY; topW = mainW; topH = marginalFrac;
rightL = mainL + mainW + gapX; rightB = mainB; rightW = marginalFrac; rightH = mainH;

axMain  = axes(f, 'Position', [mainL  mainB  mainW  mainH]);
axTop   = axes(f, 'Position', [topL   topB   topW   topH]);
axRight = axes(f, 'Position', [rightL rightB rightW rightH]);

hold(axMain, 'on');
hold(axTop, 'on');
hold(axRight, 'on');

lgd = gobjects(0,1);
lgdTxt = strings(0,1);
contourOverlays = struct('Xg', {}, 'Yg', {}, 'density', {}, 'color', {});
caseData = struct('caseOrder', {}, 'caseIndex', {}, 'xAll', {}, 'yAll', {}, 'displayMask', {});

xPdfGrid = linspace(xLim(1), xLim(2), nPdfGridX);
yPdfGrid = linspace(yLim(1), yLim(2), nPdfGridY);

% --- Draw wall region if ROI data provided ---
hasROI = isfield(plotOpts, 'roiData') && isstruct(plotOpts.roiData);
if hasROI
    draw_wall_patch(axMain, plotOpts.roiData, yExtent_mm, throatHeight_mm, ...
        xLim, yLim, doNormalize, theme);
end

for j = 1:nReCases
    ci = idxRe(j);
    xy = allLoc.(locationField){ci};
    if isempty(xy), continue; end

    % Filter out points in unwanted track area
    if hasROI && isfield(plotOpts.roiData, 'unwantedTrackMask')
        xy = filter_unwanted_points(xy, plotOpts.roiData, yExtent_mm);
        if isempty(xy), continue; end
    end

    if doNormalize
        xAll = xy(:,1) / throatHeight_mm;
        yAll = (yExtent_mm - xy(:,2)) / throatHeight_mm;
    else
        xAll = xy(:,1);
        yAll = yExtent_mm - xy(:,2);
    end

    caseData(end+1).caseOrder = j; %#ok<AGROW>
    caseData(end).caseIndex = ci;
    caseData(end).xAll = xAll;
    caseData(end).yAll = yAll;
    caseData(end).displayMask = true(size(xAll));
end

caseData = mark_display_points_by_variant(caseData, variant);

for dataIdx = 1:numel(caseData)
    j = caseData(dataIdx).caseOrder;
    ci = caseData(dataIdx).caseIndex;
    xAll = caseData(dataIdx).xAll;
    yAll = caseData(dataIdx).yAll;
    displayMask = caseData(dataIdx).displayMask;
    xPts = xAll(displayMask);
    yPts = yAll(displayMask);

    if variant.useAllForDensity
        xDensity = xAll;
        yDensity = yAll;
    else
        xDensity = xPts;
        yDensity = yPts;
    end

    marker = marker_for_case(j, variant);
    if isempty(xPts)
        h = scatter(axMain, NaN, NaN, variant.markerSize, ...
            'Marker', marker, ...
            'MarkerEdgeColor', marker_edge_color(variant), ...
            'MarkerFaceColor', cmap(j,:), ...
            'MarkerFaceAlpha', variant.markerAlpha, ...
            'LineWidth', variant.markerLineWidth);
    else
        h = scatter(axMain, xPts, yPts, variant.markerSize, ...
            'Marker', marker, ...
            'MarkerEdgeColor', marker_edge_color(variant), ...
            'MarkerFaceColor', cmap(j,:), ...
            'MarkerFaceAlpha', variant.markerAlpha, ...
            'LineWidth', variant.markerLineWidth);
    end

    lgd(end+1,1) = h; %#ok<AGROW>
    lgdTxt(end+1,1) = sprintf('k/d=%.4g', allLoc.kD(ci)); %#ok<AGROW>

    % 2D KDE density contour overlay
    if numel(xDensity) >= 10
        nGrid = 80;
        xGrid = linspace(xLim(1), xLim(2), nGrid);
        yGrid = linspace(yLim(1), yLim(2), nGrid);
        [Xg, Yg] = meshgrid(xGrid, yGrid);
        density = kde2d_simple(xDensity, yDensity, Xg, Yg);
        contourOverlays(end+1).Xg = Xg; %#ok<AGROW>
        contourOverlays(end).Yg = Yg;
        contourOverlays(end).density = density;
        contourOverlays(end).color = cmap(j,:);
    end

    % Marginal PDFs as smooth line curves for each case.
    xPdf = kde1d_simple(xDensity, xPdfGrid);
    yPdf = kde1d_simple(yDensity, yPdfGrid);
    plot(axTop, xPdfGrid, xPdf, ...
        'Color', cmap(j,:), 'LineWidth', pdfLineWidth, 'HandleVisibility', 'off');
    plot(axRight, yPdf, yPdfGrid, ...
        'Color', cmap(j,:), 'LineWidth', pdfLineWidth, 'HandleVisibility', 'off');
end

% Draw density/cluster contours last so they sit above every marker.
for ci = 1:numel(contourOverlays)
    contour(axMain, contourOverlays(ci).Xg, contourOverlays(ci).Yg, ...
        contourOverlays(ci).density, 1, ...
        'LineColor', contourOverlays(ci).color, ...
        'LineWidth', 2.0, ...
        'HandleVisibility', 'off');
end

% Main axes styling
xlim(axMain, xLim);
ylim(axMain, yLim);
xlabel(axMain, xLabel, 'Interpreter', 'latex');
ylabel(axMain, yLabel, 'Interpreter', 'latex');
set(axMain, 'XLimMode', 'manual', 'YLimMode', 'manual', 'FontName', fontName);
grid(axMain, 'off'); box(axMain, 'on');

% Top marginal styling
xlim(axTop, xLim);
set(axTop, 'XTickLabel', [], 'XLimMode', 'manual', 'FontName', fontName);
grid(axTop, 'off'); box(axTop, 'on');

% Right marginal styling
ylim(axRight, yLim);
set(axRight, 'YTickLabel', [], 'YLimMode', 'manual', 'FontName', fontName);
% Remove the 0 tick so it doesn't overlap with main plot's right x-tick
rightXTicks = get(axRight, 'XTick');
rightXTicks(rightXTicks == 0) = [];
set(axRight, 'XTick', rightXTicks);
grid(axRight, 'off'); box(axRight, 'on');

% Legend inside main axes
if ~isempty(lgd)
    leg = legend(axMain, lgd, cellstr(lgdTxt), ...
        'Location', 'northwest', 'Box', 'off', 'FontName', fontName, 'FontSize', 8);
else
    leg = [];
end

apply_plot_theme(axMain, theme);
apply_plot_theme(axTop, theme);
apply_plot_theme(axRight, theme);
set(axTop, 'XTickLabel', []);
set(axRight, 'YTickLabel', []);
style_legend_for_theme(leg, theme);

save_fig_dual_safe(f, outBase, plotOpts);
if ~isfield(plotOpts, 'keepFiguresOpen') || ~plotOpts.keepFiguresOpen
    close(f);
end
end

%% ---- theme helpers ----
function plotVariants = resolve_inception_plot_variants(plotOpts)
plotVariants = struct('fileSuffix', {}, 'maxAct', {}, 'useAllForDensity', {}, ...
    'markerByCase', {}, 'randomSeed', {}, 'markerSize', {}, ...
    'markerAlpha', {}, 'markerLineWidth', {}, 'sampleMode', {});

makeCapped = true;
if isfield(plotOpts, 'makeCappedActivationInceptionPlots')
    makeCapped = logical(plotOpts.makeCappedActivationInceptionPlots);
end

makeAll = false;
if isfield(plotOpts, 'makeAllActivationInceptionPlots')
    makeAll = logical(plotOpts.makeAllActivationInceptionPlots);
end

makeRandom100 = true;
if isfield(plotOpts, 'makeRandom100InceptionPlot')
    makeRandom100 = logical(plotOpts.makeRandom100InceptionPlot);
end

maxAct = Inf;
if isfield(plotOpts, 'maxActivationsPerCase') && isfinite(plotOpts.maxActivationsPerCase)
    maxAct = plotOpts.maxActivationsPerCase;
end

randomN = 100;
if isfield(plotOpts, 'randomInceptionTotalPoints') && isfinite(plotOpts.randomInceptionTotalPoints)
    randomN = max(0, round(plotOpts.randomInceptionTotalPoints));
elseif isfield(plotOpts, 'randomInceptionPointsPerCase') && isfinite(plotOpts.randomInceptionPointsPerCase)
    % Backward-compatible alias. The random100 plot now samples this many
    % points from the pooled Re group, not this many per k/d case.
    randomN = max(0, round(plotOpts.randomInceptionPointsPerCase));
end

randomSeed = 42;
if isfield(plotOpts, 'randomInceptionSeed') && isfinite(plotOpts.randomInceptionSeed)
    randomSeed = round(plotOpts.randomInceptionSeed);
end

if makeCapped
    plotVariants(end+1) = make_inception_plot_variant('', maxAct, false, false, randomSeed, 12, 0.50, 0.50, 'perCase'); %#ok<AGROW>
end

if makeAll
    plotVariants(end+1) = make_inception_plot_variant('_allactivations', Inf, false, false, randomSeed, 12, 0.50, 0.50, 'perCase'); %#ok<AGROW>
end

if makeRandom100
    plotVariants(end+1) = make_inception_plot_variant('_random100_inception', randomN, true, true, randomSeed, 26, 0.78, 0.75, 'total'); %#ok<AGROW>
end

if isempty(plotVariants)
    plotVariants = make_inception_plot_variant('', maxAct, false, false, randomSeed, 12, 0.50, 0.50, 'perCase');
end
end

function variant = make_inception_plot_variant(fileSuffix, maxAct, useAllForDensity, markerByCase, randomSeed, markerSize, markerAlpha, markerLineWidth, sampleMode)
variant = struct();
variant.fileSuffix = fileSuffix;
variant.maxAct = maxAct;
variant.useAllForDensity = useAllForDensity;
variant.markerByCase = markerByCase;
variant.randomSeed = randomSeed;
variant.markerSize = markerSize;
variant.markerAlpha = markerAlpha;
variant.markerLineWidth = markerLineWidth;
variant.sampleMode = sampleMode;
end

function variant = complete_inception_plot_variant(variant)
if isempty(variant) || ~isstruct(variant)
    variant = make_inception_plot_variant('', Inf, false, false, 42, 12, 0.50, 0.50, 'perCase');
    return;
end
variant = default_variant_field(variant, 'fileSuffix', '');
variant = default_variant_field(variant, 'maxAct', Inf);
variant = default_variant_field(variant, 'useAllForDensity', false);
variant = default_variant_field(variant, 'markerByCase', false);
variant = default_variant_field(variant, 'randomSeed', 42);
variant = default_variant_field(variant, 'markerSize', 12);
variant = default_variant_field(variant, 'markerAlpha', 0.50);
variant = default_variant_field(variant, 'markerLineWidth', 0.50);
variant = default_variant_field(variant, 'sampleMode', 'perCase');
end

function s = default_variant_field(s, name, value)
if ~isfield(s, name) || isempty(s.(name))
    s.(name) = value;
end
end

function caseData = mark_display_points_by_variant(caseData, variant)
if isempty(caseData)
    return;
end
maxAct = variant.maxAct;
if ~(isfinite(maxAct) && maxAct >= 0)
    return;
end
maxAct = round(maxAct);

if strcmpi(char(variant.sampleMode), 'total')
    totalCount = 0;
    for i = 1:numel(caseData)
        n = numel(caseData(i).xAll);
        caseData(i).displayMask = false(n, 1);
        totalCount = totalCount + n;
    end
    if totalCount == 0
        return;
    end
    nShow = min(maxAct, totalCount);
    rng(variant.randomSeed, 'twister');
    selectedGlobal = sort(randperm(totalCount, nShow));

    startIdx = 1;
    for i = 1:numel(caseData)
        n = numel(caseData(i).xAll);
        stopIdx = startIdx + n - 1;
        localIdx = selectedGlobal(selectedGlobal >= startIdx & selectedGlobal <= stopIdx) - startIdx + 1;
        caseData(i).displayMask(localIdx) = true;
        startIdx = stopIdx + 1;
    end
    return;
end

for i = 1:numel(caseData)
    n = numel(caseData(i).xAll);
    caseData(i).displayMask = true(n, 1);
    if n > maxAct
        rng(variant.randomSeed + 1009 * caseData(i).caseOrder, 'twister');
        idx = sort(randperm(n, maxAct));
        caseData(i).displayMask = false(n, 1);
        caseData(i).displayMask(idx) = true;
    end
end
end

function marker = marker_for_case(caseIndex, variant)
if ~variant.markerByCase
    marker = 'o';
    return;
end
markers = {'o', 'd', 's', '^', 'v', 'p', 'h', '>', '<'};
marker = markers{mod(caseIndex - 1, numel(markers)) + 1};
end

function edgeColor = marker_edge_color(variant)
if variant.markerByCase
    edgeColor = [0.05 0.05 0.05];
else
    edgeColor = 'none';
end
end

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

%% ---- distinct colormap for up to 6 cases ----
function cmap = inception_colormap(n)
% Hand-picked colors with good visual separation
colors = [ ...
    0.85  0.20  0.15;  % red
    0.00  0.45  0.75;  % blue
    0.20  0.65  0.15;  % green
    0.75  0.40  0.85;  % purple
    0.95  0.55  0.05;  % orange
    0.00  0.75  0.75;  % teal
    ];
if n <= size(colors, 1)
    cmap = colors(1:n, :);
else
    cmap = lines(n);
end
end

%% ---- simple 2D kernel density estimator ----
function density = kde2d_simple(x, y, Xg, Yg)
% Gaussian KDE on a grid using Silverman bandwidth per axis
x = x(:); y = y(:);
n = numel(x);

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

function density = kde1d_simple(x, xGrid)
% Gaussian KDE evaluated on a 1D grid for marginal location PDFs.
x = x(:);
x = x(isfinite(x));
xGrid = xGrid(:).';
n = numel(x);

if n == 0
    density = nan(size(xGrid));
    return;
end

h = silverman_bw(x);
density = zeros(size(xGrid));
for i = 1:n
    dx = (xGrid - x(i)) / h;
    density = density + exp(-0.5 * dx.^2);
end
density = density / (n * sqrt(2 * pi) * h);
end

function h = silverman_bw(x)
n = numel(x);
s = std(x, 0);
iqrVal = prctile_safe(x, 75) - prctile_safe(x, 25);
scale = min(s, iqrVal / 1.34);
if ~(isfinite(scale) && scale > 0)
    scale = max(s, iqrVal / 1.34);
end
if ~(isfinite(scale) && scale > 0)
    scale = 1;
end
h = 0.9 * scale * n^(-1/5);
end

function q = prctile_safe(x, p)
x = sort(x(isfinite(x)));
n = numel(x);
if n == 0, q = NaN; return; end
idx = 1 + (n - 1) * p / 100;
i0 = floor(idx); i1 = ceil(idx);
if i0 == i1
    q = x(i0);
else
    q = x(i0) + (idx - i0) * (x(i1) - x(i0));
end
end

%% ---- wall patch drawing ----
function draw_wall_patch(ax, roiData, yExtent_mm, throatHeight_mm, ...
    xLim, yLim, doNormalize, theme)
if ~isfield(roiData, 'wallMask'), return; end
wallMask = roiData.wallMask;
ps = roiData.maskPixelSize;

% For each column, find the topmost wall pixel (wall surface)
wallCols = find(any(wallMask, 1));
if isempty(wallCols), return; end

wallTopRow = zeros(size(wallCols));
for i = 1:numel(wallCols)
    wallTopRow(i) = find(wallMask(:, wallCols(i)), 1, 'first');
end

% Convert pixel to mm then to plot coords (y flipped)
wallX_mm = wallCols(:) * ps;
wallYSurface_mm = yExtent_mm - wallTopRow(:) * ps;  % plot y (flipped)
wallYBottom = yLim(1);  % bottom of plot

if doNormalize
    wallX = wallX_mm / throatHeight_mm;
    wallYSurface = wallYSurface_mm / throatHeight_mm;
    wallYBottom = wallYBottom / throatHeight_mm;
else
    wallX = wallX_mm;
    wallYSurface = wallYSurface_mm;
end

% Clip to x-axis limits
inRange = wallX >= xLim(1) & wallX <= xLim(2);
wallX = wallX(inRange);
wallYSurface = wallYSurface(inRange);
if isempty(wallX), return; end

% Build closed polygon: surface left→right, then bottom right→left
patchX = [wallX; flipud(wallX)];
patchY = [wallYSurface; repmat(wallYBottom, numel(wallX), 1)];

if strcmp(theme, 'poster')
    wallColor = [0.35 0.35 0.35];
else
    wallColor = [0.80 0.80 0.80];
end

patch(ax, patchX, patchY, wallColor, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.6, ...
    'HandleVisibility', 'off');
end

%% ---- filter points in unwanted track area ----
function xy = filter_unwanted_points(xy, roiData, yExtent_mm)
mask = roiData.unwantedTrackMask;
ps = roiData.maskPixelSize;
[nRows, nCols] = size(mask);

% Convert mm coordinates to pixel indices
% xy(:,1) = x_mm, xy(:,2) = y_mm (image coords, not flipped)
col_px = round(xy(:,1) / ps);
row_px = round(xy(:,2) / ps);

% Clamp to valid pixel range
col_px = max(1, min(col_px, nCols));
row_px = max(1, min(row_px, nRows));

% Check which points fall inside the unwanted mask
inUnwanted = false(size(xy, 1), 1);
for i = 1:size(xy, 1)
    inUnwanted(i) = mask(row_px(i), col_px(i));
end

xy(inUnwanted, :) = [];
end
