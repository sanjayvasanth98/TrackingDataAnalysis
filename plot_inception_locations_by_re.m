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
nBinsX = 50;
nBinsY = 40;
histAlpha = 0.40;
histFrac = 0.12; % fraction of figure for marginal panels

for theme = reshape(plotOpts.themes, 1, [])
    for r = 1:numel(ReVals)
        Rei = ReVals(r);
        idxRe = find(allLoc.Re == Rei);
        nReCases = numel(idxRe);
        cmap = inception_colormap(nReCases);

        % --- Plot 1: dimensional (mm) ---
        plot_one_inception(allLoc, idxRe, nReCases, cmap, yExtent_mm, ...
            throatHeight_mm, xLim, yLim, false, nBinsX, nBinsY, histAlpha, histFrac, ...
            '$x\;(\mathrm{mm})$', '$y\;(\mathrm{mm})$', ...
            char(theme), ...
            fullfile(outDir, sprintf('Inception2x_locations_Re_%g_%s', Rei, theme)), ...
            plotOpts);

        % --- Plot 2: normalized by throat height (x/H, y/H) ---
        plot_one_inception(allLoc, idxRe, nReCases, cmap, yExtent_mm, ...
            throatHeight_mm, xLimNorm, yLimNorm, true, nBinsX, nBinsY, histAlpha, histFrac, ...
            '$x/H$', '$y/H$', ...
            char(theme), ...
            fullfile(outDir, sprintf('Inception2x_locations_normalized_Re_%g_%s', Rei, theme)), ...
            plotOpts);
    end
end

end

%% ---- main plotting helper ----
function plot_one_inception(allLoc, idxRe, nReCases, cmap, yExtent_mm, ...
    throatHeight_mm, xLim, yLim, doNormalize, nBinsX, nBinsY, histAlpha, histFrac, ...
    xLabel, yLabel, theme, outBase, plotOpts)

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
mainW = 1 - mainL - histFrac - gapX - padR;
mainH = 1 - mainB - histFrac - gapY - padT;
topL = mainL; topB = mainB + mainH + gapY; topW = mainW; topH = histFrac;
rightL = mainL + mainW + gapX; rightB = mainB; rightW = histFrac; rightH = mainH;

axMain  = axes(f, 'Position', [mainL  mainB  mainW  mainH]);
axTop   = axes(f, 'Position', [topL   topB   topW   topH]);
axRight = axes(f, 'Position', [rightL rightB rightW rightH]);

hold(axMain, 'on');
hold(axTop, 'on');
hold(axRight, 'on');

lgd = gobjects(0,1);
lgdTxt = strings(0,1);

xEdges = linspace(xLim(1), xLim(2), nBinsX+1);
yEdges = linspace(yLim(1), yLim(2), nBinsY+1);

% --- Draw wall region if ROI data provided ---
hasROI = isfield(plotOpts, 'roiData') && isstruct(plotOpts.roiData);
if hasROI
    draw_wall_patch(axMain, plotOpts.roiData, yExtent_mm, throatHeight_mm, ...
        xLim, yLim, doNormalize, theme);
end

for j = 1:nReCases
    ci = idxRe(j);
    xy = allLoc.inception2x_xy{ci};
    if isempty(xy), continue; end

    % Filter out points in unwanted track area
    if hasROI && isfield(plotOpts.roiData, 'unwantedTrackMask')
        xy = filter_unwanted_points(xy, plotOpts.roiData, yExtent_mm);
        if isempty(xy), continue; end
    end

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
        'MarkerFaceAlpha', 0.50);

    lgd(end+1,1) = h; %#ok<AGROW>
    lgdTxt(end+1,1) = sprintf('k/d=%.4g', allLoc.kD(ci)); %#ok<AGROW>

    % 2D KDE density contour overlay
    if numel(xPts) >= 10
        nGrid = 80;
        xGrid = linspace(xLim(1), xLim(2), nGrid);
        yGrid = linspace(yLim(1), yLim(2), nGrid);
        [Xg, Yg] = meshgrid(xGrid, yGrid);
        density = kde2d_simple(xPts, yPts, Xg, Yg);
        contour(axMain, Xg, Yg, density, 1, ...
            'LineColor', cmap(j,:), 'LineWidth', 2.0);
    end

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
set(axMain, 'XLimMode', 'manual', 'YLimMode', 'manual', 'FontName', 'Times New Roman');
grid(axMain, 'off'); box(axMain, 'on');

% Top marginal styling
xlim(axTop, xLim);
set(axTop, 'XTickLabel', [], 'XLimMode', 'manual', 'FontName', 'Times New Roman');
grid(axTop, 'off'); box(axTop, 'on');

% Right marginal styling
ylim(axRight, yLim);
set(axRight, 'YTickLabel', [], 'YLimMode', 'manual', 'FontName', 'Times New Roman');
% Remove the 0 tick so it doesn't overlap with main plot's right x-tick
rightXTicks = get(axRight, 'XTick');
rightXTicks(rightXTicks == 0) = [];
set(axRight, 'XTick', rightXTicks);
grid(axRight, 'off'); box(axRight, 'on');

% Legend inside main axes
if ~isempty(lgd)
    leg = legend(axMain, lgd, cellstr(lgdTxt), ...
        'Location', 'northwest', 'Box', 'off', 'FontName', 'Times New Roman', 'FontSize', 8);
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
