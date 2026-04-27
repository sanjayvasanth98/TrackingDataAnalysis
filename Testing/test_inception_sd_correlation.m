%% test_inception_sd_correlation.m
% Standalone test for comparing activation-location density clusters with
% standard-deviation TIFF maps.
%
% Method:
%   1. Load the same allLoc structure used by test_plot_inception_locations.
%   2. Read one 8-bit SD TIFF per case in native 0-255 intensity units.
%   3. Mark high-fluctuation SD pixels as intensity >= 102 (40% of 255).
%   4. Rebuild the activation-location KDE and its HDR cluster mask.
%   5. Overlay/display the SD cluster and activation-density cluster.
%
% Edit the path blocks below if your SD TIFF files are stored elsewhere.

clear; clc;

%% Paths
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot);

matDirCandidates = [
    "C:\Users\kbsanjayvasanth\Downloads\plot_data_mat (1)"
    "E:\March Re 90,000 inception data\Processed images\results\results 33 local\plot_data_mat"
    "/mnt/e/March Re 90,000 inception data/Processed images/results/results 33 local/plot_data_mat"
    ];

sdTiffDirCandidates = [
    "E:\March Re 90,000 inception data\Processed images\Standard deviation"
    "E:\March Re 90,000 inception data\Processed images\Standard deviation tiffs"
    "E:\March Re 90,000 inception data\Processed images"
    "/mnt/e/March Re 90,000 inception data/Processed images"
    ];

roiFileCandidates = [
    "E:\March Re 90,000 inception data\Processed images\results\2000 frames\Smooth variation 2000\ROI_throat.mat"
    "/mnt/e/March Re 90,000 inception data/Processed images/results/2000 frames/Smooth variation 2000/ROI_throat.mat"
    ];

outDir = fullfile(fileparts(mfilename('fullpath')), 'test_outputs', 'InceptionSDCorrelation');
if ~isfolder(outDir), mkdir(outDir); end

%% Case-to-SD TIFF matching
% If auto-matching does not find the correct TIFF, put the exact file path
% in the File column for that case. SearchToken is used only when File is
% empty.
sdTiffByCase = table( ...
    ["5um"; "12um"; "20um"; "30um"; "53um"; "80um"], ...
    strings(6, 1), ...
    ["Smooth"; "P10S100"; "P10S70"; "P10S50"; "P10S30"; "P10S20"], ...
    'VariableNames', {'Case', 'File', 'SearchToken'});

%% Analysis options
opts = struct();
opts.locationField = 'inception2x_xy';
opts.sdMaxValue = 255;                % 8-bit TIFF brightness range
opts.sdThresholdFraction = 0.40;
opts.sdThresholdValue = round(opts.sdThresholdFraction * opts.sdMaxValue); % 102
opts.activationClusterMassFraction = 0.50; % matches inceptionDensityContourMassFraction
opts.inceptionImageSize_px = [1280 320];
opts.inceptionXLim_mm = [0 4.8];
opts.inceptionYLim_mm = [0 1.2];
opts.densityGridSize = [80 320];      % [rows cols], long enough for the 4:1 field of view
opts.excludeRoiMasks = true;
opts.makeOverlayPlots = true;
opts.keepFiguresOpen = true;
set(0, 'DefaultFigureVisible', 'on');

%% Load activation locations
matDir = first_existing_folder(matDirCandidates);
if strlength(matDir) < 1
    error('Could not find a plot_data_mat folder. Edit matDirCandidates in this test script.');
end

matFile = fullfile(matDir, "inception_locations_by_case.mat");
if ~isfile(matFile)
    error('Missing inception location MAT file: %s', matFile);
end

S = load(matFile, 'allLoc');
allLoc = S.allLoc;
fprintf('Loaded %d cases from: %s\n', numel(allLoc.caseName), matFile);

if ~isfield(allLoc, opts.locationField)
    error('allLoc.%s was not found in %s.', opts.locationField, matFile);
end

roiData = load_roi_data(roiFileCandidates, allLoc);
if isempty(roiData)
    fprintf('ROI file not found. Continuing without ROI mask exclusion.\n');
else
    fprintf('Loaded ROI masks for point/mask exclusion.\n');
end

sdTiffDirs = existing_folders(sdTiffDirCandidates);
if isempty(sdTiffDirs)
    warning('No SD TIFF search folders exist. Exact files in sdTiffByCase.File will still be used.');
end

%% Run correlation per case
nCases = numel(allLoc.caseName);
rows = table();
missingCases = strings(0, 1);

for ci = 1:nCases
    caseName = string(allLoc.caseName(ci));
    xy = allLoc.(opts.locationField){ci};
    if isempty(xy)
        warning('Skipping %s: no activation points in allLoc.%s.', caseName, opts.locationField);
        continue;
    end

    pxSize = allLoc.pixelSize(ci);
    if ~(isfinite(pxSize) && pxSize > 0)
        warning('Skipping %s: invalid pixelSize.', caseName);
        continue;
    end

    sdFile = resolve_sd_tiff_file(caseName, sdTiffByCase, sdTiffDirs);
    if strlength(sdFile) < 1
        warning('Skipping %s: no matching SD TIFF found. Set sdTiffByCase.File for this case.', caseName);
        missingCases(end+1, 1) = caseName; %#ok<SAGROW>
        continue;
    end

    sdImg = read_tiff_grayscale_8bit(sdFile);
    if isempty(sdImg)
        warning('Skipping %s: SD TIFF could not be read: %s', caseName, sdFile);
        continue;
    end

    result = analyze_one_case_sd_correlation( ...
        allLoc, ci, xy, sdImg, sdFile, opts, roiData);

    rows = append_table_compat_local(rows, result.summaryRow);

    fprintf(['%s: n=%d, SD>=%.0f area=%.3f, point high-SD=%.3f, ', ...
        'cluster/high-SD Jaccard=%.3f, density-SD Pearson=%.3f\n'], ...
        caseName, result.summaryRow.nActivationPoints, result.summaryRow.sdThresholdValue, ...
        result.summaryRow.highSdAreaFraction, ...
        result.summaryRow.pointHighSdFraction, ...
        result.summaryRow.jaccardHighSdVsActivationCluster, ...
        result.summaryRow.densitySdPearsonR);

    if opts.makeOverlayPlots
        make_sd_correlation_overlay(result, allLoc, ci, opts, outDir);
    end
end

if ~isempty(rows)
    rows = sortrows(rows, {'Re', 'kD'});
    outCsv = fullfile(outDir, 'inception_sd_correlation_summary_test.csv');
    write_table_csv_compat_local(rows, outCsv);
    save(fullfile(outDir, 'inception_sd_correlation_test_input.mat'), ...
        'allLoc', 'opts', 'sdTiffByCase', 'rows');
    fprintf('\nWrote SD correlation summary: %s\n', outCsv);
else
    error(['No cases were analyzed. Most likely the SD TIFF names did not auto-match. ', ...
        'Edit sdTiffByCase.File in this script with the exact TIFF paths.']);
end

if ~isempty(missingCases)
    fprintf('\nCases missing SD TIFFs: %s\n', strjoin(cellstr(missingCases), ', '));
end

fprintf('Done. Output in: %s\n', outDir);

%% ------------------------------------------------------------------------
function result = analyze_one_case_sd_correlation(allLoc, ci, xyRaw, sdImg, sdFile, opts, roiData)

caseName = string(allLoc.caseName(ci));
pxSize = allLoc.pixelSize(ci);
yExtent_mm = opts.inceptionImageSize_px(2) * median(allLoc.pixelSize(isfinite(allLoc.pixelSize) & allLoc.pixelSize > 0));

if opts.excludeRoiMasks && ~isempty(roiData) && isfield(roiData, 'unwantedTrackMask')
    xy = filter_unwanted_points_local(xyRaw, roiData, yExtent_mm);
else
    xy = xyRaw;
end

[nRows, nCols] = size(sdImg);
gridRows = opts.densityGridSize(1);
gridCols = opts.densityGridSize(2);
[Xg, Yg] = meshgrid(linspace(1, nCols, gridCols), linspace(1, nRows, gridRows));

xPlot_mm = Xg * pxSize;
yPlot_mm = yExtent_mm - Yg * pxSize;
analysisMask = xPlot_mm >= opts.inceptionXLim_mm(1) & xPlot_mm <= opts.inceptionXLim_mm(2) & ...
    yPlot_mm >= opts.inceptionYLim_mm(1) & yPlot_mm <= opts.inceptionYLim_mm(2);

if opts.excludeRoiMasks && ~isempty(roiData)
    analysisMask = analysisMask & roi_keep_mask_on_grid(roiData, Xg, Yg, pxSize);
end

sdGrid = interp2(sdImg, Xg, Yg, 'linear', NaN);
highSdGrid = sdGrid >= opts.sdThresholdValue;
highSdGrid(~analysisMask) = false;

xy = xy(isfinite(xy(:,1)) & isfinite(xy(:,2)), :);
xPlotPts = xy(:,1);
yPlotPts = yExtent_mm - xy(:,2);
inPlot = xPlotPts >= opts.inceptionXLim_mm(1) & xPlotPts <= opts.inceptionXLim_mm(2) & ...
    yPlotPts >= opts.inceptionYLim_mm(1) & yPlotPts <= opts.inceptionYLim_mm(2);
xy = xy(inPlot, :);

colPts = xy(:,1) / pxSize;
rowPts = xy(:,2) / pxSize;
[colPts, rowPts] = keep_points_inside_mask(colPts, rowPts, analysisMask, nRows, nCols, gridRows, gridCols);

density = kde2d_grid(colPts, rowPts, Xg, Yg);
density(~analysisMask) = 0;

level = hdr_level_from_density(density, analysisMask, opts.activationClusterMassFraction);
activationCluster = density >= level & analysisMask;

validMask = analysisMask & isfinite(sdGrid) & isfinite(density);
sdVals = sdGrid(validMask);
densityVals = density(validMask);
highSdVals = highSdGrid(validMask);
clusterVals = activationCluster(validMask);

pointHighSdFraction = sample_high_sd_fraction(sdImg, xy, pxSize, opts.sdThresholdValue);
highSdAreaFraction = nnz(highSdGrid & validMask) / max(nnz(validMask), 1);
clusterAreaFraction = nnz(activationCluster & validMask) / max(nnz(validMask), 1);
pointHighSdEnrichment = pointHighSdFraction / max(highSdAreaFraction, eps);

intersectionCount = nnz(highSdGrid & activationCluster & validMask);
unionCount = nnz((highSdGrid | activationCluster) & validMask);
clusterCount = nnz(activationCluster & validMask);
highSdCount = nnz(highSdGrid & validMask);

jaccardVal = intersectionCount / max(unionCount, 1);
diceVal = 2 * intersectionCount / max(clusterCount + highSdCount, 1);
clusterCoveredByHighSd = intersectionCount / max(clusterCount, 1);
highSdCoveredByCluster = intersectionCount / max(highSdCount, 1);

densitySdPearson = pearson_corr_local(densityVals, sdVals);
densitySdSpearman = spearman_corr_local(densityVals, sdVals);
clusterHighSdPearson = pearson_corr_local(double(clusterVals), double(highSdVals));

if sum(densityVals) > 0
    densityWeightedMeanSd = sum(sdVals .* densityVals) / sum(densityVals);
else
    densityWeightedMeanSd = NaN;
end
roiMeanSd = mean_omitnan_local(sdVals);

summaryRow = table( ...
    caseName, allLoc.Re(ci), allLoc.kD(ci), string(sdFile), ...
    size(xy, 1), opts.sdThresholdFraction, opts.sdThresholdValue, opts.activationClusterMassFraction, ...
    highSdAreaFraction, clusterAreaFraction, clusterCoveredByHighSd, ...
    highSdCoveredByCluster, jaccardVal, diceVal, pointHighSdFraction, ...
    pointHighSdEnrichment, densitySdPearson, densitySdSpearman, ...
    clusterHighSdPearson, densityWeightedMeanSd, roiMeanSd, ...
    'VariableNames', {'Case', 'Re', 'kD', 'sdTiffFile', ...
    'nActivationPoints', 'sdThresholdFraction', 'sdThresholdValue', 'activationClusterMassFraction', ...
    'highSdAreaFraction', 'activationClusterAreaFraction', 'activationClusterCoveredByHighSdFraction', ...
    'highSdCoveredByActivationClusterFraction', 'jaccardHighSdVsActivationCluster', ...
    'diceHighSdVsActivationCluster', 'pointHighSdFraction', 'pointHighSdEnrichment', ...
    'densitySdPearsonR', 'densitySdSpearmanR', 'clusterMaskHighSdPearsonR', ...
    'densityWeightedMeanSd', 'roiMeanSd'});

result = struct();
result.caseName = caseName;
result.sdFile = string(sdFile);
result.sdGrid = sdGrid;
result.highSdGrid = highSdGrid;
result.density = density;
result.activationCluster = activationCluster;
result.analysisMask = analysisMask;
result.level = level;
result.Xg = Xg;
result.Yg = Yg;
result.xPlot_mm = xPlot_mm;
result.yPlot_mm = yPlot_mm;
result.xy = xy;
result.yExtent_mm = yExtent_mm;
result.summaryRow = summaryRow;
end

function make_sd_correlation_overlay(result, allLoc, ci, opts, outDir)
caseName = result.caseName;
token = sanitize_case_token_local(caseName);
pxSize = allLoc.pixelSize(ci);

xVec = result.Xg(1, :) * pxSize;
yVec = result.yExtent_mm - result.Yg(:, 1) * pxSize;
yVecAsc = flipud(yVec(:));

f = figure('Color', 'w', 'Position', [120 120 1200 360]);
if ~opts.keepFiguresOpen
    set(f, 'Visible', 'off');
end
ax = axes(f);

imagesc(ax, xVec, yVecAsc, flipud(result.sdGrid));
set(ax, 'YDir', 'normal');
colormap(ax, gray(256));
caxis(ax, [0 opts.sdMaxValue]);
hold(ax, 'on');

hSd = result.highSdGrid;
hAct = result.activationCluster;
hOverlap = hSd & hAct;
hSdOnly = hSd & ~hAct;
hActOnly = hAct & ~hSd;

add_cluster_overlay(ax, xVec, yVecAsc, hSdOnly, [0.00 0.45 1.00], 0.28);
add_cluster_overlay(ax, xVec, yVecAsc, hActOnly, [1.00 0.10 0.05], 0.28);
add_cluster_overlay(ax, xVec, yVecAsc, hOverlap, [0.70 0.10 0.95], 0.50);

contour(ax, xVec, yVecAsc, flipud(double(result.highSdGrid)), [0.5 0.5], ...
    'Color', [0.00 0.35 0.95], 'LineWidth', 1.4, 'DisplayName', sprintf('SD cluster >= %.0f', opts.sdThresholdValue));
contour(ax, xVec, yVecAsc, flipud(double(result.activationCluster)), [0.5 0.5], ...
    'Color', [0.90 0.05 0.02], 'LineWidth', 2.0, 'DisplayName', 'Activation density cluster');

xPts = result.xy(:,1);
yPts = result.yExtent_mm - result.xy(:,2);
scatter(ax, xPts, yPts, 10, 'MarkerFaceColor', [1.00 0.75 0.05], ...
    'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.35, 'DisplayName', 'Activation points');

xlim(ax, opts.inceptionXLim_mm);
ylim(ax, opts.inceptionYLim_mm);
xlabel(ax, 'x (mm)');
ylabel(ax, 'y (mm)');
title(ax, sprintf('%s: activation density vs SD map', caseName), 'Interpreter', 'none');
hSdOnlyL = patch(ax, NaN, NaN, [0.00 0.45 1.00], 'FaceAlpha', 0.28, ...
    'EdgeColor', 'none', 'DisplayName', 'SD cluster only');
hActOnlyL = patch(ax, NaN, NaN, [1.00 0.10 0.05], 'FaceAlpha', 0.28, ...
    'EdgeColor', 'none', 'DisplayName', 'Activation cluster only');
hOverlapL = patch(ax, NaN, NaN, [0.70 0.10 0.95], 'FaceAlpha', 0.50, ...
    'EdgeColor', 'none', 'DisplayName', 'Overlap');
legend(ax, [hSdOnlyL hActOnlyL hOverlapL], ...
    {'SD cluster only', 'Activation cluster only', 'Overlap'}, ...
    'Location', 'northoutside', 'Orientation', 'horizontal', 'Box', 'off');
box(ax, 'on');
cb = colorbar(ax);
ylabel(cb, 'SD intensity (0-255)');

outFile = fullfile(outDir, sprintf('Inception_SD_correlation_%s.png', token));
save_figure_png_local(f, outFile);
if ~opts.keepFiguresOpen
    close(f);
end
end

function add_cluster_overlay(ax, xVec, yVecAsc, mask, rgbColor, alphaVal)
maskForPlot = flipud(logical(mask));
rgb = zeros([size(maskForPlot), 3]);
for k = 1:3
    rgb(:,:,k) = rgbColor(k);
end
alpha = alphaVal * double(maskForPlot);
image(ax, xVec, yVecAsc, rgb, 'AlphaData', alpha, 'HitTest', 'off');
set(ax, 'YDir', 'normal');
end

function sdFile = resolve_sd_tiff_file(caseName, sdTiffByCase, sdTiffDirs)
sdFile = "";
idx = find(strcmpi(string(sdTiffByCase.Case), caseName), 1, 'first');
if ~isempty(idx)
    explicitFile = string(sdTiffByCase.File(idx));
    if strlength(explicitFile) > 0
        if isfile(explicitFile)
            sdFile = explicitFile;
            return;
        end
        for di = 1:numel(sdTiffDirs)
            candidate = fullfile(sdTiffDirs(di), explicitFile);
            if isfile(candidate)
                sdFile = string(candidate);
                return;
            end
        end
        warning('Explicit SD TIFF for %s was not found: %s', caseName, explicitFile);
    end
    searchToken = string(sdTiffByCase.SearchToken(idx));
else
    searchToken = caseName;
end

files = strings(0, 1);
for di = 1:numel(sdTiffDirs)
    files = [files; recursive_tiff_files(sdTiffDirs(di))]; %#ok<AGROW>
end
files = unique(files);
if isempty(files)
    return;
end

tokenLower = lower(char(searchToken));
sdTerms = ["std", "stdev", "standard", "deviation", "sd"];
bestScore = -Inf;
bestFile = "";
for fi = 1:numel(files)
    [~, base, ext] = fileparts(files(fi));
    baseLower = lower(char(base));
    hasToken = contains(baseLower, tokenLower);
    hasSdTerm = false;
    for ti = 1:numel(sdTerms)
        hasSdTerm = hasSdTerm || contains(baseLower, char(sdTerms(ti)));
    end
    isTiff = any(strcmpi(ext, {'.tif', '.tiff'}));
    if ~(hasToken && hasSdTerm && isTiff)
        continue;
    end

    score = 100 * double(hasToken) + 20 * double(hasSdTerm) - 0.001 * strlength(files(fi));
    if score > bestScore
        bestScore = score;
        bestFile = files(fi);
    end
end

sdFile = bestFile;
end

function files = recursive_tiff_files(rootDir)
files = strings(0, 1);
if strlength(rootDir) < 1 || ~isfolder(rootDir)
    return;
end

items = dir(rootDir);
for i = 1:numel(items)
    name = string(items(i).name);
    if items(i).isdir
        if name == "." || name == ".."
            continue;
        end
        files = [files; recursive_tiff_files(fullfile(rootDir, name))]; %#ok<AGROW>
    else
        [~, ~, ext] = fileparts(name);
        if strcmpi(ext, '.tif') || strcmpi(ext, '.tiff')
            files(end+1, 1) = string(fullfile(rootDir, name)); %#ok<AGROW>
        end
    end
end
end

function img8 = read_tiff_grayscale_8bit(fileName)
img8 = [];
try
    [raw, cmap] = imread(fileName);
catch
    return;
end

if exist('cmap', 'var') && ~isempty(cmap)
    raw = ind2rgb(raw, cmap);
end

if ndims(raw) == 3
    rawClass = class(raw);
    raw = double(raw);
    if strcmp(rawClass, 'logical')
        raw = 255 * raw;
    elseif isinteger(cast(0, rawClass))
        raw = raw / double(intmax(rawClass)) * 255;
    else
        finiteRaw = raw(isfinite(raw));
        if ~isempty(finiteRaw) && max(finiteRaw) <= 1
            raw = raw * 255;
        end
    end
    img8 = 0.2989 * raw(:,:,1) + 0.5870 * raw(:,:,2) + 0.1140 * raw(:,:,3);
elseif islogical(raw)
    img8 = 255 * double(raw);
elseif isinteger(raw)
    img8 = double(raw) / double(intmax(class(raw))) * 255;
else
    img8 = double(raw);
    finiteRaw = img8(isfinite(img8));
    if ~isempty(finiteRaw) && max(finiteRaw) <= 1
        img8 = img8 * 255;
    end
end

finiteVals = img8(isfinite(img8));
if isempty(finiteVals)
    img8 = [];
    return;
end
img8 = min(max(img8, 0), 255);
end

function density = kde2d_grid(x, y, Xg, Yg)
x = x(:);
y = y(:);
valid = isfinite(x) & isfinite(y);
x = x(valid);
y = y(valid);
n = numel(x);

density = zeros(size(Xg));
if n < 2
    return;
end

hx = silverman_bw_local(x);
hy = silverman_bw_local(y);
for i = 1:n
    dx = (Xg - x(i)) / hx;
    dy = (Yg - y(i)) / hy;
    density = density + exp(-0.5 * (dx.^2 + dy.^2));
end
density = density / (n * 2 * pi * hx * hy);
end

function level = hdr_level_from_density(density, mask, massFraction)
vals = density(mask & isfinite(density) & density > 0);
if isempty(vals)
    level = Inf;
    return;
end

vals = sort(vals(:), 'descend');
cumMass = cumsum(vals);
target = min(max(massFraction, 0.01), 0.99) * cumMass(end);
idx = find(cumMass >= target, 1, 'first');
if isempty(idx)
    idx = numel(vals);
end
level = vals(idx);
end

function h = silverman_bw_local(x)
x = x(isfinite(x));
n = numel(x);
if n < 2
    h = 1;
    return;
end
s = std(x, 0);
iqrVal = prctile_local(x, 75) - prctile_local(x, 25);
scale = min(s, iqrVal / 1.34);
if ~(isfinite(scale) && scale > 0)
    scale = max(s, iqrVal / 1.34);
end
if ~(isfinite(scale) && scale > 0)
    scale = 1;
end
h = max(0.9 * scale * n^(-1/5), eps);
end

function q = prctile_local(x, p)
x = sort(x(isfinite(x)));
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

function frac = sample_high_sd_fraction(sdImg, xy, pxSize, sdThreshold)
if isempty(xy)
    frac = NaN;
    return;
end
[nRows, nCols] = size(sdImg);
cols = round(xy(:,1) / pxSize);
rows = round(xy(:,2) / pxSize);
valid = rows >= 1 & rows <= nRows & cols >= 1 & cols <= nCols;
if ~any(valid)
    frac = NaN;
    return;
end
idx = sub2ind(size(sdImg), rows(valid), cols(valid));
frac = mean(sdImg(idx) >= sdThreshold);
end

function [colPts, rowPts] = keep_points_inside_mask(colPts, rowPts, mask, nRows, nCols, gridRows, gridCols)
if isempty(colPts)
    return;
end
gridCol = round(1 + (colPts - 1) * (gridCols - 1) / max(nCols - 1, 1));
gridRow = round(1 + (rowPts - 1) * (gridRows - 1) / max(nRows - 1, 1));
gridCol = max(1, min(gridCols, gridCol));
gridRow = max(1, min(gridRows, gridRow));
inside = false(size(colPts));
for i = 1:numel(colPts)
    inside(i) = mask(gridRow(i), gridCol(i));
end
colPts = colPts(inside);
rowPts = rowPts(inside);
end

function keepMask = roi_keep_mask_on_grid(roiData, Xg, Yg, sdPxSize)
keepMask = true(size(Xg));
if ~isfield(roiData, 'maskPixelSize') || ~(isfinite(roiData.maskPixelSize) && roiData.maskPixelSize > 0)
    return;
end

roiMask = false(size(roiData.unwantedTrackMask));
if isfield(roiData, 'wallMask') && isequal(size(roiData.wallMask), size(roiMask))
    roiMask = roiMask | logical(roiData.wallMask);
end
if isfield(roiData, 'unwantedTrackMask')
    roiMask = roiMask | logical(roiData.unwantedTrackMask);
end

[nRows, nCols] = size(roiMask);
xMm = Xg * sdPxSize;
yMm = Yg * sdPxSize;
roiCols = round(xMm / roiData.maskPixelSize);
roiRows = round(yMm / roiData.maskPixelSize);
roiCols = max(1, min(nCols, roiCols));
roiRows = max(1, min(nRows, roiRows));

bad = false(size(Xg));
for i = 1:numel(Xg)
    bad(i) = roiMask(roiRows(i), roiCols(i));
end
keepMask = ~bad;
end

function xy = filter_unwanted_points_local(xy, roiData, yExtent_mm) %#ok<INUSD>
if isempty(xy) || ~isfield(roiData, 'unwantedTrackMask') || ~isfield(roiData, 'maskPixelSize')
    return;
end
mask = roiData.unwantedTrackMask;
ps = roiData.maskPixelSize;
[nRows, nCols] = size(mask);
cols = round(xy(:,1) / ps);
rows = round(xy(:,2) / ps);
cols = max(1, min(nCols, cols));
rows = max(1, min(nRows, rows));

inUnwanted = false(size(xy, 1), 1);
for i = 1:size(xy, 1)
    inUnwanted(i) = mask(rows(i), cols(i));
end
xy(inUnwanted, :) = [];
end

function r = pearson_corr_local(x, y)
x = x(:);
y = y(:);
valid = isfinite(x) & isfinite(y);
x = x(valid);
y = y(valid);
if numel(x) < 2
    r = NaN;
    return;
end
x = x - mean(x);
y = y - mean(y);
den = sqrt(sum(x.^2) * sum(y.^2));
if den <= 0
    r = NaN;
else
    r = sum(x .* y) / den;
end
end

function r = spearman_corr_local(x, y)
x = x(:);
y = y(:);
valid = isfinite(x) & isfinite(y);
x = x(valid);
y = y(valid);
if numel(x) < 2
    r = NaN;
    return;
end
r = pearson_corr_local(tied_rank_local(x), tied_rank_local(y));
end

function mu = mean_omitnan_local(x)
x = x(isfinite(x));
if isempty(x)
    mu = NaN;
else
    mu = mean(x);
end
end

function ranks = tied_rank_local(x)
[xs, order] = sort(x(:));
n = numel(xs);
ranks = nan(n, 1);
i = 1;
while i <= n
    j = i;
    while j < n && xs(j + 1) == xs(i)
        j = j + 1;
    end
    ranks(order(i:j)) = (i + j) / 2;
    i = j + 1;
end
end

function roiData = load_roi_data(roiFileCandidates, allLoc)
roiData = [];
roiFile = first_existing_file(roiFileCandidates);
if strlength(roiFile) < 1
    return;
end
R = load(roiFile);
if ~isfield(R, 'wallMask') || ~isfield(R, 'unwantedTrackMask')
    return;
end
roiData = struct();
roiData.wallMask = R.wallMask;
roiData.unwantedTrackMask = R.unwantedTrackMask;
if isfield(R, 'ROI_throat') && isfield(R.ROI_throat, 'throat_xy_px')
    roiData.throat_xy_px = R.ROI_throat.throat_xy_px;
end
pixelVals = allLoc.pixelSize(isfinite(allLoc.pixelSize) & allLoc.pixelSize > 0);
roiData.maskPixelSize = median(pixelVals);
end

function folder = first_existing_folder(candidates)
folder = "";
for i = 1:numel(candidates)
    if isfolder(candidates(i))
        folder = string(candidates(i));
        return;
    end
end
end

function folders = existing_folders(candidates)
folders = strings(0, 1);
for i = 1:numel(candidates)
    if isfolder(candidates(i))
        folders(end+1, 1) = string(candidates(i)); %#ok<AGROW>
    end
end
folders = unique(folders);
end

function fileName = first_existing_file(candidates)
fileName = "";
for i = 1:numel(candidates)
    if isfile(candidates(i))
        fileName = string(candidates(i));
        return;
    end
end
end

function rows = append_table_compat_local(rows, row)
if isempty(rows)
    rows = row;
else
    rows = [rows; row]; %#ok<AGROW>
end
end

function write_table_csv_compat_local(T, fileName)
if exist('write_table_csv_compat', 'file') == 2
    write_table_csv_compat(T, fileName);
else
    writetable(T, fileName);
end
end

function save_figure_png_local(f, fileName)
if exist('exportgraphics', 'file') == 2
    exportgraphics(f, fileName, 'Resolution', 200);
else
    print(f, fileName, '-dpng', '-r200');
end
end

function token = sanitize_case_token_local(caseName)
token = char(string(caseName));
token = regexprep(token, '[^\w.-]+', '_');
if isempty(token)
    token = 'case';
end
end
