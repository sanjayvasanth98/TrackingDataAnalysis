%% test_inception_sd_correlation.m
% Standalone test for comparing activation-location density clusters with
% standard-deviation TIFF maps.
%
% Method:
%   1. Load the same allLoc structure used by test_plot_inception_locations.
%   2. Read one 8-bit SD TIFF per case in native 0-255 intensity units.
%   3. Mark high-fluctuation SD pixels using opts.sdThresholdValue.
%   4. Rebuild the activation-location KDE and its HDR cluster mask.
%   5. Overlay/display the SD cluster and activation-density cluster.
%
% Edit the path blocks below if your SD TIFF files are stored elsewhere.

clear; clc;

%% Paths
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot);
addpath(fileparts(mfilename('fullpath')));

matDir = test_plotmat_location("inception_locations_by_case.mat");

sdTiffDirCandidates = [ "E:\March Re 90,000 inception data\plots\Standard deviation plots\STD_Smoothtiff.tif", ...
    "E:\March Re 90,000 inception data\plots\Standard deviation plots\STD_P10S100.tif", ...
    "E:\March Re 90,000 inception data\plots\Standard deviation plots\STD_P10S70.tif", ...
    "E:\March Re 90,000 inception data\plots\Standard deviation plots\STD_P10S50.tif", ...
    "E:\March Re 90,000 inception data\plots\Standard deviation plots\STD_P10S30.tif", ...
    "E:\March Re 90,000 inception data\plots\Standard deviation plots\STD_P10S20.tif" ];

roiFileCandidates = ["E:\March Re 90,000 inception data\Processed images\results\2000 frames\Smooth variation 2000\ROI_throat.mat"];

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
opts.sharedClusterFraction = 0.50;    % use same displayed fraction for SD and activation cluster
opts.sdThresholdFraction = opts.sharedClusterFraction;
opts.sdThresholdValue = round(opts.sdThresholdFraction * opts.sdMaxValue);
opts.sdColormap = 'YlOrRd';
opts.activationClusterMassFraction = opts.sharedClusterFraction; % HDR density mass fraction
opts.inceptionImageSize_px = [1280 320];
opts.inceptionXLim_mm = [0 4.8];
opts.inceptionYLim_mm = [0 1.2];
opts.densityGridSize = [80 320];      % [rows cols], long enough for the 4:1 field of view
opts.excludeRoiMasks = true;
opts.makeOverlayPlots = true;
opts.keepFiguresOpen = true;
opts.overlayMaxInceptionPoints = 750;
opts.overlayRandomSeed = 42;
opts.keepLargestClusterOnly = true;
opts.makeThresholdSensitivitySweep = true;
opts.sdThresholdSweepFractions = 0.20:0.05:0.80;
opts.activationMassSweepFractions = 0.20:0.05:0.80;
opts.nullNumShuffles = 500;
opts.nullRandomSeed = 42;
opts.allowManualSdTiffSelection = true; % opens file picker when auto-match fails
opts.allowTokenOnlyTiffFallback = false; % set true only if SD TIFF names lack SD/std/deviation words
set(0, 'DefaultFigureVisible', 'on');

%% Load activation locations
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

sdTiffSearchCandidates = sdTiffDirCandidates(:);
explicitSdTiffFiles = existing_files(sdTiffSearchCandidates);
if isempty(explicitSdTiffFiles)
    sdTiffSearchCandidates = [sdTiffSearchCandidates; matDir; string(fileparts(matDir))];
end
sdTiffDirs = existing_folders(sdTiffSearchCandidates);
if isempty(sdTiffDirs) && isempty(explicitSdTiffFiles)
    warning('No SD TIFF search folders/files exist. Exact files in sdTiffByCase.File will still be used.');
    sdTiffFiles = strings(0, 1);
else
    if ~isempty(explicitSdTiffFiles)
        fprintf('Using %d explicit SD TIFF file(s):\n', numel(explicitSdTiffFiles));
        for fi = 1:numel(explicitSdTiffFiles)
            fprintf('  %s\n', explicitSdTiffFiles(fi));
        end
    end
    fprintf('Searching %d folder(s) for TIFF files:\n', numel(sdTiffDirs));
    for di = 1:numel(sdTiffDirs)
        fprintf('  %s\n', sdTiffDirs(di));
    end
    sdTiffFiles = unique([explicitSdTiffFiles; collect_tiff_files(sdTiffDirs)]);
    fprintf('Found %d available TIFF file(s).\n', numel(sdTiffFiles));
end

%% Run correlation per case
nCases = numel(allLoc.caseName);
rows = table();
sensitivityRows = table();
missingCases = strings(0, 1);
caseResults = struct([]);

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

    sdFile = resolve_sd_tiff_file(caseName, sdTiffByCase, sdTiffDirs, sdTiffFiles, opts);
    if strlength(sdFile) < 1 && opts.allowManualSdTiffSelection
        sdFile = prompt_for_sd_tiff_file(caseName, sdTiffDirs);
    end
    if strlength(sdFile) < 1
        warning(['Skipping %s: no matching SD TIFF found. Set sdTiffByCase.File for this case, ', ...
            'or set opts.allowTokenOnlyTiffFallback=true if the SD TIFF names do not contain SD/std/deviation.'], caseName);
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
    if isfield(result, 'sensitivityTable') && ~isempty(result.sensitivityTable)
        sensitivityRows = append_table_compat_local(sensitivityRows, result.sensitivityTable);
    end
    caseResults = append_struct_local(caseResults, result);

    fprintf(['%s: n=%d, SD>=%.0f area=%.3f, point high-SD=%.3f, ', ...
        'cluster/high-SD Jaccard=%.3f, density-SD Pearson=%.3f, streamwise r_x=%.3f (null p=%.3f)\n'], ...
        caseName, result.summaryRow.nActivationPoints, result.summaryRow.sdThresholdValue, ...
        result.summaryRow.highSdAreaFraction, ...
        result.summaryRow.pointHighSdFraction, ...
        result.summaryRow.jaccardHighSdVsActivationCluster, ...
        result.summaryRow.densitySdPearsonR, ...
        result.summaryRow.streamwiseProfilePearsonR, ...
        result.summaryRow.streamwiseNullPearsonP);

    if opts.makeOverlayPlots
        make_sd_correlation_overlay(result, allLoc, ci, opts, outDir, roiData);
    end
end

if ~isempty(rows)
    rows = sortrows(rows, {'Re', 'kD'});
    outCsv = fullfile(outDir, 'inception_sd_correlation_summary_test.csv');
    write_table_csv_compat_local(rows, outCsv);
    if ~isempty(sensitivityRows)
        sensitivityCsv = fullfile(outDir, 'inception_sd_threshold_sensitivity_test.csv');
        write_table_csv_compat_local(sensitivityRows, sensitivityCsv);
        make_threshold_sensitivity_heatmap(sensitivityRows, opts, outDir);
        fprintf('Wrote threshold sensitivity table: %s\n', sensitivityCsv);
    end
    save(fullfile(outDir, 'inception_sd_correlation_test_input.mat'), ...
        'allLoc', 'opts', 'sdTiffByCase', 'rows', 'sensitivityRows', 'caseResults');
    fprintf('\nWrote SD correlation summary: %s\n', outCsv);
    make_streamwise_profiles_all_cases(caseResults, opts, outDir);
    make_streamwise_xcorr_all_cases(caseResults, opts, outDir);
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
if isfield(opts, 'keepLargestClusterOnly') && logical(opts.keepLargestClusterOnly)
    highSdGrid = keep_largest_mask_component(highSdGrid);
end

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
if isfield(opts, 'keepLargestClusterOnly') && logical(opts.keepLargestClusterOnly)
    activationCluster = keep_largest_mask_component(activationCluster);
end

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
[streamX_mm, streamSdProfile, streamActivationProfile, streamPearson, streamSpearman, ...
    streamLag_mm, streamXcorr, streamXcorrPeakR, streamXcorrPeakLag_mm] = ...
    compute_streamwise_profiles(sdGrid, density, analysisMask, xPlot_mm, yPlot_mm);
[nullPearsonMean, nullPearsonStd, nullPearsonP, nullPeakMean, nullPeakStd, nullPeakP] = ...
    streamwise_circular_shift_null(streamSdProfile, streamActivationProfile, ...
    streamPearson, streamXcorrPeakR, opts, ci);
sensitivityTable = table();
if isfield(opts, 'makeThresholdSensitivitySweep') && logical(opts.makeThresholdSensitivitySweep)
    sensitivityTable = compute_threshold_sensitivity_table( ...
        caseName, allLoc.Re(ci), allLoc.kD(ci), sdGrid, density, analysisMask, opts);
end

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
    clusterHighSdPearson, streamPearson, streamSpearman, streamXcorrPeakR, streamXcorrPeakLag_mm, ...
    nullPearsonMean, nullPearsonStd, nullPearsonP, nullPeakMean, nullPeakStd, nullPeakP, ...
    densityWeightedMeanSd, roiMeanSd, ...
    'VariableNames', {'Case', 'Re', 'kD', 'sdTiffFile', ...
    'nActivationPoints', 'sdThresholdFraction', 'sdThresholdValue', 'activationClusterMassFraction', ...
    'highSdAreaFraction', 'activationClusterAreaFraction', 'activationClusterCoveredByHighSdFraction', ...
    'highSdCoveredByActivationClusterFraction', 'jaccardHighSdVsActivationCluster', ...
    'diceHighSdVsActivationCluster', 'pointHighSdFraction', 'pointHighSdEnrichment', ...
    'densitySdPearsonR', 'densitySdSpearmanR', 'clusterMaskHighSdPearsonR', ...
    'streamwiseProfilePearsonR', 'streamwiseProfileSpearmanR', ...
    'streamwiseXcorrPeakR', 'streamwiseXcorrPeakLag_mm', ...
    'streamwiseNullPearsonMean', 'streamwiseNullPearsonStd', 'streamwiseNullPearsonP', ...
    'streamwiseNullXcorrPeakMean', 'streamwiseNullXcorrPeakStd', 'streamwiseNullXcorrPeakP', ...
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
result.streamwiseX_mm = streamX_mm;
result.streamwiseSdProfile = streamSdProfile;
result.streamwiseActivationProfile = streamActivationProfile;
result.streamwiseProfilePearsonR = streamPearson;
result.streamwiseProfileSpearmanR = streamSpearman;
result.streamwiseLag_mm = streamLag_mm;
result.streamwiseXcorr = streamXcorr;
result.streamwiseXcorrPeakR = streamXcorrPeakR;
result.streamwiseXcorrPeakLag_mm = streamXcorrPeakLag_mm;
result.streamwiseNullPearsonMean = nullPearsonMean;
result.streamwiseNullPearsonStd = nullPearsonStd;
result.streamwiseNullPearsonP = nullPearsonP;
result.streamwiseNullXcorrPeakMean = nullPeakMean;
result.streamwiseNullXcorrPeakStd = nullPeakStd;
result.streamwiseNullXcorrPeakP = nullPeakP;
result.sensitivityTable = sensitivityTable;
result.summaryRow = summaryRow;
end

function make_sd_correlation_overlay(result, allLoc, ci, opts, outDir, roiData)
caseName = result.caseName;
token = sanitize_case_token_local(caseName);
pxSize = allLoc.pixelSize(ci);

xVec = result.Xg(1, :) * pxSize;
yVec = result.yExtent_mm - result.Yg(:, 1) * pxSize;
yVecAsc = flipud(yVec(:));

f = figure('Color', 'w', 'Position', [120 120 1160 340]);
if ~opts.keepFiguresOpen
    set(f, 'Visible', 'off');
end
ax = axes(f, 'Position', [0.065 0.17 0.82 0.72]);

imagesc(ax, xVec, yVecAsc, flipud(result.sdGrid));
set(ax, 'YDir', 'normal');
colormap(ax, sd_ylo_rd_colormap(256));
caxis(ax, [0 opts.sdMaxValue]);
hold(ax, 'on');

if nargin >= 6 && ~isempty(roiData)
    draw_wall_patch_overlay(ax, roiData, result.yExtent_mm, opts.inceptionXLim_mm, opts.inceptionYLim_mm);
end

draw_largest_binary_contour(ax, xVec, yVecAsc, result.highSdGrid, [0.05 0.05 0.05], 1.5);
draw_largest_binary_contour(ax, xVec, yVecAsc, result.activationCluster, [0.00 0.15 0.95], 2.0);

[xPts, yPts] = random_overlay_points(result, opts);
hPts = plot(ax, xPts, yPts, '+', ...
    'Color', [0.05 0.05 0.05], ...
    'MarkerSize', 4.5, ...
    'LineWidth', 0.75, ...
    'DisplayName', 'Activation points');

xlim(ax, opts.inceptionXLim_mm);
ylim(ax, opts.inceptionYLim_mm);
xlabel(ax, 'x (mm)');
ylabel(ax, 'y (mm)');
title(ax, sprintf('%s: activation density vs SD map', caseName), 'Interpreter', 'none');
hSdL = plot(ax, NaN, NaN, '-', 'Color', [0.05 0.05 0.05], ...
    'LineWidth', 1.5, 'DisplayName', 'SD cluster');
hActL = plot(ax, NaN, NaN, '-', 'Color', [0.00 0.15 0.95], ...
    'LineWidth', 2.0, 'DisplayName', 'Activation cluster');
legend(ax, [hSdL hActL hPts], ...
    {'SD cluster', 'Activation cluster', 'Activation points'}, ...
    'Location', 'northwest', 'Box', 'off', 'FontSize', 8);
box(ax, 'on');
cb = colorbar(ax);
ylabel(cb, 'SD intensity (0-255)');
set(ax, 'Position', [0.065 0.17 0.82 0.72]);
set(cb, 'Units', 'normalized', 'Position', [0.905 0.17 0.018 0.72]);
style_publication_axes(ax);
style_publication_colorbar(cb);

outFile = fullfile(outDir, sprintf('Inception_SD_correlation_%s.png', token));
save_figure_png_local(f, outFile);
if ~opts.keepFiguresOpen
    close(f);
end
end

function [xPts, yPts] = random_overlay_points(result, opts)
xy = result.xy;
if isempty(xy)
    xPts = [];
    yPts = [];
    return;
end

n = size(xy, 1);
nShow = n;
if isfield(opts, 'overlayMaxInceptionPoints') && isfinite(opts.overlayMaxInceptionPoints)
    nShow = min(n, max(0, round(opts.overlayMaxInceptionPoints)));
end

if nShow < n
    seed = 42;
    if isfield(opts, 'overlayRandomSeed') && isfinite(opts.overlayRandomSeed)
        seed = round(opts.overlayRandomSeed);
    end
    rng(seed + stable_string_seed(result.caseName), 'twister');
    idx = sort(randperm(n, nShow));
else
    idx = 1:n;
end

xPts = xy(idx, 1);
yPts = result.yExtent_mm - xy(idx, 2);
end

function draw_largest_binary_contour(ax, xVec, yVecAsc, mask, lineColor, lineWidth)
maskForPlot = flipud(double(logical(mask)));
contourMatrix = contourc(xVec, yVecAsc, maskForPlot, [0.5 0.5]);
segments = contour_matrix_to_segments_local(contourMatrix);
if isempty(segments)
    return;
end

bestIdx = 1;
bestScore = -Inf;
for si = 1:numel(segments)
    xy = segments{si};
    x = xy(1, :).';
    y = xy(2, :).';
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);
    if numel(x) < 2
        continue;
    end
    areaVal = 0;
    if numel(x) >= 3
        areaVal = abs(polyarea(x, y));
    end
    perimeter = sum(hypot(diff(x), diff(y)));
    score = areaVal + eps * perimeter;
    if score > bestScore
        bestScore = score;
        bestIdx = si;
    end
end

xy = segments{bestIdx};
plot(ax, xy(1, :), xy(2, :), '-', ...
    'Color', lineColor, ...
    'LineWidth', lineWidth, ...
    'HandleVisibility', 'off');
end

function segments = contour_matrix_to_segments_local(contourMatrix)
segments = {};
if isempty(contourMatrix)
    return;
end

col = 1;
while col < size(contourMatrix, 2)
    nPts = round(contourMatrix(2, col));
    startCol = col + 1;
    endCol = col + nPts;
    if nPts >= 2 && endCol <= size(contourMatrix, 2)
        segments{end+1} = contourMatrix(:, startCol:endCol); %#ok<AGROW>
    end
    col = endCol + 1;
end
end

function make_streamwise_profiles_all_cases(caseResults, opts, outDir)
if isempty(caseResults)
    return;
end

nCases = numel(caseResults);
colors = case_line_colormap(nCases);
f = figure('Color', 'w', 'Position', [140 120 1180 430]);
if ~opts.keepFiguresOpen
    set(f, 'Visible', 'off');
end
ax = axes(f, 'Position', [0.08 0.18 0.86 0.66]);
hold(ax, 'on');

handles = gobjects(0, 1);
labels = strings(0, 1);
for i = 1:nCases
    x = caseResults(i).streamwiseX_mm(:);
    sdProfile = normalize_profile_for_plot(caseResults(i).streamwiseSdProfile(:));
    actProfile = normalize_profile_for_plot(caseResults(i).streamwiseActivationProfile(:));
    c = colors(i, :);

    hS = plot(ax, x, sdProfile, '-', ...
        'Color', c, 'LineWidth', 2.0);
    plot(ax, x, actProfile, '--', ...
        'Color', c, 'LineWidth', 2.0);
    handles(end+1, 1) = hS; %#ok<AGROW>
    labels(end+1, 1) = sprintf('%s', char(caseResults(i).caseName)); %#ok<AGROW>
end

xlim(ax, opts.inceptionXLim_mm);
ylim(ax, [0 1.05]);
xlabel(ax, 'x (mm)');
ylabel(ax, 'Normalized profile');
title(ax, 'Streamwise profiles: solid = SD fluctuation, dashed = activation density', ...
    'Interpreter', 'none');
legend(ax, handles, cellstr(labels), ...
    'Location', 'northeast', 'Box', 'off', 'FontSize', 8);
style_publication_axes(ax);

save_figure_png_local(f, fullfile(outDir, 'Streamwise_profiles_all_cases.png'));
if ~opts.keepFiguresOpen
    close(f);
end
end

function make_streamwise_xcorr_all_cases(caseResults, opts, outDir)
if isempty(caseResults)
    return;
end

nCases = numel(caseResults);
colors = case_line_colormap(nCases);
f = figure('Color', 'w', 'Position', [160 120 620 620]);
if ~opts.keepFiguresOpen
    set(f, 'Visible', 'off');
end
ax = axes(f, 'Position', [0.14 0.13 0.78 0.78]);
hold(ax, 'on');

handles = gobjects(0, 1);
labels = strings(0, 1);
for i = 1:nCases
    lag = caseResults(i).streamwiseLag_mm(:);
    r = caseResults(i).streamwiseXcorr(:);
    c = colors(i, :);
    h = plot(ax, lag, r, '-', 'Color', c, 'LineWidth', 2.0);
    handles(end+1, 1) = h; %#ok<AGROW>
    labels(end+1, 1) = sprintf('%s peak %.2f at %.2f mm', ...
        char(caseResults(i).caseName), ...
        caseResults(i).streamwiseXcorrPeakR, ...
        caseResults(i).streamwiseXcorrPeakLag_mm); %#ok<AGROW>
end

xlim(ax, [-2 2]);
ylim(ax, [-1 1]);
pbaspect(ax, [1 1 1]);
plot(ax, [-2 2], [0 0], '-', 'Color', [0.65 0.65 0.65], ...
    'LineWidth', 0.8, 'HandleVisibility', 'off');
plot(ax, [0 0], [-1 1], '-', 'Color', [0.80 0.80 0.80], ...
    'LineWidth', 0.8, 'HandleVisibility', 'off');
xlabel(ax, 'Streamwise lag, \Delta x (mm)');
ylabel(ax, 'Normalized cross-correlation');
title(ax, 'Streamwise normalized cross-correlation', 'Interpreter', 'none');
legend(ax, handles, cellstr(labels), ...
    'Location', 'northeast', 'Box', 'off', 'FontSize', 8);
style_publication_axes(ax);

save_figure_png_local(f, fullfile(outDir, 'Streamwise_normalized_xcorr_all_cases.png'));
if ~opts.keepFiguresOpen
    close(f);
end
end

function sensitivityTable = compute_threshold_sensitivity_table(caseName, ReVal, kDVal, sdGrid, density, analysisMask, opts)
sdFracs = opts.sdThresholdSweepFractions(:);
actFracs = opts.activationMassSweepFractions(:);
nRows = numel(sdFracs) * numel(actFracs);

Case = strings(nRows, 1);
Re = nan(nRows, 1);
kD = nan(nRows, 1);
sdThresholdFraction = nan(nRows, 1);
sdThresholdValue = nan(nRows, 1);
activationClusterMassFraction = nan(nRows, 1);
highSdAreaFraction = nan(nRows, 1);
activationClusterAreaFraction = nan(nRows, 1);
activationClusterCoveredByHighSdFraction = nan(nRows, 1);
highSdCoveredByActivationClusterFraction = nan(nRows, 1);
jaccardHighSdVsActivationCluster = nan(nRows, 1);
diceHighSdVsActivationCluster = nan(nRows, 1);

validMask = analysisMask & isfinite(sdGrid) & isfinite(density);
row = 0;
for si = 1:numel(sdFracs)
    sdFrac = sdFracs(si);
    highSdGrid = sdGrid >= round(sdFrac * opts.sdMaxValue);
    highSdGrid(~analysisMask) = false;
    if isfield(opts, 'keepLargestClusterOnly') && logical(opts.keepLargestClusterOnly)
        highSdGrid = keep_largest_mask_component(highSdGrid);
    end

    for ai = 1:numel(actFracs)
        actFrac = actFracs(ai);
        level = hdr_level_from_density(density, analysisMask, actFrac);
        activationCluster = density >= level & analysisMask;
        if isfield(opts, 'keepLargestClusterOnly') && logical(opts.keepLargestClusterOnly)
            activationCluster = keep_largest_mask_component(activationCluster);
        end

        intersectionCount = nnz(highSdGrid & activationCluster & validMask);
        unionCount = nnz((highSdGrid | activationCluster) & validMask);
        highSdCount = nnz(highSdGrid & validMask);
        clusterCount = nnz(activationCluster & validMask);
        validCount = max(nnz(validMask), 1);

        row = row + 1;
        Case(row) = caseName;
        Re(row) = ReVal;
        kD(row) = kDVal;
        sdThresholdFraction(row) = sdFrac;
        sdThresholdValue(row) = round(sdFrac * opts.sdMaxValue);
        activationClusterMassFraction(row) = actFrac;
        highSdAreaFraction(row) = highSdCount / validCount;
        activationClusterAreaFraction(row) = clusterCount / validCount;
        activationClusterCoveredByHighSdFraction(row) = intersectionCount / max(clusterCount, 1);
        highSdCoveredByActivationClusterFraction(row) = intersectionCount / max(highSdCount, 1);
        jaccardHighSdVsActivationCluster(row) = intersectionCount / max(unionCount, 1);
        diceHighSdVsActivationCluster(row) = 2 * intersectionCount / max(clusterCount + highSdCount, 1);
    end
end

sensitivityTable = table(Case, Re, kD, sdThresholdFraction, sdThresholdValue, ...
    activationClusterMassFraction, highSdAreaFraction, activationClusterAreaFraction, ...
    activationClusterCoveredByHighSdFraction, highSdCoveredByActivationClusterFraction, ...
    jaccardHighSdVsActivationCluster, diceHighSdVsActivationCluster);
end

function make_threshold_sensitivity_heatmap(sensitivityRows, opts, outDir)
if isempty(sensitivityRows)
    return;
end

sdFracs = unique(sensitivityRows.sdThresholdFraction);
actFracs = unique(sensitivityRows.activationClusterMassFraction);
meanJaccard = nan(numel(actFracs), numel(sdFracs));
meanDice = nan(numel(actFracs), numel(sdFracs));

for si = 1:numel(sdFracs)
    for ai = 1:numel(actFracs)
        idx = sensitivityRows.sdThresholdFraction == sdFracs(si) & ...
            sensitivityRows.activationClusterMassFraction == actFracs(ai);
        meanJaccard(ai, si) = mean_omitnan_local(sensitivityRows.jaccardHighSdVsActivationCluster(idx));
        meanDice(ai, si) = mean_omitnan_local(sensitivityRows.diceHighSdVsActivationCluster(idx));
    end
end

finiteOverlapVals = [meanJaccard(isfinite(meanJaccard)); meanDice(isfinite(meanDice)); eps];
overlapCLim = [0 max(finiteOverlapVals)];

f = figure('Color', 'w', 'Position', [180 160 880 370]);
if ~opts.keepFiguresOpen
    set(f, 'Visible', 'off');
end

ax1 = axes(f, 'Position', [0.075 0.17 0.33 0.72]);
imagesc(ax1, sdFracs, actFracs, meanJaccard);
set(ax1, 'YDir', 'normal');
colormap(ax1, parula(256));
caxis(ax1, overlapCLim);
hold(ax1, 'on');
plot(ax1, opts.sdThresholdFraction, opts.activationClusterMassFraction, '+', ...
    'Color', [0 0 0], 'MarkerSize', 9, 'LineWidth', 1.4);
xlabel(ax1, 'SD threshold fraction');
ylabel(ax1, 'Activation HDR mass fraction');
title(ax1, 'Mean Jaccard');
cb1 = colorbar(ax1);
set(ax1, 'Position', [0.075 0.17 0.33 0.72]);
set(cb1, 'Units', 'normalized', 'Position', [0.415 0.17 0.018 0.72]);
style_publication_axes(ax1);
style_publication_colorbar(cb1);

ax2 = axes(f, 'Position', [0.555 0.17 0.33 0.72]);
imagesc(ax2, sdFracs, actFracs, meanDice);
set(ax2, 'YDir', 'normal');
colormap(ax2, parula(256));
caxis(ax2, overlapCLim);
hold(ax2, 'on');
plot(ax2, opts.sdThresholdFraction, opts.activationClusterMassFraction, '+', ...
    'Color', [0 0 0], 'MarkerSize', 9, 'LineWidth', 1.4);
xlabel(ax2, 'SD threshold fraction');
ylabel(ax2, 'Activation HDR mass fraction');
title(ax2, 'Mean Dice');
cb2 = colorbar(ax2);
set(ax2, 'Position', [0.555 0.17 0.33 0.72]);
set(cb2, 'Units', 'normalized', 'Position', [0.895 0.17 0.018 0.72]);
style_publication_axes(ax2);
style_publication_colorbar(cb2);

save_figure_png_local(f, fullfile(outDir, 'Threshold_sensitivity_mean_overlap.png'));
if ~opts.keepFiguresOpen
    close(f);
end
end

function [nullPearsonMean, nullPearsonStd, nullPearsonP, nullPeakMean, nullPeakStd, nullPeakP] = ...
    streamwise_circular_shift_null(sdProfile, activationProfile, observedPearson, observedPeak, opts, caseIndex)
nNull = 500;
if isfield(opts, 'nullNumShuffles') && isfinite(opts.nullNumShuffles)
    nNull = max(0, round(opts.nullNumShuffles));
end

nullPearson = nan(nNull, 1);
nullPeak = nan(nNull, 1);
sdProfile = sdProfile(:);
activationProfile = activationProfile(:);
n = min(numel(sdProfile), numel(activationProfile));
sdProfile = sdProfile(1:n);
activationProfile = activationProfile(1:n);

if nNull < 1 || n < 3
    nullPearsonMean = NaN;
    nullPearsonStd = NaN;
    nullPearsonP = NaN;
    nullPeakMean = NaN;
    nullPeakStd = NaN;
    nullPeakP = NaN;
    return;
end

seed = 42;
if isfield(opts, 'nullRandomSeed') && isfinite(opts.nullRandomSeed)
    seed = round(opts.nullRandomSeed);
end
rng(seed + 1009 * caseIndex, 'twister');

for k = 1:nNull
    shift = randi(max(n - 1, 1));
    shiftedActivation = circshift(activationProfile, shift);
    valid = isfinite(sdProfile) & isfinite(shiftedActivation);
    if nnz(valid) >= max(10, round(0.20 * n))
        nullPearson(k) = pearson_corr_local(sdProfile(valid), shiftedActivation(valid));
    end
    [~, rLag] = normalized_cross_correlation_profile((1:n).', sdProfile, shiftedActivation);
    if any(isfinite(rLag))
        nullPeak(k) = max(rLag(isfinite(rLag)));
    end
end

nullPearsonMean = mean_omitnan_local(nullPearson);
nullPearsonStd = std_omitnan_local(nullPearson);
nullPeakMean = mean_omitnan_local(nullPeak);
nullPeakStd = std_omitnan_local(nullPeak);
nullPearsonP = one_sided_upper_tail_p(nullPearson, observedPearson);
nullPeakP = one_sided_upper_tail_p(nullPeak, observedPeak);
end

function p = normalize_profile_for_plot(p)
p = p(:);
finiteVals = p(isfinite(p));
if isempty(finiteVals)
    p(:) = NaN;
    return;
end
minVal = min(finiteVals);
maxVal = max(finiteVals);
if maxVal > minVal
    p = (p - minVal) ./ (maxVal - minVal);
else
    p = zeros(size(p));
end
end

function draw_wall_patch_overlay(ax, roiData, yExtent_mm, xLim, yLim)
if ~isfield(roiData, 'wallMask') || ~isfield(roiData, 'maskPixelSize')
    return;
end

wallMask = roiData.wallMask;
ps = roiData.maskPixelSize;
wallCols = find(any(wallMask, 1));
if isempty(wallCols)
    return;
end

wallTopRow = zeros(size(wallCols));
for i = 1:numel(wallCols)
    wallTopRow(i) = find(wallMask(:, wallCols(i)), 1, 'first');
end

wallX = wallCols(:) * ps;
wallYSurface = yExtent_mm - wallTopRow(:) * ps;
inRange = wallX >= xLim(1) & wallX <= xLim(2);
wallX = wallX(inRange);
wallYSurface = wallYSurface(inRange);
if isempty(wallX)
    return;
end

patchX = [wallX; flipud(wallX)];
patchY = [wallYSurface; repmat(yLim(1), numel(wallX), 1)];
patch(ax, patchX, patchY, [0.80 0.80 0.80], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.60, ...
    'HandleVisibility', 'off');
end

function [xProfile_mm, sdProfile, activationProfile, pearsonR, spearmanR, ...
    lag_mm, xcorrProfile, xcorrPeakR, xcorrPeakLag_mm] = ...
    compute_streamwise_profiles(sdGrid, density, analysisMask, xPlot_mm, yPlot_mm)
nCols = size(sdGrid, 2);
xProfile_mm = xPlot_mm(1, :).';
sdProfile = nan(nCols, 1);
activationProfile = nan(nCols, 1);

dyVals = abs(diff(yPlot_mm(:, 1)));
dyVals = dyVals(isfinite(dyVals) & dyVals > 0);
if isempty(dyVals)
    dyMm = 1;
else
    dyMm = median(dyVals);
end

for j = 1:nCols
    keep = analysisMask(:, j) & isfinite(sdGrid(:, j)) & isfinite(density(:, j));
    if ~any(keep)
        continue;
    end
    sdProfile(j) = mean(sdGrid(keep, j));
    activationProfile(j) = sum(density(keep, j)) * dyMm;
end

valid = isfinite(sdProfile) & isfinite(activationProfile);
pearsonR = pearson_corr_local(sdProfile(valid), activationProfile(valid));
spearmanR = spearman_corr_local(sdProfile(valid), activationProfile(valid));
[lag_mm, xcorrProfile] = normalized_cross_correlation_profile(xProfile_mm, sdProfile, activationProfile);
if any(isfinite(xcorrProfile))
    finiteIdx = find(isfinite(xcorrProfile));
    [xcorrPeakR, relPeak] = max(xcorrProfile(finiteIdx));
    idxPeak = finiteIdx(relPeak);
    xcorrPeakLag_mm = lag_mm(idxPeak);
else
    xcorrPeakR = NaN;
    xcorrPeakLag_mm = NaN;
end
end

function [lag_mm, rLag] = normalized_cross_correlation_profile(x, s, a)
x = x(:);
s = s(:);
a = a(:);
n = min([numel(x), numel(s), numel(a)]);
x = x(1:n);
s = s(1:n);
a = a(1:n);

dxVals = diff(x);
dxVals = dxVals(isfinite(dxVals) & dxVals > 0);
if isempty(dxVals)
    dx = 1;
else
    dx = median(dxVals);
end

lags = (-(n - 1):(n - 1)).';
lag_mm = lags * dx;
rLag = nan(size(lags));
minOverlap = max(10, round(0.20 * n));

for ii = 1:numel(lags)
    lag = lags(ii);
    if lag >= 0
        sUse = s(1:n-lag);
        aUse = a(1+lag:n);
    else
        shift = -lag;
        sUse = s(1+shift:n);
        aUse = a(1:n-shift);
    end
    valid = isfinite(sUse) & isfinite(aUse);
    if nnz(valid) < minOverlap
        continue;
    end
    rLag(ii) = pearson_corr_local(sUse(valid), aUse(valid));
end
end

function cmap = sd_ylo_rd_colormap(n)
if nargin < 1 || isempty(n)
    n = 256;
end
ensure_cbrewer2_on_path_local();
if exist('cbrewer2', 'file') == 2
    try
        cmap = cbrewer2('seq', 'YlOrRd', n, 'pchip', 'rgb');
        cmap = normalize_colormap_rgb_local(cmap);
        return;
    catch
        try
            cmap = cbrewer2('YlOrRd', n);
            cmap = normalize_colormap_rgb_local(cmap);
            return;
        catch
            % Fall back to the built-in ColorBrewer anchor table below.
        end
    end
end

anchors = [ ...
    255 255 204
    255 237 160
    254 217 118
    254 178 76
    253 141 60
    252 78 42
    227 26 28
    189 0 38
    128 0 38] ./ 255;
x = linspace(0, 1, size(anchors, 1));
xi = linspace(0, 1, n);
cmap = interp1(x, anchors, xi, 'pchip');
cmap = normalize_colormap_rgb_local(cmap);
end

function ensure_cbrewer2_on_path_local()
if exist('cbrewer2', 'file') == 2
    return;
end
thisFile = mfilename('fullpath');
repoRootLocal = fileparts(fileparts(thisFile));
cbrewer2Dir = fullfile(repoRootLocal, 'external', 'cbrewer2', 'cbrewer2');
if isfolder(cbrewer2Dir)
    addpath(cbrewer2Dir);
end
end

function cmap = normalize_colormap_rgb_local(cmap)
cmap = double(cmap);
if isempty(cmap) || size(cmap, 2) ~= 3
    cmap = parula(256);
    return;
end
if max(cmap(:)) > 1
    cmap = cmap ./ 255;
end
cmap = min(max(cmap, 0), 1);
end

function cmap = case_line_colormap(n)
base = [ ...
    0.00 0.24 1.00
    0.95 0.10 0.05
    0.00 0.65 0.18
    0.62 0.18 0.95
    1.00 0.55 0.00
    0.00 0.72 0.95
    0.90 0.00 0.55
    0.15 0.15 0.15];
if n <= size(base, 1)
    cmap = base(1:n, :);
else
    cmap = lines(n);
end
end

function style_publication_axes(ax)
fontName = 'Arial';
if exist('resolve_plot_font_name', 'file') == 2
    try
        fontName = resolve_plot_font_name();
    catch
        fontName = 'Arial';
    end
end
set(ax, ...
    'FontName', fontName, ...
    'FontSize', 10, ...
    'LineWidth', 0.85, ...
    'TickDir', 'in', ...
    'Layer', 'top', ...
    'Box', 'on');
grid(ax, 'off');
end

function style_publication_colorbar(cb)
fontName = 'Arial';
if exist('resolve_plot_font_name', 'file') == 2
    try
        fontName = resolve_plot_font_name();
    catch
        fontName = 'Arial';
    end
end
set(cb, ...
    'FontName', fontName, ...
    'FontSize', 10, ...
    'LineWidth', 0.85);
if isprop(cb, 'TickDirection')
    set(cb, 'TickDirection', 'in');
end
end

function sdFile = resolve_sd_tiff_file(caseName, sdTiffByCase, sdTiffDirs, files, opts)
sdFile = "";
if nargin < 5 || isempty(opts)
    opts = struct();
end
allowTokenOnly = isfield(opts, 'allowTokenOnlyTiffFallback') && logical(opts.allowTokenOnlyTiffFallback);

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

files = unique(files);
if isempty(files)
    return;
end

tokenLower = lower(char(searchToken));
sdTerms = ["std", "stdev", "standard", "deviation", "sd"];
bestScore = -Inf;
bestFile = "";
tokenOnlyFiles = strings(0, 1);
tokenOnlyScores = [];
for fi = 1:numel(files)
    [~, base, ext] = fileparts(files(fi));
    baseLower = lower(char(base));
    hasToken = contains(baseLower, tokenLower);
    hasSdTerm = false;
    for ti = 1:numel(sdTerms)
        hasSdTerm = hasSdTerm || contains(baseLower, char(sdTerms(ti)));
    end
    isTiff = any(strcmpi(ext, {'.tif', '.tiff'}));
    if ~(hasToken && isTiff)
        continue;
    end

    score = 100 * double(hasToken) + 20 * double(hasSdTerm) - 0.001 * strlength(files(fi));
    if hasSdTerm && score > bestScore
        bestScore = score;
        bestFile = files(fi);
    end
    tokenOnlyFiles(end+1, 1) = files(fi); %#ok<AGROW>
    tokenOnlyScores(end+1, 1) = score; %#ok<AGROW>
end

sdFile = bestFile;
if strlength(sdFile) < 1 && ~isempty(tokenOnlyFiles)
    if allowTokenOnly
        [~, bestIdx] = max(tokenOnlyScores);
        sdFile = tokenOnlyFiles(bestIdx);
        warning(['Using token-only TIFF match for %s because opts.allowTokenOnlyTiffFallback=true: %s\n', ...
            'Confirm this is the SD TIFF, not a median/minimum image.'], caseName, sdFile);
    else
        fprintf(['Case %s: found %d TIFF(s) matching token "%s", but none contained ', ...
            'SD/std/deviation in the file name.\n'], caseName, numel(tokenOnlyFiles), searchToken);
        print_tiff_candidates(tokenOnlyFiles, 8);
    end
end
end

function files = collect_tiff_files(folders)
files = strings(0, 1);
for di = 1:numel(folders)
    files = [files; recursive_tiff_files(folders(di))]; %#ok<AGROW>
end
files = unique(files);
end

function print_tiff_candidates(files, maxToPrint)
if nargin < 2 || isempty(maxToPrint)
    maxToPrint = 8;
end
nPrint = min(numel(files), maxToPrint);
for i = 1:nPrint
    fprintf('  candidate: %s\n', files(i));
end
if numel(files) > nPrint
    fprintf('  ... %d more candidate(s)\n', numel(files) - nPrint);
end
end

function sdFile = prompt_for_sd_tiff_file(caseName, sdTiffDirs)
sdFile = "";
if exist('uigetfile', 'file') ~= 2
    return;
end

startDir = pwd;
if ~isempty(sdTiffDirs)
    startDir = char(sdTiffDirs(1));
end

try
    [fileName, pathName] = uigetfile( ...
        {'*.tif;*.tiff', 'TIFF files (*.tif, *.tiff)'; '*.*', 'All files'}, ...
        sprintf('Select 8-bit SD TIFF for %s', char(caseName)), ...
        startDir);
catch
    return;
end

if isequal(fileName, 0) || isequal(pathName, 0)
    return;
end

sdFile = string(fullfile(pathName, fileName));
fprintf('Selected SD TIFF for %s: %s\n', caseName, sdFile);
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

function maskOut = keep_largest_mask_component(maskIn)
maskIn = logical(maskIn);
maskOut = false(size(maskIn));
if ~any(maskIn(:))
    return;
end

if exist('bwconncomp', 'file') == 2
    cc = bwconncomp(maskIn, 8);
    if cc.NumObjects < 1
        return;
    end
    sizes = cellfun(@numel, cc.PixelIdxList);
    [~, bestIdx] = max(sizes);
    maskOut(cc.PixelIdxList{bestIdx}) = true;
    return;
end

maskOut = largest_component_fallback(maskIn);
end

function bestMask = largest_component_fallback(maskIn)
[nRows, nCols] = size(maskIn);
visited = false(nRows, nCols);
bestIdx = [];

for r0 = 1:nRows
    for c0 = 1:nCols
        if ~maskIn(r0, c0) || visited(r0, c0)
            continue;
        end

        queueR = zeros(nnz(maskIn), 1);
        queueC = zeros(nnz(maskIn), 1);
        componentIdx = zeros(nnz(maskIn), 1);
        head = 1;
        tail = 1;
        nComp = 0;
        queueR(tail) = r0;
        queueC(tail) = c0;
        visited(r0, c0) = true;

        while head <= tail
            r = queueR(head);
            c = queueC(head);
            head = head + 1;
            nComp = nComp + 1;
            componentIdx(nComp) = sub2ind([nRows nCols], r, c);

            for dr = -1:1
                for dc = -1:1
                    if dr == 0 && dc == 0
                        continue;
                    end
                    rr = r + dr;
                    cc = c + dc;
                    if rr < 1 || rr > nRows || cc < 1 || cc > nCols
                        continue;
                    end
                    if maskIn(rr, cc) && ~visited(rr, cc)
                        tail = tail + 1;
                        queueR(tail) = rr;
                        queueC(tail) = cc;
                        visited(rr, cc) = true;
                    end
                end
            end
        end

        componentIdx = componentIdx(1:nComp);
        if numel(componentIdx) > numel(bestIdx)
            bestIdx = componentIdx;
        end
    end
end

bestMask = false(size(maskIn));
bestMask(bestIdx) = true;
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

function s = std_omitnan_local(x)
x = x(isfinite(x));
if numel(x) < 2
    s = NaN;
else
    s = std(x, 0);
end
end

function p = one_sided_upper_tail_p(nullVals, observedVal)
nullVals = nullVals(isfinite(nullVals));
if isempty(nullVals) || ~isfinite(observedVal)
    p = NaN;
else
    p = (1 + nnz(nullVals >= observedVal)) / (numel(nullVals) + 1);
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

function files = existing_files(candidates)
files = strings(0, 1);
for i = 1:numel(candidates)
    if isfile(candidates(i))
        files(end+1, 1) = string(candidates(i)); %#ok<AGROW>
    end
end
files = unique(files);
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

function s = append_struct_local(s, value)
if isempty(s)
    s = value;
else
    s(end+1, 1) = value; %#ok<AGROW>
end
end

function seed = stable_string_seed(value)
chars = char(string(value));
if isempty(chars)
    seed = 0;
else
    weights = 1:numel(chars);
    seed = mod(sum(double(chars) .* weights), 100000);
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
