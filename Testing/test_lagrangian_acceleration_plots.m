%% test_lagrangian_acceleration_plots.m
% Standalone test for Lagrangian acceleration publication plots.
% Preferred use: point matDir or matFile to a completed main run and
% regenerate plots from lagrangian_acceleration_by_case.mat. If no saved
% .mat file is found, the script falls back to synthetic smoke-test data.

clear; clc;

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot);
cbrewer2Dir = fullfile(repoRoot, 'external', 'cbrewer2', 'cbrewer2');
if isfolder(cbrewer2Dir)
    addpath(cbrewer2Dir);
end

%% Paths
% Option 1: point this to your plot_data_mat folder from a completed run.
% Example:
% matDir = "E:\March Re 90,000 inception data\Processed images\results\results 33 local\plot_data_mat";
matDir = "C:\Users\kbsanjayvasanth\Downloads\plot_data_mat (1)";

% Option 2: or point directly to the .mat file. If this is non-empty it
% takes priority over matDir.
% Valid files are either:
%   resultsDir\plot_data_mat\lagrangian_acceleration_by_case.mat
%   resultsDir\Figures_PNG_SVG\lagrangian acceleration\lagrangian_acceleration_by_case.mat
matFile = "";

% Optional: point this to the ROI_throat.mat used by the main batch run so
% the heatmap test also shows the fixed-color wall overlay.
% If this file is missing, the test creates a synthetic wall mask so the
% plotting path is still exercised.
roiFile = "E:\March Re 90,000 inception data\Processed images\results\2000 frames\Smooth variation 2000\ROI_throat.mat";

outDir = fullfile(fileparts(mfilename('fullpath')), 'test_outputs', 'LagrangianAcceleration');
if ~isfolder(outDir), mkdir(outDir); end

set(0, 'DefaultFigureVisible', 'on');
plotOpts = struct();
plotOpts.enableNormalTheme = true;
plotOpts.enablePosterTheme = false;
plotOpts.savePNG = true;
plotOpts.saveSVG = false;
plotOpts.themes = "normal";
plotOpts.keepFiguresOpen = true;
plotOpts.inceptionImageSize_px = [1280 320];
plotOpts.lagrangianHeatmapImageSize_px = [1280 512];

%% Plot selection
% Set any of these to false when you only want to regenerate specific plots.
makeAllFramePdfPlot = false;
makeTriggerWindowPdfPlot = false;
makePeakTriggerGrowthPlot = false;
makeHeatmapPlots = true;
makeSanityCheckPlots = false;

%% Heatmap appearance overrides
% Edit these top-level values for testing. They are applied after loading
% the saved .mat, so they override any lagAccelOpts stored in the run file.
% With square bins on, the 2nd value controls square size:
%   larger 2nd value = more y-bins = smaller square cells
%   smaller 2nd value = fewer y-bins = larger square cells
heatmapGridSizeOverride = [12 12];
heatmapStatsOverride = "mean";
makePooledHeatmapPlotsOverride = false;
makeCaseHeatmapPlotsOverride = true;
heatmapMinSamplesPerBinOverride = 8;
heatmapSquareBinsOverride = true;
heatmapXLimNormOverride = [0 0.48];
heatmapShowActivationContoursOverride = false;
heatmapShowAccelerationRidgeOverride = true;
heatmapAccelerationRidgeStatsOverride = "mean";
heatmapAccelerationRidgePercentileOverride = 0;
heatmapAccelerationRidgeSmoothWindowOverride = 9;
heatmapColormapOverride = "cbrewer2:YlGnBu";
heatmapColormapByStatOverride = struct( ...
    'mean', "cbrewer2:GnBu", ...
    'median', "cbrewer2:PuBu", ...
    'p90', "cbrewer2:Blues");
activationOverlayEdgeMarginNormOverride = 0;
heatmapColorPercentileRangeOverride = [5 95];
heatmapActivationMarkerSizeOverride = 23;
activationOverlayPerCaseOverride = 250; % set Inf to plot every available non-wall activation marker
heatmapLegendLocationOverride = 'northwest';
heatmapLegendFontSizeOverride = 10;
triggerPdfLegendLocationOverride = 'southeast';
makeTriggerPdfZoomPlotOverride = true;
triggerPdfZoomXLimOverride = [1e-1 1];
makeTriggerSurvivalPlotOverride = false;
triggerSurvivalXLimOverride = [1e-1 1];
triggerSurvivalLegendLocationOverride = 'southwest';
showTriggerSurvivalLegendOverride = false;
makeTriggerThresholdSummaryPlotOverride = false;
triggerSurvivalThresholdsOverride = [0.3 0.5 0.75];
allFramePdfLegendLocationOverride = 'southwest';
peakGrowthLegendLocationOverride = 'northwest';
peakGrowthLegendFontSizeOverride = 10;

%% Load saved result .mat if available; otherwise build synthetic demo data
[allLagAccel, lagAccelOpts, dataLabel] = load_or_make_lagrangian_data(matDir, matFile);
lagAccelOpts.makeAllFramePdfPlot = makeAllFramePdfPlot;
lagAccelOpts.makeTriggerWindowPdfPlot = makeTriggerWindowPdfPlot;
lagAccelOpts.makePeakTriggerGrowthPlot = makePeakTriggerGrowthPlot;
lagAccelOpts.makeHeatmapPlots = makeHeatmapPlots;
lagAccelOpts.makeSanityCheckPlots = makeSanityCheckPlots;
lagAccelOpts.makePooledHeatmapPlots = makePooledHeatmapPlotsOverride;
lagAccelOpts.makeCaseHeatmapPlots = makeCaseHeatmapPlotsOverride;
lagAccelOpts.xLimNorm = heatmapXLimNormOverride;
lagAccelOpts.heatmapGridSize = heatmapGridSizeOverride;
lagAccelOpts.heatmapStats = heatmapStatsOverride;
lagAccelOpts.heatmapMinSamplesPerBin = heatmapMinSamplesPerBinOverride;
lagAccelOpts.heatmapSquareBins = heatmapSquareBinsOverride;
lagAccelOpts.heatmapShowActivationContours = heatmapShowActivationContoursOverride;
lagAccelOpts.heatmapShowAccelerationRidge = heatmapShowAccelerationRidgeOverride;
lagAccelOpts.heatmapAccelerationRidgeStats = heatmapAccelerationRidgeStatsOverride;
lagAccelOpts.heatmapAccelerationRidgePercentile = heatmapAccelerationRidgePercentileOverride;
lagAccelOpts.heatmapAccelerationRidgeSmoothWindow = heatmapAccelerationRidgeSmoothWindowOverride;
lagAccelOpts.heatmapColormap = heatmapColormapOverride;
lagAccelOpts.heatmapColormapByStat = heatmapColormapByStatOverride;
lagAccelOpts.activationOverlayEdgeMarginNorm = activationOverlayEdgeMarginNormOverride;
lagAccelOpts.heatmapColorPercentileRange = heatmapColorPercentileRangeOverride;
lagAccelOpts.heatmapActivationMarkerSize = heatmapActivationMarkerSizeOverride;
lagAccelOpts.activationOverlayPerCase = activationOverlayPerCaseOverride;
lagAccelOpts.heatmapLegendLocation = heatmapLegendLocationOverride;
lagAccelOpts.heatmapLegendFontSize = heatmapLegendFontSizeOverride;
lagAccelOpts.triggerPdfLegendLocation = triggerPdfLegendLocationOverride;
lagAccelOpts.makeTriggerPdfZoomPlot = makeTriggerPdfZoomPlotOverride;
lagAccelOpts.triggerPdfZoomXLim = triggerPdfZoomXLimOverride;
lagAccelOpts.makeTriggerSurvivalPlot = makeTriggerSurvivalPlotOverride;
lagAccelOpts.triggerSurvivalXLim = triggerSurvivalXLimOverride;
lagAccelOpts.triggerSurvivalLegendLocation = triggerSurvivalLegendLocationOverride;
lagAccelOpts.showTriggerSurvivalLegend = showTriggerSurvivalLegendOverride;
lagAccelOpts.makeTriggerThresholdSummaryPlot = makeTriggerThresholdSummaryPlotOverride;
lagAccelOpts.triggerSurvivalThresholds = triggerSurvivalThresholdsOverride;
lagAccelOpts.allFramePdfLegendLocation = allFramePdfLegendLocationOverride;
lagAccelOpts.peakGrowthLegendLocation = peakGrowthLegendLocationOverride;
lagAccelOpts.peakGrowthLegendFontSize = peakGrowthLegendFontSizeOverride;
plotOpts.roiData = load_or_make_lagrangian_roi(roiFile, allLagAccel, lagAccelOpts);

fprintf('Using %s Lagrangian acceleration data.\n', dataLabel);
if ~isempty(allLagAccel)
    fprintf('Loaded %d case(s).\n', numel(allLagAccel));
    for i = 1:numel(allLagAccel)
        nSamples = 0;
        nTrigger = 0;
        if isfield(allLagAccel(i), 'sampleAstar')
            nSamples = numel(allLagAccel(i).sampleAstar);
        end
        if isfield(allLagAccel(i), 'activatedTriggerAstar')
            nTrigger = numel(allLagAccel(i).activatedTriggerAstar);
        end
        fprintf('  %s: Re=%g, k/d=%.4g, all-frame samples=%d, trigger samples=%d\n', ...
            string(allLagAccel(i).caseName), allLagAccel(i).Re, allLagAccel(i).kD, nSamples, nTrigger);
    end
end

summaryTable = lagrangian_acceleration_to_table(allLagAccel);
write_table_csv_compat(summaryTable, fullfile(outDir, 'lagrangian_acceleration_summary_test.csv'));
save(fullfile(outDir, 'lagrangian_acceleration_test_input.mat'), 'allLagAccel', 'lagAccelOpts', 'dataLabel');

fprintf('Generating Lagrangian acceleration plots...\n');
plot_lagrangian_acceleration_analysis(allLagAccel, outDir, plotOpts, lagAccelOpts);
fprintf('Done. Output in: %s\n', outDir);


% =========================================================================
function [allLagAccel, lagAccelOpts, dataLabel] = load_or_make_lagrangian_data(matDir, matFile)
matDir = string(matDir);
matFile = string(matFile);

if strlength(strtrim(matFile)) == 0 && strlength(strtrim(matDir)) > 0
    matFile = fullfile(matDir, "lagrangian_acceleration_by_case.mat");
end

if strlength(strtrim(matFile)) > 0 && isfile(matFile)
    S = load(matFile);
    if ~isfield(S, 'allLagAccel')
        error('MAT file does not contain allLagAccel: %s', matFile);
    end
    allLagAccel = S.allLagAccel;
    if isfield(S, 'lagAccelOpts')
        lagAccelOpts = S.lagAccelOpts;
    else
        lagAccelOpts = default_lag_accel_plot_opts();
    end
    lagAccelOpts = complete_lag_accel_plot_opts(lagAccelOpts);
    dataLabel = "saved .mat";
    fprintf('Loaded: %s\n', matFile);
    return;
elseif strlength(strtrim(matFile)) > 0
    warning('Requested Lagrangian acceleration .mat file not found: %s\nFalling back to synthetic data.', matFile);
end

[allLagAccel, lagAccelOpts] = make_synthetic_dataset();
dataLabel = "synthetic fallback";
end


% =========================================================================
function lagAccelOpts = default_lag_accel_plot_opts()
lagAccelOpts = struct();
lagAccelOpts.throatHeight_mm = 10;
lagAccelOpts.xLimNorm = [0 0.5];
lagAccelOpts.yLimNorm = [0 0.12];
lagAccelOpts.makeAllFramePdfPlot = true;
lagAccelOpts.makeTriggerWindowPdfPlot = true;
lagAccelOpts.makePeakTriggerGrowthPlot = true;
lagAccelOpts.makeHeatmapPlots = true;
lagAccelOpts.makeSanityCheckPlots = false;
lagAccelOpts.heatmapGridSize = [20 20]; % with square bins, 2nd value sets y-bin count; x-bin count is auto
lagAccelOpts.heatmapStats = ["median", "p90"];
lagAccelOpts.heatmapColormap = "cbrewer2:YlGnBu";
lagAccelOpts.heatmapPreserveSpatialAspect = true;
lagAccelOpts.heatmapSquareBins = true;
lagAccelOpts.heatmapMinSamplesPerBin = 2;
lagAccelOpts.activationOverlayPerCase = 250;
lagAccelOpts.randomSeed = 42;
lagAccelOpts.minSpearmanN = 8;
lagAccelOpts.maxSanityTracks = 250;
lagAccelOpts.maxSanityTracksToPlot = 250;
lagAccelOpts.imageSize_px = [1280 320];
end


% =========================================================================
function lagAccelOpts = complete_lag_accel_plot_opts(lagAccelOpts)
defaults = default_lag_accel_plot_opts();
names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(lagAccelOpts, names{i}) || isempty(lagAccelOpts.(names{i}))
        lagAccelOpts.(names{i}) = defaults.(names{i});
    end
end
end


% =========================================================================
function roiData = load_or_make_lagrangian_roi(roiFile, allLagAccel, lagAccelOpts)
roiFile = string(roiFile);
pixelSize = infer_lagrangian_pixel_size(allLagAccel);

if strlength(strtrim(roiFile)) > 0 && isfile(roiFile)
    R = load(roiFile);
    roiData = struct();
    roiData.wallMask = extract_roi_mask(R, 'wallMask');
    roiData.unwantedTrackMask = extract_roi_mask(R, 'unwantedTrackMask');
    roiData.throat_xy_px = extract_roi_throat(R);
    roiData.maskPixelSize = pixelSize;
    fprintf('Loaded ROI wall overlay for heatmap test: %s\n', roiFile);
    fprintf('  Wall mask: %d pixels\n', sum(roiData.wallMask(:)));
    return;
end

roiData = make_synthetic_lagrangian_roi(pixelSize, lagAccelOpts);
fprintf('ROI file not found for heatmap test; using synthetic wall overlay.\n');
end


% =========================================================================
function pixelSize = infer_lagrangian_pixel_size(allLagAccel)
pixelSize = 0.00375009375;
if isempty(allLagAccel)
    return;
end

vals = nan(0, 1);
for i = 1:numel(allLagAccel)
    if isfield(allLagAccel(i), 'pixelSize_mm') && isfinite(allLagAccel(i).pixelSize_mm) && allLagAccel(i).pixelSize_mm > 0
        vals(end+1, 1) = allLagAccel(i).pixelSize_mm; %#ok<AGROW>
    end
end
if ~isempty(vals)
    pixelSize = median(vals);
end
end


% =========================================================================
function mask = extract_roi_mask(R, fieldName)
mask = [];
if isfield(R, fieldName)
    mask = R.(fieldName);
elseif isfield(R, 'ROI_throat') && isfield(R.ROI_throat, fieldName)
    mask = R.ROI_throat.(fieldName);
end

if isempty(mask)
    mask = false(320, 1400);
end
mask = logical(mask);
end


% =========================================================================
function throat_xy_px = extract_roi_throat(R)
throat_xy_px = [NaN NaN];
if isfield(R, 'ROI_throat') && isfield(R.ROI_throat, 'throat_xy_px')
    throat_xy_px = R.ROI_throat.throat_xy_px;
elseif isfield(R, 'x_throat') && isfield(R, 'y_throat')
    throat_xy_px = [R.x_throat R.y_throat];
end
end


% =========================================================================
function roiData = make_synthetic_lagrangian_roi(pixelSize, lagAccelOpts)
if isfield(lagAccelOpts, 'imageSize_px') && numel(lagAccelOpts.imageSize_px) >= 2
    nRows = max(32, round(lagAccelOpts.imageSize_px(2)));
    nCols = max(32, round(lagAccelOpts.imageSize_px(1)));
else
    nRows = 320;
    nCols = 1400;
end

wallMask = false(nRows, nCols);
yExtent_mm = nRows * pixelSize;
x = linspace(0, 1, nCols);
surface_mm = 0.10 + 0.035 * sin(2*pi*x) + 0.018 * exp(-((x - 0.55) / 0.18).^2);
surface_mm = max(0.04, min(0.22, surface_mm));
topRows = round((yExtent_mm - surface_mm) ./ pixelSize);
topRows = max(1, min(nRows, topRows));
for c = 1:nCols
    wallMask(topRows(c):end, c) = true;
end

roiData = struct();
roiData.wallMask = wallMask;
roiData.unwantedTrackMask = false(nRows, nCols);
roiData.throat_xy_px = [NaN NaN];
roiData.maskPixelSize = pixelSize;
end


% =========================================================================
function [allLagAccel, lagAccelOpts] = make_synthetic_dataset()
lagAccelOpts = default_lag_accel_plot_opts();

rng(7, 'twister');
caseSpecs = { ...
    "synthetic_5um",  95000, 0.005, 0.55, 1.00
    "synthetic_15um", 95000, 0.015, 0.80, 1.22
    "synthetic_5um",  120000, 0.005, 0.70, 1.08
    "synthetic_15um", 120000, 0.015, 1.05, 1.35};

allLagAccel = repmat(make_synthetic_case_template(), size(caseSpecs, 1), 1);
for i = 1:size(caseSpecs, 1)
    allLagAccel(i) = make_synthetic_case(caseSpecs{i,1}, caseSpecs{i,2}, caseSpecs{i,3}, ...
        caseSpecs{i,4}, caseSpecs{i,5}, lagAccelOpts);
end
end


% =========================================================================
function d = make_synthetic_case(caseName, Re, kD, scaleFactor, growthFactor, opts)
d = make_synthetic_case_template();
d.caseName = string(caseName);
d.Re = Re;
d.kD = kD;
d.pixelSize_mm = 0.00375;
d.dMean_m = 90e-6 * growthFactor;
d.U_ref_m_s = 0.42 + 0.05 * growthFactor;
d.opts = struct('makeSanityPlots', true);

nAll = 2400;
nAct = 900;
nNon = 1500;
x = 0.5 * rand(nAll, 1);
y = 0.12 * rand(nAll, 1);
hotspot = exp(-((x - 0.20 - 0.035*growthFactor).^2 / 0.006 + (y - 0.055).^2 / 0.0006));
baseA = exp(log(0.08 * scaleFactor) + 0.75 * randn(nAll, 1));
d.sampleAstar = baseA .* (1 + 2.5 * hotspot);
d.sampleRawAstar = d.sampleAstar .* exp(0.65 * randn(nAll, 1));
d.sampleXNorm = x;
d.sampleYNorm = y;
d.sampleIsActivatedTrack = false(nAll, 1);
d.sampleIsActivatedTrack(1:nAct) = true;

d.activatedAllAstar = d.sampleAstar(1:nAct) .* exp(0.25 + 0.35 * randn(nAct, 1));
d.nonActivatedAllAstar = d.sampleAstar((nAct+1):end) .* exp(-0.15 + 0.30 * randn(nNon, 1));
d.activatedTriggerAstar = exp(log(0.18 * scaleFactor) + 0.60 * randn(420, 1));
d.nonActivatedRandomAstar = exp(log(0.075 * scaleFactor) + 0.55 * randn(420, 1));

nPairs = 80;
d.peakTriggerAstar = exp(log(0.18 * scaleFactor) + 0.70 * randn(nPairs, 1));
d.growthRatio = max(1.05, 1.15 + 1.2 * rescale01(log(d.peakTriggerAstar)) + 0.18 * randn(nPairs, 1));

nActPts = 180;
d.activationXYNorm = [ ...
    min(max(0.20 + 0.035 * growthFactor + 0.045 * randn(nActPts, 1), opts.xLimNorm(1)), opts.xLimNorm(2)), ...
    min(max(0.055 + 0.015 * randn(nActPts, 1), opts.yLimNorm(1)), opts.yLimNorm(2))];
d.activationXY_mm = [d.activationXYNorm(:,1) * opts.throatHeight_mm, ...
    opts.throatHeight_mm * (opts.yLimNorm(2) - d.activationXYNorm(:,2))];

d.stationaryAstar = exp(log(0.012) + 0.40 * randn(300, 1));
d.nStationaryTracks = 12;
d.sanityTracks = make_sanity_tracks(opts);
d.summary = make_summary(d);
end


% =========================================================================
function tpl = make_synthetic_case_template()
tpl = struct( ...
    'caseName', "", 'Re', NaN, 'kD', NaN, 'pixelSize_mm', NaN, ...
    'dMean_m', NaN, 'U_ref_m_s', NaN, ...
    'sampleAstar', nan(0,1), 'sampleRawAstar', nan(0,1), ...
    'sampleXNorm', nan(0,1), 'sampleYNorm', nan(0,1), ...
    'sampleIsActivatedTrack', false(0,1), ...
    'activatedAllAstar', nan(0,1), 'nonActivatedAllAstar', nan(0,1), ...
    'activatedTriggerAstar', nan(0,1), 'nonActivatedRandomAstar', nan(0,1), ...
    'peakTriggerAstar', nan(0,1), 'growthRatio', nan(0,1), ...
    'activationXY_mm', zeros(0,2), 'activationXYNorm', zeros(0,2), ...
    'stationaryAstar', nan(0,1), 'nStationaryTracks', 0, ...
    'tracks', struct([]), 'sanityTracks', make_sanity_tracks(struct('xLimNorm',[0 0.5],'yLimNorm',[0 0.12])), ...
    'skipCounts', struct(), 'summary', struct(), 'opts', struct('makeSanityPlots', true));
end


% =========================================================================
function tracks = make_sanity_tracks(opts)
if ~isfield(opts, 'xLimNorm'), opts.xLimNorm = [0 0.5]; end
if ~isfield(opts, 'yLimNorm'), opts.yLimNorm = [0 0.12]; end
tracks = repmat(struct('TRACK_ID', NaN, 'isActivated', false, ...
    'xRawNorm', nan(0,1), 'yRawNorm', nan(0,1), ...
    'xSmoothNorm', nan(0,1), 'ySmoothNorm', nan(0,1), ...
    'activationXYNorm', [NaN NaN]), 5, 1);
for i = 1:5
    t = linspace(0, 1, 40).';
    xSmooth = 0.45 - 0.25 * t + 0.01 * sin(2*pi*t + i);
    ySmooth = 0.025 + 0.015 * i + 0.004 * sin(3*pi*t + 0.4*i);
    tracks(i).TRACK_ID = i;
    tracks(i).isActivated = mod(i, 2) == 0;
    tracks(i).xSmoothNorm = min(max(xSmooth, opts.xLimNorm(1)), opts.xLimNorm(2));
    tracks(i).ySmoothNorm = min(max(ySmooth, opts.yLimNorm(1)), opts.yLimNorm(2));
    tracks(i).xRawNorm = tracks(i).xSmoothNorm + 0.004 * randn(size(t));
    tracks(i).yRawNorm = tracks(i).ySmoothNorm + 0.003 * randn(size(t));
    if tracks(i).isActivated
        tracks(i).activationXYNorm = [tracks(i).xSmoothNorm(27), tracks(i).ySmoothNorm(27)];
    end
end
end


% =========================================================================
function summary = make_summary(d)
summary = struct();
summary.Case = d.caseName;
summary.Re = d.Re;
summary.kD = d.kD;
summary.dMean_um = d.dMean_m * 1e6;
summary.U_ref_m_s = d.U_ref_m_s;
summary.nTracksUsable = 100;
summary.nActivatedTracksUsable = 42;
summary.nNonActivatedTracksUsable = 58;
summary.nAllAccelSamples = numel(d.sampleAstar);
summary.nActivatedAllSamples = numel(d.activatedAllAstar);
summary.nNonActivatedAllSamples = numel(d.nonActivatedAllAstar);
summary.nActivatedTriggerSamples = numel(d.activatedTriggerAstar);
summary.nNonActivatedRandomSamples = numel(d.nonActivatedRandomAstar);
summary.nPeakGrowthPairs = numel(d.peakTriggerAstar);
summary.nStationaryTracks = d.nStationaryTracks;
summary.nStationarySamples = numel(d.stationaryAstar);
summary.activatedAllMedianAstar = median(d.activatedAllAstar);
summary.activatedAllP90Astar = prctile_local(d.activatedAllAstar, 90);
summary.nonActivatedAllMedianAstar = median(d.nonActivatedAllAstar);
summary.nonActivatedAllP90Astar = prctile_local(d.nonActivatedAllAstar, 90);
summary.activatedTriggerMedianAstar = median(d.activatedTriggerAstar);
summary.activatedTriggerP90Astar = prctile_local(d.activatedTriggerAstar, 90);
summary.nonActivatedRandomMedianAstar = median(d.nonActivatedRandomAstar);
summary.nonActivatedRandomP90Astar = prctile_local(d.nonActivatedRandomAstar, 90);
summary.stationaryP95Astar = prctile_local(d.stationaryAstar, 95);
summary.sgWindowFrames = 7;
summary.sgPolyOrder = 3;
summary.triggerWindowFrames = 3;
summary.trackPopulation = "strictPrimary";
end


% =========================================================================
function y = rescale01(x)
x = x(:);
xMin = min(x);
xMax = max(x);
if xMax <= xMin
    y = zeros(size(x));
else
    y = (x - xMin) ./ (xMax - xMin);
end
end


% =========================================================================
function q = prctile_local(x, p)
x = sort(x(isfinite(x(:))));
n = numel(x);
if n == 0, q = NaN; return; end
idx = 1 + (n - 1) * p / 100;
i0 = floor(idx);
i1 = ceil(idx);
if i0 == i1
    q = x(i0);
else
    q = x(i0) + (idx - i0) * (x(i1) - x(i0));
end
end
