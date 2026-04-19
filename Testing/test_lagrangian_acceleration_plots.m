%% test_lagrangian_acceleration_plots.m
% Standalone synthetic test for Lagrangian acceleration publication plots.
% This does not require TrackMate XML data; it exercises the plotting and
% summary-table paths with realistic-shaped synthetic distributions.

clear; clc;

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot);

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

lagAccelOpts = struct();
lagAccelOpts.throatHeight_mm = 10;
lagAccelOpts.xLimNorm = [0 0.5];
lagAccelOpts.yLimNorm = [0 0.12];
lagAccelOpts.heatmapGridSize = [25 25];
lagAccelOpts.heatmapStats = ["median", "p90"];
lagAccelOpts.heatmapMinSamplesPerBin = 2;
lagAccelOpts.activationOverlayPerCase = 75;
lagAccelOpts.randomSeed = 42;
lagAccelOpts.minSpearmanN = 8;

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

summaryTable = lagrangian_acceleration_to_table(allLagAccel);
write_table_csv_compat(summaryTable, fullfile(outDir, 'lagrangian_acceleration_summary_synthetic.csv'));
save(fullfile(outDir, 'lagrangian_acceleration_synthetic.mat'), 'allLagAccel', 'lagAccelOpts');

fprintf('Generating synthetic Lagrangian acceleration plots...\n');
plot_lagrangian_acceleration_analysis(allLagAccel, outDir, plotOpts, lagAccelOpts);
fprintf('Done. Output in: %s\n', outDir);


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
summary.triggerWindowFrames = 5;
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
