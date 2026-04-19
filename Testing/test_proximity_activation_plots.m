%% test_proximity_activation_plots.m
% Standalone test for proximity-activation publication plots.
% Preferred use: point matDir or matFile to a completed main run and
% regenerate plots from proximity_activation_by_case.mat. If no saved .mat
% file is found, the script falls back to synthetic smoke-test data.

clear; clc;

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot);

%% Paths
% Option 1: point this to your plot_data_mat folder from a completed run.
% Example:
% matDir = "E:\March Re 90,000 inception data\Processed images\results\results 33 local\plot_data_mat";
matDir = "";

% Option 2: or point directly to the .mat file. If this is non-empty it
% takes priority over matDir.
% Valid files are either:
%   resultsDir\plot_data_mat\proximity_activation_by_case.mat
%   resultsDir\Figures_PNG_SVG\Proximity activation\proximity_activation_by_case.mat
matFile = "";

outDir = fullfile(fileparts(mfilename('fullpath')), 'test_outputs', 'ProximityActivation');
if ~isfolder(outDir), mkdir(outDir); end

set(0, 'DefaultFigureVisible', 'on');
plotOpts = struct();
plotOpts.enableNormalTheme = true;
plotOpts.enablePosterTheme = false;
plotOpts.savePNG = true;
plotOpts.saveSVG = false;
plotOpts.themes = "normal";
plotOpts.keepFiguresOpen = true;

%% Load saved result .mat if available; otherwise build synthetic demo data
[allProximityActivation, opts, dataLabel] = load_or_make_proximity_data(matDir, matFile);

fprintf('Using %s proximity activation data.\n', dataLabel);
if ~isempty(allProximityActivation)
    fprintf('Loaded %d case(s).\n', numel(allProximityActivation));
    for i = 1:numel(allProximityActivation)
        nPairs = 0;
        nEvents = 0;
        if isfield(allProximityActivation(i), 'pairTable') && istable(allProximityActivation(i).pairTable)
            nPairs = height(allProximityActivation(i).pairTable);
        end
        if isfield(allProximityActivation(i), 'eventTable') && istable(allProximityActivation(i).eventTable)
            nEvents = height(allProximityActivation(i).eventTable);
        end
        fprintf('  %s: Re=%g, k/d=%.4g, primary events=%d, neighbor pairs=%d\n', ...
            string(allProximityActivation(i).caseName), allProximityActivation(i).Re, ...
            allProximityActivation(i).kD, nEvents, nPairs);
    end
end

[pairTable, eventTable, binnedStats, summaryTable] = proximity_activation_to_tables(allProximityActivation);
write_table_csv_compat(pairTable, fullfile(outDir, 'proximity_activation_pairs_test.csv'));
write_table_csv_compat(eventTable, fullfile(outDir, 'proximity_activation_events_test.csv'));
write_table_csv_compat(binnedStats, fullfile(outDir, 'proximity_activation_binned_test.csv'));
write_table_csv_compat(summaryTable, fullfile(outDir, 'proximity_activation_summary_test.csv'));
save(fullfile(outDir, 'proximity_activation_test_input.mat'), 'allProximityActivation', 'opts', 'dataLabel');

fprintf('Generating proximity activation plots...\n');
plot_proximity_activation_analysis(allProximityActivation, outDir, plotOpts, opts);
fprintf('Done. Output in: %s\n', outDir);


% =========================================================================
function [allProximityActivation, opts, dataLabel] = load_or_make_proximity_data(matDir, matFile)
matDir = string(matDir);
matFile = string(matFile);

if strlength(strtrim(matFile)) == 0 && strlength(strtrim(matDir)) > 0
    matFile = fullfile(matDir, "proximity_activation_by_case.mat");
end

if strlength(strtrim(matFile)) > 0 && isfile(matFile)
    S = load(matFile);
    if ~isfield(S, 'allProximityActivation')
        error('MAT file does not contain allProximityActivation: %s', matFile);
    end
    allProximityActivation = S.allProximityActivation;
    if isfield(S, 'proximityActivationOpts')
        opts = S.proximityActivationOpts;
    elseif isfield(S, 'opts')
        opts = S.opts;
    else
        opts = default_plot_opts();
    end
    opts = complete_plot_opts(opts);
    dataLabel = "saved .mat";
    fprintf('Loaded: %s\n', matFile);
    return;
elseif strlength(strtrim(matFile)) > 0
    warning('Requested proximity .mat file not found: %s\nFalling back to synthetic data.', matFile);
end

[allProximityActivation, opts] = make_synthetic_dataset();
dataLabel = "synthetic fallback";
end


% =========================================================================
function opts = default_plot_opts()
opts = struct();
opts.maxGamma = 20;
opts.extendedMaxGamma = 40;
opts.gammaBins = [0 2 5 10 20 40];
opts.minPairsForTrend = 2;
end


% =========================================================================
function opts = complete_plot_opts(opts)
defaults = default_plot_opts();
names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(opts, names{i}) || isempty(opts.(names{i}))
        opts.(names{i}) = defaults.(names{i});
    end
end
end


% =========================================================================
function [allProximityActivation, opts] = make_synthetic_dataset()
opts = default_plot_opts();

rng(11, 'twister');
caseSpecs = { ...
    "synthetic_5um",  95000, 0.030, 0.050
    "synthetic_30um", 95000, 0.267, 0.075
    "synthetic_5um",  124000, 0.030, 0.045
    "synthetic_30um", 124000, 0.267, 0.065};

allProximityActivation = repmat(make_case_template(), size(caseSpecs, 1), 1);
for i = 1:size(caseSpecs, 1)
    allProximityActivation(i) = make_synthetic_case(caseSpecs{i,1}, caseSpecs{i,2}, ...
        caseSpecs{i,3}, caseSpecs{i,4}, opts);
end
end


% =========================================================================
function d = make_synthetic_case(caseName, Re, kD, baseResponse, opts)
d = make_case_template();
d.caseName = string(caseName);
d.Re = Re;
d.kD = kD;
d.pixelSize_mm = 0.00375;
d.opts = opts;

nEvents = 35;
nPairs = 260;
gamma = min(opts.extendedMaxGamma, exprnd_fallback(9, nPairs, 1));
gamma = max(0.35, gamma);
secondaryProb = 0.04 * exp(-gamma / 5) .* (1 + 0.4 * (kD > 0.1));
isSecondary = rand(nPairs, 1) < secondaryProb;
response = abs(baseResponse + 0.12 * exp(-gamma / 4) .* randn(nPairs, 1) + 0.035 * randn(nPairs, 1));
response(isSecondary) = response(isSecondary) + 0.18 + 0.08 * rand(sum(isSecondary), 1);
randomResponse = abs(baseResponse * 0.55 + 0.025 * randn(nPairs, 1));
rdotRatio = exp(log(0.07 + 0.28 * exp(-gamma / 5)) + 0.45 * randn(nPairs, 1));

primaryIds = randi(nEvents, nPairs, 1);
neighborIds = (1:nPairs).' + 10000;
primaryActFrame = 1000 + primaryIds * 20;
primaryPeakFrame = primaryActFrame + 6;
primaryCollapseFrame = primaryActFrame + 14;
neighborActFrame = nan(nPairs, 1);
neighborActFrame(isSecondary) = primaryActFrame(isSecondary) + randi([1 19], sum(isSecondary), 1);

d.pairTable = table( ...
    repmat(string(caseName), nPairs, 1), repmat(Re, nPairs, 1), repmat(kD, nPairs, 1), ...
    ones(nPairs,1), repmat("synthetic.xml", nPairs, 1), ...
    primaryIds, primaryActFrame, primaryPeakFrame, primaryCollapseFrame, ...
    repmat(0.012, nPairs, 1), repmat(0.075, nPairs, 1), repmat(1250, nPairs, 1), ...
    neighborIds, repmat(0.011, nPairs, 1), 0.65 + 0.7 * rand(nPairs, 1), ...
    gamma .* 0.075 .* (1.05 + 0.1 * randn(nPairs, 1)), ...
    gamma .* 0.075, gamma .* 0.075 .* (0.95 + 0.1 * randn(nPairs, 1)), ...
    gamma .* 0.075 .* (0.85 + 0.08 * randn(nPairs, 1)), ...
    gamma .* (1.05 + 0.10 * randn(nPairs, 1)), gamma, gamma .* (0.95 + 0.1 * randn(nPairs, 1)), ...
    gamma .* (0.85 + 0.08 * randn(nPairs, 1)), gamma, ...
    response, response .* (0.7 + 0.4 * rand(nPairs, 1)), -response .* (0.2 + 0.3 * rand(nPairs, 1)), ...
    response .* (0.15 * randn(nPairs, 1)), randomResponse, ...
    isSecondary, neighborActFrame, isSecondary, isSecondary & rand(nPairs,1) > 0.25, isSecondary & rand(nPairs,1) > 0.72, ...
    neighborActFrame - primaryActFrame, neighborActFrame - primaryPeakFrame, neighborActFrame - primaryCollapseFrame, ...
    0.08 * randn(nPairs, 1), 0.15 * randn(nPairs, 1), randi([-5 5], nPairs, 1), rdotRatio, ...
    randi([5 18], nPairs, 1), randi([4 15], nPairs, 1), ...
    'VariableNames', { ...
    'Case','Re','kD','XMLSample','XMLFile', ...
    'primaryTrackId','primaryActivationFrame','primaryPeakFrame','primaryCollapseFrame', ...
    'primaryR0_mm','primaryRmax_mm','primaryPeakArea_px2', ...
    'neighborTrackId','neighborR0_mm','sizeRatioNeighborToPrimary', ...
    'distanceAtActivation_mm','distanceAtPeak_mm','distanceAtCollapse_mm','minDistanceDuringEvent_mm', ...
    'gammaAtActivation','gammaAtPeak','gammaAtCollapse','gammaMinDuringEvent','gammaForPlot', ...
    'maxAbsDeltaR_over_R0','maxExpansion_over_R0','maxCompression_over_R0','endDeltaR_over_R0','randomMaxAbsDeltaR_over_R0', ...
    'neighborActivated','neighborActivationFrame','neighborSecondaryActivated','neighborActivatedAfterPrimaryPeak','neighborActivatedAfterPrimaryCollapse', ...
    'activationLagFromPrimaryActivation_frames','activationLagFromPrimaryPeak_frames','activationLagFromPrimaryCollapse_frames', ...
    'rhoZeroLag','rhoMaxLag','lagAtMax_frames','rdotRmsRatio','nCorrelationFrames','nResponseFrames'});

d.eventTable = synthetic_event_table(caseName, Re, kD, nEvents, d.pairTable);
d.binnedStats = synthetic_binned_table(caseName, Re, kD, d.pairTable, opts);
d.summaryTable = table(string(caseName), Re, kD, nEvents, sum(d.eventTable.nViableNeighbors > 0), ...
    nPairs, sum(isSecondary), sum(isSecondary) / nPairs, ...
    sum(d.pairTable.gammaMinDuringEvent <= 10), ...
    sum(d.pairTable.gammaMinDuringEvent <= 10 & d.pairTable.neighborSecondaryActivated), ...
    sum(d.pairTable.gammaMinDuringEvent <= 10 & d.pairTable.neighborSecondaryActivated) / max(sum(d.pairTable.gammaMinDuringEvent <= 10), 1), ...
    finite_median(d.eventTable.nearestStableGamma), median(response), median(randomResponse), median(rdotRatio), ...
    opts.maxGamma, opts.extendedMaxGamma, 5, ...
    'VariableNames', {'Case','Re','kD','nPrimaryEvents','nPrimaryEventsWithNeighbors', ...
    'nNeighborPairs','nSecondaryActivations','secondaryActivationProbability', ...
    'nPairsWithin10Rmax','nSecondaryWithin10Rmax','secondaryActivationProbabilityWithin10Rmax', ...
    'medianNearestStableGamma','medianMaxAbsDeltaR_over_R0','medianRandomMaxAbsDeltaR_over_R0','medianRdotRmsRatio', ...
    'maxGamma','extendedMaxGamma','secondaryPostCollapseFrames'});
end


% =========================================================================
function T = synthetic_event_table(caseName, Re, kD, nEvents, pairTable)
rows = cell(nEvents, 22);
for i = 1:nEvents
    mask = pairTable.primaryTrackId == i;
    gamma = pairTable.gammaMinDuringEvent(mask);
    secondary = pairTable.neighborSecondaryActivated(mask);
    stableGamma = gamma(~secondary);
    activatedGamma = gamma(secondary);
    rows(i,:) = { ...
        string(caseName), Re, kD, 1, "synthetic.xml", ...
        i, 1000 + i*20, 1000 + i*20 + 6, 1000 + i*20 + 14, ...
        0.012, 0.075, 1250, ...
        sum(mask), sum(gamma <= 5), sum(gamma <= 10), sum(gamma <= 20), ...
        sum(gamma <= 5 & secondary), sum(gamma <= 10 & secondary), sum(gamma <= 20 & secondary), ...
        finite_min(gamma), finite_min(stableGamma), finite_min(activatedGamma)};
end
T = cell2table(rows, 'VariableNames', { ...
    'Case','Re','kD','XMLSample','XMLFile', ...
    'primaryTrackId','primaryActivationFrame','primaryPeakFrame','primaryCollapseFrame', ...
    'primaryR0_mm','primaryRmax_mm','primaryPeakArea_px2', ...
    'nViableNeighbors','nViableWithin5Rmax','nViableWithin10Rmax','nViableWithin20Rmax', ...
    'nSecondaryActivatedWithin5Rmax','nSecondaryActivatedWithin10Rmax','nSecondaryActivatedWithin20Rmax', ...
    'nearestViableGamma','nearestStableGamma','nearestActivatedGamma'});
end


% =========================================================================
function T = synthetic_binned_table(caseName, Re, kD, pairTable, opts)
rows = {};
for bi = 1:(numel(opts.gammaBins) - 1)
    lo = opts.gammaBins(bi);
    hi = opts.gammaBins(bi + 1);
    gamma = pairTable.gammaMinDuringEvent;
    mask = gamma >= lo & gamma < hi;
    if bi == (numel(opts.gammaBins) - 1)
        mask = gamma >= lo & gamma <= hi;
    end
    nPairs = sum(mask);
    if nPairs == 0
        continue;
    end
    nSecondary = sum(pairTable.neighborSecondaryActivated(mask));
    [ciLow, ciHigh] = wilson_ci_local(nSecondary, nPairs);
    resp = pairTable.maxAbsDeltaR_over_R0(mask);
    randomResp = pairTable.randomMaxAbsDeltaR_over_R0(mask);
    rdotRatio = pairTable.rdotRmsRatio(mask);
    rows(end+1,:) = {string(caseName), Re, kD, lo, hi, (lo+hi)/2, nPairs, nSecondary, ...
        nSecondary/nPairs, ciLow, ciHigh, median(resp), prctile_fallback(resp, 90), ...
        median(randomResp), prctile_fallback(randomResp, 90), median(rdotRatio), prctile_fallback(rdotRatio, 90)}; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', { ...
    'Case','Re','kD','gammaBinLow','gammaBinHigh','gammaBinCenter', ...
    'nPairs','nSecondaryActivated','secondaryActivationProbability','secondaryActivationCI_low','secondaryActivationCI_high', ...
    'medianMaxAbsDeltaR_over_R0','p90MaxAbsDeltaR_over_R0', ...
    'medianRandomMaxAbsDeltaR_over_R0','p90RandomMaxAbsDeltaR_over_R0', ...
    'medianRdotRmsRatio','p90RdotRmsRatio'});
end


% =========================================================================
function tpl = make_case_template()
tpl = struct('caseName',"", 'Re',NaN, 'kD',NaN, 'pixelSize_mm',NaN, ...
    'pairTable',table(), 'eventTable',table(), 'binnedStats',table(), ...
    'summaryTable',table(), 'skipCounts',struct(), 'opts',struct());
end


% =========================================================================
function x = exprnd_fallback(mu, m, n)
x = -mu .* log(max(rand(m, n), realmin));
end


% =========================================================================
function v = finite_min(x)
x = x(isfinite(x(:)));
if isempty(x)
    v = NaN;
else
    v = min(x);
end
end


% =========================================================================
function v = finite_median(x)
x = x(isfinite(x(:)));
if isempty(x)
    v = NaN;
else
    v = median(x);
end
end


% =========================================================================
function q = prctile_fallback(x, p)
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
function [ciLow, ciHigh] = wilson_ci_local(k, n)
z = 1.95996398454005;
phat = k / n;
denom = 1 + z^2 / n;
center = (phat + z^2 / (2 * n)) / denom;
halfWidth = z * sqrt((phat * (1 - phat) / n) + (z^2 / (4 * n^2))) / denom;
ciLow = max(0, center - halfWidth);
ciHigh = min(1, center + halfWidth);
end
