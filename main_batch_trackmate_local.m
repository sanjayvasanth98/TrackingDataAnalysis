%% main_batch_trackmate_arc.m
% Batch TrackMate XML analysis (ARC + local)
% - Computes framewise Injection (left-moving availability) and Activation/Growth events
% - Supports one or multiple Reynolds-number sets for the same roughness ladder
% - Plots A/I vs k/D_h, inception (2x growth) locations, Tau vs k/D_h,
%   and upstream moving microbubble size distributions
%
% Output:
%   resultsDir/
%     Summary_tracks.csv
%     fit_AI_vs_kDh_by_Re.txt
%     Figures_PNG_SVG/(normal|poster)/...

clear; clc;
clear functions;
rehash;
bootstrap_user_matlab_path();
rehash;

%% ---------------- USER SETTINGS ----------------
% Toggle run mode by commenting one line
runMode = "local"; % <---edit
% runMode = "arc"; % <---edit

isArc = strcmpi(runMode, "arc");
visMode = ["on", "off"];
set(0,'DefaultFigureVisible', visMode(1 + isArc));

% Local quick-test parsing limit (total parsed tracks per XML).
% Keep Inf to parse all tracks; this is separate from left-moving analysis cap below.
localParseTrackLimit = Inf; % <---edit: keep Inf to avoid truncating true left-moving population
maxTracksChoices = [localParseTrackLimit, Inf];
maxTracksToParse = maxTracksChoices(1 + isArc);

% Local quick-test analysis cap (left-moving tracks kept for diagnostics/inception).
% For 100 left-moving tracks, use: localParseTrackLimit = Inf and localEventLimit = 100. <---edit
localEventLimit = 100; % <---edit: max number of left-moving tracks to analyze per case
leftMovingTrackLimitChoices = [localEventLimit, Inf];
maxLeftMovingTracks = leftMovingTrackLimitChoices(1 + isArc);

% Case selection:
% - "all"            -> run every case below
% - 1                -> run one case by index
% - [1 3 6]          -> run multiple cases by index
% - "5um"            -> run one case by name
% - ["5um","30um"]   -> run multiple cases by name
caseSelection = 1; % <---edit

% Cache parsed outputs (.mat) to avoid re-parsing XML on reruns
useMatCache = true; % <---edit
forceReparse = false; % <---edit

% Where to save all results
resultsDir = "E:\March Re 90,000 inception data\Processed images\results\results 4"; % <---edit
figDir     = fullfile(resultsDir, "Figures_PNG_SVG");
if ~isfolder(resultsDir), mkdir(resultsDir); end
if ~isfolder(figDir), mkdir(figDir); end
cacheFile = fullfile(resultsDir, "all_cases_cache.mat");

cacheDB = struct('key', strings(0,1), 'out', {cell(0,1)});
cacheUpdated = false;
if useMatCache && ~forceReparse && isfile(cacheFile)
    S = load(cacheFile, 'cacheDB');
    if isfield(S, 'cacheDB') && isfield(S.cacheDB, 'key') && isfield(S.cacheDB, 'out')
        cacheDB = S.cacheDB;
        fprintf("Loaded cache DB: %s (%d entries)\n", cacheFile, numel(cacheDB.key));
    end
elseif useMatCache && ~forceReparse
    legacyCacheFile = fullfile(resultsDir, "cache_mat", "all_cases_cache.mat");
    if isfile(legacyCacheFile)
        S = load(legacyCacheFile, 'cacheDB');
        if isfield(S, 'cacheDB') && isfield(S.cacheDB, 'key') && isfield(S.cacheDB, 'out')
            cacheDB = S.cacheDB;
            fprintf("Loaded legacy cache DB: %s (%d entries)\n", legacyCacheFile, numel(cacheDB.key));
        end
    end
end

% Parser + metric policy (included in cache key).
parserOpts = struct();
parserOpts.parserVersion = 2;
parserOpts.parseTrackedSpotsOnly = true;
parserOpts.parseFilteredTracksOnly = true;

% Bulk flow in this dataset is left->right; recirculating injection is counterflow (right->left).
flowOpts = struct();
flowOpts.bulkDirection = "left_to_right"; % <---edit if dataset orientation changes
flowOpts.minNetDxCounterflow_mm = 0.08;
flowOpts.minNegativeStepFraction = 0.65;
flowOpts.maxPositiveStepFraction = 0.30;
flowOpts.requireRightOrigin = true;
flowOpts.rightOriginFrac = 0.60; % rightmost 40% source band

qcOpts = struct();
qcOpts.minTrackSpots = 5;
qcOpts.maxTrackGaps = 1;
qcOpts.maxLeftMovingTracks = maxLeftMovingTracks;
qcOpts.rejectSplitMergeComplex = true;
qcOpts.wallBandEnabled = false;
qcOpts.wallBandYLimits_mm = [];

activationOpts = struct();
activationOpts.areaJumpFactor = 2.0;
activationOpts.preWindowFrames = 2;
activationOpts.minPrePoints = 2;
activationOpts.requiredPostFrames = 2;
activationOpts.postMedianFactor = 1.35;
activationOpts.postMaxFactor = 1.5;
activationOpts.enableBurstFallback = true;
activationOpts.burstJumpFactor = 2.2;
activationOpts.burstPostFactor = 1.4;

% Density plot bin size (physical units)
binSize_phys = 0.02;    % <---edit if needed based on mm/px and field size

% Centralized plot/export configuration
plotOpts = struct();
plotOpts.enableNormalTheme = true; % <---edit
plotOpts.enablePosterTheme = false; % <---edit
plotOpts.savePNG = true; % <---edit
plotOpts.saveSVG = false; % <---edit
plotOpts.makeTrackDiagnostics = true; % <---edit
plotOpts.saveDiagnosticGifs = true; % <---edit
plotOpts.makeVideoOverlayGifs = false; % <---edit: overlay tracks on source AVI and export GIF
plotOpts.saveVideoOverlayGifs = true; % <---edit
plotOpts.diagnosticGifTrailLength = 10;
plotOpts.diagnosticGifDelayTime = 0.08;
plotOpts.videoOverlayTrailLength = 10; % <---edit
plotOpts.videoOverlayDelayTime = 0.08; % <---edit
plotOpts.videoOverlayFadeHalfLifeFrames = 30; % <---edit
plotOpts.videoOverlayUseContiguousRange = true; % <---edit
plotOpts.videoOverlayMaxFrames = 600; % <---edit (Inf for all)
plotOpts.videoOverlayMarkerSize = 26; % <---edit
plotOpts.upstreamSizeXLim_um = []; % <---edit
plotOpts.inceptionImageSize_px = [1280 320]; % [width height]
plotOpts.inceptionXLim_mm = [0 5]; % <---edit
plotOpts.inceptionYLim_mm = [0 1.2]; % <---edit
plotOpts.themes = enabled_plot_themes(plotOpts);

%% ---------------- DEFINE CASES (ONE OR MULTIPLE RE) ----------------
% Required fields per case:
%   name, Re, kDh, xmlFile, pixelSize, dt
% If you include one Re only, everything still runs.
% Edit case definitions below for your runs. <---edit

cases = struct([]);

cases(1).name      = "5um";
cases(1).Re        = 95000;
cases(1).kDh       = 5; % <-- set your k/D_h
cases(1).xmlFile   = "E:\March Re 90,000 inception data\Processed images\Smooth variation 2\smoothvar2_48lit.xml"; % <-- set XML path
cases(1).pixelSize = 0.00375009375;  % mm/px
cases(1).dt        = 1/102247;

cases(2).name      = "12um";
cases(2).Re        = 95000;
cases(2).kDh       = 12;
cases(2).xmlFile   = "E:\March Re 90,000 inception data\Processed images\P10S100\P10S100_48lit.xml";
cases(2).pixelSize = 0.00375009375;
cases(2).dt        = 1/102247;

cases(3).name      = "20um";
cases(3).Re        = 95000;
cases(3).kDh       = 20;
cases(3).xmlFile   = "E:\March Re 90,000 inception data\Processed images\P10S70\P10S70_48lit.xml";
cases(3).pixelSize = 0.00375009375;
cases(3).dt        = 1/102247;

cases(4).name      = "30um";
cases(4).Re        = 95000;
cases(4).kDh       = 30;
cases(4).xmlFile   = "E:\March Re 90,000 inception data\Processed images\P10S50\P10S50_48lit.xml";
cases(4).pixelSize = 0.00375009375;
cases(4).dt        = 1/102247;

cases(5).name      = "53um";
cases(5).Re        = 95000;
cases(5).kDh       = 53;
cases(5).xmlFile   = "E:\March Re 90,000 inception data\Processed images\P10S30\P10S30_48lit.xml";
cases(5).pixelSize = 0.00375009375;
cases(5).dt        = 1/102247;

cases(6).name      = "80um";
cases(6).Re        = 95000;
cases(6).kDh       = 80;
cases(6).xmlFile   = "E:\March Re 90,000 inception data\Processed images\P10S20\P10S20_48lit.xml";
cases(6).pixelSize = 0.00375009375;
cases(6).dt        = 1/102247;

for ci = 1:numel(cases)
    cases(ci).diagnosticTrackIds = []; % empty = all parsed tracks in this case; otherwise list TRACK_IDs
    if ~isfield(cases, 'videoFile')
        cases(ci).videoFile = ""; % <---edit: optional AVI path for overlay GIF export
    elseif strlength(string(cases(ci).videoFile)) < 1
        cases(ci).videoFile = ""; % <---edit: optional AVI path for overlay GIF export
    end
end

% Example:
% cases(2).diagnosticTrackIds = [123 456];
% cases(1).videoFile = "E:\path\to\case1.avi"; % <---edit

% Example second Reynolds set (uncomment/edit as needed)
% cases(7).name      = "5um";
% cases(7).Re        = 124000;
% cases(7).kDh       = 0.0005;
% cases(7).xmlFile   = "<path_to_Re124k_case_5um.xml>";
% cases(7).pixelSize = 0.00375009375;
% cases(7).dt        = 1/102247;

[cases, selectedCaseIdx, totalCaseCount] = select_cases(cases, caseSelection);
selectedLabels = strings(numel(cases), 1);
for si = 1:numel(cases)
    selectedLabels(si) = sprintf('%d:%s', selectedCaseIdx(si), char(cases(si).name));
end
fprintf('Selected %d/%d case(s): %s\n', numel(cases), totalCaseCount, strjoin(cellstr(selectedLabels), ', '));

%% ---------------- RUN ALL CASES ----------------
summaryRows = table();

allLoc = struct();
allLoc.caseName = strings(0,1);
allLoc.Re       = nan(0,1);
allLoc.kDh      = nan(0,1);
allLoc.pixelSize = nan(0,1);
allLoc.inception2x_xy = cell(0,1);  % [x y] activation points on left-moving tracks

trackFigOutDir = fullfile(figDir, "TrackDiagnostics");
if ~isfolder(trackFigOutDir), mkdir(trackFigOutDir); end

gifOutDir = fullfile(resultsDir, "diagnostic gifs");
if plotOpts.saveDiagnosticGifs && ~isfolder(gifOutDir), mkdir(gifOutDir); end

videoGifOutDir = fullfile(resultsDir, "video overlay gifs");
if plotOpts.makeVideoOverlayGifs && plotOpts.saveVideoOverlayGifs && ~isfolder(videoGifOutDir), mkdir(videoGifOutDir); end

allSize = struct();
allSize.caseName = strings(0,1);
allSize.Re       = nan(0,1);
allSize.kDh      = nan(0,1);
allSize.size_eqd = cell(0,1);

gateSummaryRows = table();
cachePolicyTag = build_cache_policy_tag(parserOpts, qcOpts, flowOpts, activationOpts);
runTimer = tic;

for i = 1:numel(cases)
    caseTimer = tic;
    fprintf("\n=== Case %d/%d: %s (Re=%g, k/D_h=%.6g) ===\n", ...
        i, numel(cases), cases(i).name, cases(i).Re, cases(i).kDh);

    caseKey = build_case_cache_key(cases(i), maxTracksToParse, cachePolicyTag);
    idxCache = find(cacheDB.key == caseKey, 1, 'first');
    useCacheEntry = false;

    if useMatCache && ~forceReparse && ~isempty(idxCache)
        cachedOut = cacheDB.out{idxCache};
        if is_cache_entry_compatible(cachedOut, parserOpts)
            out = cachedOut;
            useCacheEntry = true;
            fprintf("Using cache entry %d/%d\n", idxCache, numel(cacheDB.key));
        else
            fprintf("Cache entry %d/%d incompatible with parser policy/version. Reparsing XML.\n", ...
                idxCache, numel(cacheDB.key));
        end
    end

    if ~useCacheEntry
        % Parse XML (NO plotting inside)
        out = analyze_trackmate_xml_arc( ...
            cases(i).xmlFile, ...
            pixelSize = cases(i).pixelSize, ...
            dt        = cases(i).dt, ...
            maxTracks = maxTracksToParse, ...
            parseROI  = false, ...
            verbose   = ~isArc, ...
            makePlots = false, ...
            parseTrackedSpotsOnly = parserOpts.parseTrackedSpotsOnly, ...
            parseFilteredTracksOnly = parserOpts.parseFilteredTracksOnly);

        if useMatCache
            if isempty(idxCache)
                cacheDB.key(end+1,1) = caseKey;
                cacheDB.out{end+1,1} = out;
            else
                cacheDB.out{idxCache,1} = out;
            end
            cacheUpdated = true;
        end
    end

    % Compute injection/activation metrics + locations
    metrics = trackmate_case_metrics(out, qcOpts, flowOpts, activationOpts);
    g = metrics.gateStats;
    fprintf(['Strict gate summary: total=%d, injected=%d, activated=%d, ', ...
        'reject(short=%d nonfinite=%d nonmonoT=%d topo=%d flow=%d origin=%d noAct=%d wall=%d)\n'], ...
        g.nTracksTotal, g.nInjected, g.nActivated, ...
        g.nRejectedTooShort, g.nRejectedNonFinite, g.nRejectedNonMonotonicTime, ...
        g.nRejectedTopology, g.nRejectedFlow, g.nRejectedOrigin, ...
        g.nRejectedNoActivation, g.nRejectedWallBand);

    nValidTracks = metrics.nBasicValidTracks;
    nStrictRecirculationTracks = metrics.nStrictPrimaryTracks;
    nStrictActivatedTracks = metrics.nStrictActivatedTracks;
    strictRecirculationFrac_total = nStrictRecirculationTracks / max(metrics.nTracksTotal, 1);
    strictActivationFrac_valid = nStrictActivatedTracks / max(nValidTracks, 1);
    tau_mean_val = metrics.tau_mean;
    tau_std_val = metrics.tau_std;
    nTauVals = numel(metrics.tau_values);
    A_over_I_err_low = metrics.A_over_I - metrics.A_over_I_ci_low;
    A_over_I_err_high = metrics.A_over_I_ci_high - metrics.A_over_I;
    if nTauVals > 1
        tau_sem = tau_std_val / sqrt(nTauVals);
    else
        tau_sem = NaN;
    end

    % Accumulate inception (2x growth) locations for plotting
    allLoc.caseName(end+1,1) = string(cases(i).name);
    allLoc.Re(end+1,1)       = cases(i).Re;
    allLoc.kDh(end+1,1)      = cases(i).kDh;
    allLoc.pixelSize(end+1,1) = cases(i).pixelSize;
    allLoc.inception2x_xy{end+1,1} = choose_inception_activation_xy(metrics);

    % Accumulate upstream-size samples for distribution plot
    allSize.caseName(end+1,1) = string(cases(i).name);
    allSize.Re(end+1,1)       = cases(i).Re;
    allSize.kDh(end+1,1)      = cases(i).kDh;
    allSize.size_eqd{end+1,1} = metrics.upstreamSize_eqd;

    gateRow = table( ...
        string(cases(i).name), cases(i).Re, cases(i).kDh, ...
        g.nTracksTotal, g.nInjected, g.nActivated, ...
        g.nRejectedTooShort, g.nRejectedNonFinite, g.nRejectedNonMonotonicTime, ...
        g.nRejectedTopology, g.nRejectedFlow, g.nRejectedOrigin, ...
        g.nRejectedNoActivation, g.nRejectedWallBand, ...
        g.originThreshold, g.xStartMin, g.xStartMax, ...
        'VariableNames', {'Case','Re','kDh','nTracksTotal','nInjected','nActivated', ...
        'nRejectedTooShort','nRejectedNonFinite','nRejectedNonMonotonicTime', ...
        'nRejectedTopology','nRejectedFlow','nRejectedOrigin', ...
        'nRejectedNoActivation','nRejectedWallBand', ...
        'originThreshold','xStartMin','xStartMax'});
    gateSummaryRows = [gateSummaryRows; gateRow]; %#ok<AGROW>

    write_framewise_case_csv(resultsDir, cases(i), metrics);
    run_diagnostic_parity_checks(cases(i), metrics);

    if plotOpts.makeTrackDiagnostics
        plot_verification_tracks_for_case(cases(i), metrics, trackFigOutDir, plotOpts);
    end

    if plotOpts.saveDiagnosticGifs
        save_diagnostic_track_gifs(cases(i), metrics, gifOutDir, plotOpts);
    end

    if plotOpts.makeVideoOverlayGifs
        save_video_overlay_gif_from_avi(cases(i), metrics, videoGifOutDir, plotOpts);
    end

    elapsed_case_sec = toc(caseTimer);
    fprintf('Case elapsed: %.2f s (%s)\n', elapsed_case_sec, format_elapsed_hms(elapsed_case_sec));

    % Store per-case summary
    row = table( ...
        string(cases(i).name), cases(i).Re, cases(i).kDh, ...
        metrics.nTracksTotal, nValidTracks, ...
        nStrictRecirculationTracks, nStrictActivatedTracks, ...
        metrics.strictTrackFrameExposure, metrics.strictActivationEventsTotal, ...
        strictRecirculationFrac_total, strictActivationFrac_valid, ...
        metrics.A_over_I, metrics.A_over_I_ci_low, metrics.A_over_I_ci_high, A_over_I_err_low, A_over_I_err_high, ...
        tau_mean_val, tau_std_val, tau_sem, nTauVals, elapsed_case_sec, ...
        'VariableNames', {'Case','Re','kDh','nTracksTotal','nValidTracks', ...
        'nStrictRecirculationTracks','nStrictActivatedTracks', ...
        'strictTrackFrameExposure','strictActivationEventsTotal', ...
        'strictRecirculationFrac_total','strictActivationFrac_valid','A_over_I','A_over_I_ci_low','A_over_I_ci_high','A_over_I_err_low','A_over_I_err_high', ...
        'tau_mean','tau_std','tau_sem','nTau','elapsed_case_sec'});

    summaryRows = [summaryRows; row]; %#ok<AGROW>

    % Close any accidental figures
    close all force;
end

% Sort by Re then k/D_h
summaryRows = sortrows(summaryRows, {'Re','kDh'});

% Save summary table
summaryCsv = fullfile(resultsDir, "Summary_tracks.csv");
write_table_csv_compat(summaryRows, summaryCsv);
fprintf("\nSaved: %s\n", summaryCsv);

if ~isempty(gateSummaryRows)
    gateSummaryRows = sortrows(gateSummaryRows, {'Re','kDh'});
    gateSummaryCsv = fullfile(resultsDir, "track_gate_summary.csv");
    write_table_csv_compat(gateSummaryRows, gateSummaryCsv);
    fprintf("Saved: %s\n", gateSummaryCsv);
end

%% ---------------- PLOT 1: A/I vs k/D_h (per Re) ----------------
fitTxtFile = fullfile(resultsDir, "fit_AI_vs_kDh_by_Re.txt");
plot_ai_vs_kdh_re(summaryRows, figDir, fitTxtFile, plotOpts);

%% ---------------- PLOT 2: Inception (2x growth) locations ----------------
locFigOutDir = fullfile(figDir, "InceptionLocations");
if ~isfolder(locFigOutDir), mkdir(locFigOutDir); end
plot_inception_locations_by_re(allLoc, locFigOutDir, plotOpts);

%% ---------------- PLOT 3: mean residence time Tau vs k/D_h ----------------
plot_tau_vs_kdh_re(summaryRows, figDir, plotOpts);

%% ---------------- PLOT 4: upstream moving microbubble size distributions ----------------
distFigOutDir = fullfile(figDir, "UpstreamSizeDistributions");
if ~isfolder(distFigOutDir), mkdir(distFigOutDir); end
plot_upstream_size_distribution_by_re(allSize, distFigOutDir, binSize_phys, plotOpts);

if useMatCache && cacheUpdated
    save(cacheFile, 'cacheDB', '-v7.3');
    fprintf("Saved cache DB: %s (%d entries)\n", cacheFile, numel(cacheDB.key));
end

totalElapsedSec = toc(runTimer);
fprintf("Total elapsed time (selected cases): %.2f s (%s)\n", totalElapsedSec, format_elapsed_hms(totalElapsedSec));
fprintf("\nAll done. Results in: %s\n", resultsDir);

function [xMean, xStd, yMean, yStd, nPts] = xy_stats(xy)
nPts = 0;
xMean = NaN;
xStd = NaN;
yMean = NaN;
yStd = NaN;

if isempty(xy) || size(xy,2) < 2
    return;
end

xy = xy(isfinite(xy(:,1)) & isfinite(xy(:,2)), :);
nPts = size(xy,1);
if nPts == 0
    return;
end

xVals = xy(:,1);
yVals = xy(:,2);
xMean = mean(xVals, 'omitnan');
xStd = std(xVals, 0, 'omitnan');
yMean = mean(yVals, 'omitnan');
yStd = std(yVals, 0, 'omitnan');
end

function [m, s, md, p90, n] = vector_stats(v)
n = 0;
m = NaN;
s = NaN;
md = NaN;
p90 = NaN;

if isempty(v)
    return;
end

v = v(:);
v = v(isfinite(v));
n = numel(v);
if n == 0
    return;
end

m = mean(v, 'omitnan');
s = std(v, 0, 'omitnan');
md = median(v, 'omitnan');
p90 = percentile_legacy(v, 90);
end

function q = percentile_legacy(v, p)
q = NaN;
if isempty(v) || ~isfinite(p)
    return;
end

v = sort(v(:));
n = numel(v);
if n == 0
    return;
end
if n == 1
    q = v(1);
    return;
end

p = min(100, max(0, p));
idx = 1 + (n - 1) * (p / 100);
i0 = floor(idx);
i1 = ceil(idx);
if i0 == i1
    q = v(i0);
else
    frac = idx - i0;
    q = v(i0) + frac * (v(i1) - v(i0));
end
end

function write_framewise_case_csv(resultsDir, caseDef, metrics)
strictAxis = get_metric_vector(metrics, 'strict_frame_axis', metrics.frame_axis);
strictVisible = get_metric_vector(metrics, 'strict_frame_nVisible', metrics.frame_nLeftMovingVisible);
strictAct = get_metric_vector(metrics, 'strict_frame_nActivationEvents', metrics.frame_nActivationEvents);
strictCumExposure = get_metric_vector(metrics, 'strict_frame_cumExposure', metrics.frame_cumExposure);
strictCumAct = get_metric_vector(metrics, 'strict_frame_cumActivationEvents', metrics.frame_cumActivationEvents);

legacyAxis = get_metric_vector(metrics, 'frame_axis_netLeftLegacy', nan(0,1));
legacyVisible = get_metric_vector(metrics, 'frame_nLeftMovingVisible_netLeftLegacy', nan(0,1));
legacyAct = get_metric_vector(metrics, 'frame_nActivationEvents_netLeftLegacy', nan(0,1));
legacyCumExposure = get_metric_vector(metrics, 'frame_cumExposure_netLeftLegacy', nan(0,1));
legacyCumAct = get_metric_vector(metrics, 'frame_cumActivationEvents_netLeftLegacy', nan(0,1));

frameAxis = unique([strictAxis(:); legacyAxis(:)]);
frameAxis = frameAxis(isfinite(frameAxis));
if isempty(frameAxis)
    frameAxis = nan(0,1);
end

n = numel(frameAxis);
nStrictVisible = zeros(n,1);
nStrictAct = zeros(n,1);
cumStrictExposure = zeros(n,1);
cumStrictAct = zeros(n,1);
nLegacyVisible = zeros(n,1);
nLegacyAct = zeros(n,1);
cumLegacyExposure = zeros(n,1);
cumLegacyAct = zeros(n,1);

for ii = 1:n
    f = frameAxis(ii);
    nStrictVisible(ii) = lookup_exact(strictAxis, strictVisible, f);
    nStrictAct(ii) = lookup_exact(strictAxis, strictAct, f);
    cumStrictExposure(ii) = lookup_cumulative(strictAxis, strictCumExposure, f);
    cumStrictAct(ii) = lookup_cumulative(strictAxis, strictCumAct, f);

    nLegacyVisible(ii) = lookup_exact(legacyAxis, legacyVisible, f);
    nLegacyAct(ii) = lookup_exact(legacyAxis, legacyAct, f);
    cumLegacyExposure(ii) = lookup_cumulative(legacyAxis, legacyCumExposure, f);
    cumLegacyAct(ii) = lookup_cumulative(legacyAxis, legacyCumAct, f);
end

frameTbl = table( ...
    frameAxis, ...
    nStrictVisible, ...
    nStrictAct, ...
    cumStrictExposure, ...
    cumStrictAct, ...
    nLegacyVisible, ...
    nLegacyAct, ...
    cumLegacyExposure, ...
    cumLegacyAct, ...
    'VariableNames', {'frame', ...
    'nStrictVisible','nStrictActivationEvents','cumStrictExposure','cumStrictActivationEvents', ...
    'nLeftMovingVisible_netLeftLegacy','nActivationEvents_netLeftLegacy', ...
    'cumExposure_netLeftLegacy','cumActivationEvents_netLeftLegacy'});

caseToken = sanitize_case_token(caseDef.name);
outCsv = fullfile(resultsDir, sprintf('framewise_counts_%s_Re_%g_kDh_%g.csv', caseToken, caseDef.Re, caseDef.kDh));
write_table_csv_compat(frameTbl, outCsv);
fprintf("Saved: %s\n", outCsv);
end

function v = get_metric_vector(metrics, fieldName, fallback)
v = fallback;
if isfield(metrics, fieldName)
    raw = metrics.(fieldName);
    if ~isempty(raw)
        v = raw(:);
    end
end
if isempty(v)
    v = nan(0,1);
end
end

function v = lookup_exact(axisVals, dataVals, frameVal)
v = 0;
if isempty(axisVals) || isempty(dataVals)
    return;
end
axisVals = axisVals(:);
dataVals = dataVals(:);
if numel(axisVals) ~= numel(dataVals)
    return;
end
idx = find(axisVals == frameVal, 1, 'first');
if isempty(idx) || ~isfinite(dataVals(idx))
    return;
end
v = dataVals(idx);
end

function v = lookup_cumulative(axisVals, dataVals, frameVal)
v = 0;
if isempty(axisVals) || isempty(dataVals)
    return;
end
axisVals = axisVals(:);
dataVals = dataVals(:);
if numel(axisVals) ~= numel(dataVals)
    return;
end
idx = find(axisVals <= frameVal, 1, 'last');
if isempty(idx) || ~isfinite(dataVals(idx))
    return;
end
v = dataVals(idx);
end

function token = sanitize_case_token(caseName)
token = char(string(caseName));
token = regexprep(token, '[^A-Za-z0-9]+', '_');
token = regexprep(token, '^_+', '');
token = regexprep(token, '_+$', '');
if isempty(token)
    token = 'case';
end
end

function xy = choose_inception_activation_xy(metrics)
xy = zeros(0,2);
if isfield(metrics, 'activationEvent_xy_netLeftLegacy') && ~isempty(metrics.activationEvent_xy_netLeftLegacy)
    xy = metrics.activationEvent_xy_netLeftLegacy;
elseif isfield(metrics, 'activationEvent_xy') && ~isempty(metrics.activationEvent_xy)
    xy = metrics.activationEvent_xy;
elseif isfield(metrics, 'inception2x_xy') && ~isempty(metrics.inception2x_xy)
    xy = metrics.inception2x_xy;
end
end

function s = format_elapsed_hms(secondsVal)
if ~isfinite(secondsVal) || secondsVal < 0
    s = 'n/a';
    return;
end
totalSeconds = round(secondsVal);
h = floor(totalSeconds / 3600);
m = floor(mod(totalSeconds, 3600) / 60);
sSec = mod(totalSeconds, 60);
s = sprintf('%02d:%02d:%02d', h, m, sSec);
end

function v = get_metric_field(metrics, fieldName, defaultValue)
v = defaultValue;
if isfield(metrics, fieldName)
    raw = metrics.(fieldName);
    if ~isempty(raw)
        v = raw;
    end
end
end

function run_diagnostic_parity_checks(caseDef, metrics)
if ~isstruct(metrics) || ~isfield(metrics, 'trackCatalog') || isempty(metrics.trackCatalog)
    return;
end

trackRequests = [];
if isfield(caseDef, 'diagnosticTrackIds') && ~isempty(caseDef.diagnosticTrackIds)
    trackRequests = caseDef.diagnosticTrackIds(:).';
end
[selectedIdx, selectedTrackIds] = resolve_diagnostic_track_indices(metrics.trackCatalog, trackRequests, caseDef, 'Parity check');
if isempty(selectedIdx)
    return;
end

strictMask = false(numel(selectedIdx), 1);
for i = 1:numel(selectedIdx)
    tr = metrics.trackCatalog(selectedIdx(i));
    strictMask(i) = is_catalog_true(tr, 'isStrictPrimary');
end
strictIds = unique(selectedTrackIds(strictMask));

if isempty(trackRequests) && isfield(metrics, 'nStrictPrimaryTracks')
    if numel(strictIds) ~= metrics.nStrictPrimaryTracks
        warning(['Parity check failed for %s: strict diagnostic ID count (%d) ', ...
            '!= metrics.nStrictPrimaryTracks (%d).'], ...
            char(caseDef.name), numel(strictIds), metrics.nStrictPrimaryTracks);
    end
end

eventIds = nan(0,1);
if isfield(metrics, 'strictActivationEvent_trackId')
    eventIds = metrics.strictActivationEvent_trackId(:);
elseif isfield(metrics, 'activationEvent_trackId')
    eventIds = metrics.activationEvent_trackId(:);
end
eventIds = eventIds(isfinite(eventIds));

if isempty(trackRequests) && isfield(metrics, 'strictActivationEventsTotal')
    if numel(eventIds) ~= metrics.strictActivationEventsTotal
        warning(['Parity check failed for %s: strict activation events (%d) ', ...
            '!= metrics.strictActivationEventsTotal (%d).'], ...
            char(caseDef.name), numel(eventIds), metrics.strictActivationEventsTotal);
    end
end
end

function tf = is_catalog_true(tr, fieldName)
tf = false;
if ~isstruct(tr) || ~isfield(tr, fieldName)
    return;
end
v = tr.(fieldName);
if islogical(v)
    tf = any(v(:));
elseif isnumeric(v)
    tf = any(v(:) ~= 0);
end
end

function caseKey = build_case_cache_key(caseDef, maxTracks, policyTag)
if isfile(caseDef.xmlFile)
    d = dir(caseDef.xmlFile);
    fileSig = sprintf("%d_%0.0f", d.bytes, d.datenum * 1e6);
else
    fileSig = "missing";
end

maxTracksStr = "all";
if isfinite(maxTracks)
    maxTracksStr = string(max(0, floor(maxTracks)));
end

if nargin < 3 || isempty(policyTag)
    policyTag = "default";
end

caseKey = sprintf("xml=%s|sig=%s|px=%.12g|dt=%.12g|max=%s|policy=%s", ...
    string(caseDef.xmlFile), fileSig, caseDef.pixelSize, caseDef.dt, maxTracksStr, string(policyTag));
end

function tf = is_cache_entry_compatible(out, parserOpts)
tf = false;
if ~isstruct(out) || ~isfield(out, 'meta') || ~isstruct(out.meta)
    return;
end

meta = out.meta;
if ~isfield(meta, 'parserVersion') || ~isfinite(meta.parserVersion)
    return;
end
if meta.parserVersion < parserOpts.parserVersion
    return;
end

if parserOpts.parseTrackedSpotsOnly
    if ~isfield(meta, 'parseTrackedSpotsOnly') || ~is_true_flag(meta.parseTrackedSpotsOnly)
        return;
    end
end

if parserOpts.parseFilteredTracksOnly
    if ~isfield(meta, 'parseFilteredTracksOnly') || ~is_true_flag(meta.parseFilteredTracksOnly)
        return;
    end
end

tf = true;
end

function tf = is_true_flag(v)
tf = false;
if islogical(v)
    tf = any(v(:));
elseif isnumeric(v)
    tf = any(v(:) ~= 0);
elseif isstring(v) || ischar(v)
    s = string(v);
    tf = any(strcmpi(s, "true") | strcmpi(s, "1"));
end
end

function policyTag = build_cache_policy_tag(parserOpts, qcOpts, flowOpts, activationOpts)
policyTag = sprintf([ ...
    'pv=%d|pts=%d|pft=%d|bulk=%s|dx=%.6g|neg=%.4g|pos=%.4g|origin=%d|originFrac=%.4g|', ...
    'spots=%d|gaps=%.4g|rsmc=%d|wall=%d|jump=%.4g|pre=%d|minpre=%d|post=%d|', ...
    'postMed=%.4g|postMax=%.4g|burst=%d|burstJump=%.4g|burstPost=%.4g'], ...
    parserOpts.parserVersion, double(logical(parserOpts.parseTrackedSpotsOnly)), double(logical(parserOpts.parseFilteredTracksOnly)), ...
    char(string(flowOpts.bulkDirection)), flowOpts.minNetDxCounterflow_mm, ...
    flowOpts.minNegativeStepFraction, flowOpts.maxPositiveStepFraction, ...
    double(logical(flowOpts.requireRightOrigin)), flowOpts.rightOriginFrac, ...
    qcOpts.minTrackSpots, qcOpts.maxTrackGaps, double(logical(qcOpts.rejectSplitMergeComplex)), ...
    double(logical(qcOpts.wallBandEnabled)), activationOpts.areaJumpFactor, ...
    activationOpts.preWindowFrames, activationOpts.minPrePoints, activationOpts.requiredPostFrames, ...
    activationOpts.postMedianFactor, activationOpts.postMaxFactor, ...
    double(logical(activationOpts.enableBurstFallback)), activationOpts.burstJumpFactor, ...
    activationOpts.burstPostFactor);
end

function [casesOut, selectedIdx, totalCaseCount] = select_cases(casesIn, selection)
totalCaseCount = numel(casesIn);
allIdx = (1:totalCaseCount).';

if nargin < 2 || isempty(selection)
    selection = "all";
end

if islogical(selection)
    mask = selection(:);
    if numel(mask) ~= totalCaseCount
        error('Logical caseSelection must have %d elements.', totalCaseCount);
    end
    selectedIdx = find(mask);
elseif isnumeric(selection)
    vals = selection(:);
    vals = vals(isfinite(vals));
    if isempty(vals)
        selectedIdx = allIdx;
    else
        vals = round(vals);
        if any(vals < 1) || any(vals > totalCaseCount)
            error('Numeric caseSelection indices must be in [1, %d].', totalCaseCount);
        end
        selectedIdx = unique(vals, 'stable');
    end
else
    if iscell(selection)
        selection = string(selection);
    end
    tokens = string(selection);
    tokens = tokens(:);
    tokens = strtrim(tokens);
    tokens = tokens(strlength(tokens) > 0);
    if isempty(tokens)
        selectedIdx = allIdx;
    elseif any(strcmpi(tokens, "all"))
        selectedIdx = allIdx;
    else
        caseNames = strings(totalCaseCount, 1);
        for i = 1:totalCaseCount
            caseNames(i) = string(casesIn(i).name);
        end
        selectedIdx = zeros(0,1);
        for ti = 1:numel(tokens)
            idx = find(strcmpi(caseNames, tokens(ti)));
            if isempty(idx)
                error("Unknown caseSelection '%s'. Valid names: %s", ...
                    char(tokens(ti)), strjoin(cellstr(caseNames), ", "));
            end
            selectedIdx = [selectedIdx; idx(:)]; %#ok<AGROW>
        end
        selectedIdx = unique(selectedIdx, 'stable');
    end
end

if isempty(selectedIdx)
    error('caseSelection resolved to 0 cases. Please select at least one case.');
end

casesOut = casesIn(selectedIdx);
end

function themes = enabled_plot_themes(plotOpts)
themes = strings(0,1);
if plotOpts.enableNormalTheme
    themes(end+1,1) = "normal"; %#ok<AGROW>
end
if plotOpts.enablePosterTheme
    themes(end+1,1) = "poster"; %#ok<AGROW>
end
if isempty(themes)
    themes = "normal";
end
end

function bootstrap_user_matlab_path()
homeDir = getenv('HOME');
if isempty(homeDir)
    homeDir = getenv('USERPROFILE');
end
if isempty(homeDir)
    warning('Could not resolve HOME/USERPROFILE. Skipping user MATLAB path bootstrap.');
    return;
end

userPathFile = fullfile(homeDir, 'matlab', 'pathdef.m');
userCodeRoot = fullfile(homeDir, 'matlab');

if usejava('desktop')
    if isfile(userPathFile)
        try
            run(userPathFile);
            fprintf('Interactive mode: loaded custom MATLAB path:\n  %s\n', userPathFile);
        catch ME
            warning('Interactive mode: failed to run pathdef.m (%s).', ME.message);
        end
    else
        warning('Interactive mode: pathdef.m not found. Using default MATLAB path.');
    end
    return;
end

loadedPathdef = false;
if isfile(userPathFile)
    try
        run(userPathFile);
        loadedPathdef = true;
        fprintf('Headless mode: loaded custom MATLAB path:\n  %s\n', userPathFile);
    catch ME
        warning('Headless mode: failed to run pathdef.m (%s).', ME.message);
    end
end

if isfolder(userCodeRoot)
    addpath(genpath(userCodeRoot));
    fprintf('Headless mode: added user MATLAB root:\n  %s\n', userCodeRoot);
elseif ~loadedPathdef
    warning('Headless mode: user MATLAB root not found: %s', userCodeRoot);
end
end
