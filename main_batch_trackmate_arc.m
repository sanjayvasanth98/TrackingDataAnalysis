%% main_batch_trackmate_arc.m
% Batch TrackMate XML analysis (ARC + local)
% - Computes Injection (upstream right->left) and Activation/Growth (AREA jump >= factor)
% - Supports one or multiple Reynolds-number sets for the same roughness ladder
% - Plots A/I vs k/D_h, inception (2x growth) locations, Tau vs k/D_h,
%   and upstream moving microbubble size distributions
%
% Output:
%   resultsDir/
%     summary_AI_tau_vs_kDh_Re.csv
%     fit_AI_vs_kDh_by_Re.txt
%     Figures_PNG_SVG/(normal|poster)/...

clear; clc;

%% ---------------- USER SETTINGS ----------------
% Toggle run mode by commenting one line
runMode = "local";
% runMode = "arc";

isArc = strcmpi(runMode, "arc");
visMode = ["on", "off"];
set(0,'DefaultFigureVisible', visMode(1 + isArc));

% Local quick-test limit (tracks/events parsed per XML)
% - Use 100 (or any number) for fast local testing
% - Use Inf to parse all tracks
localEventLimit = 200;
maxTracksChoices = [localEventLimit, Inf];
maxTracksToParse = maxTracksChoices(1 + isArc);

% Cache parsed outputs (.mat) to avoid re-parsing XML on reruns
useMatCache = true;
forceReparse = false;

% Where to save all results
resultsDir = "E:\March Re 90,000 inception data\Processed images\results\results 2";
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
flowOpts.bulkDirection = "left_to_right";
flowOpts.minNetDxCounterflow_mm = 0.08;
flowOpts.minNegativeStepFraction = 0.65;
flowOpts.maxPositiveStepFraction = 0.30;
flowOpts.requireRightOrigin = true;
flowOpts.rightOriginFrac = 0.60; % rightmost 40% source band

qcOpts = struct();
qcOpts.minTrackSpots = 5;
qcOpts.maxTrackGaps = 1;
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

% Legacy scalar knobs kept for diagnostic GIF helper signatures.
minNetDx_phys = flowOpts.minNetDxCounterflow_mm;
minTrackSpots = qcOpts.minTrackSpots;
areaJumpFactor = activationOpts.areaJumpFactor;

% Density plot bin size (physical units)
binSize_phys = 0.02;    % adjust based on your mm/px and field size (e.g., 0.02 mm)

% Centralized plot/export configuration
plotOpts = struct();
plotOpts.enableNormalTheme = true;
plotOpts.enablePosterTheme = false;
plotOpts.savePNG = true;
plotOpts.saveSVG = false;
plotOpts.makeTrackDiagnostics = true;
plotOpts.saveDiagnosticGifs = true;
plotOpts.diagnosticGifTrailLength = 10;
plotOpts.diagnosticGifDelayTime = 0.08;
plotOpts.upstreamSizeXLim_um = [];
plotOpts.inceptionImageSize_px = [1280 320]; % [width height]
plotOpts.inceptionXLim_mm = [0 5];
plotOpts.themes = enabled_plot_themes(plotOpts);

%% ---------------- DEFINE CASES (ONE OR MULTIPLE RE) ----------------
% Required fields per case:
%   name, Re, kDh, xmlFile, pixelSize, dt
% If you include one Re only, everything still runs.

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
end

% Example:
% cases(2).diagnosticTrackIds = [123 456];

% Example second Reynolds set (uncomment/edit as needed)
% cases(7).name      = "5um";
% cases(7).Re        = 124000;
% cases(7).kDh       = 0.0005;
% cases(7).xmlFile   = "<path_to_Re124k_case_5um.xml>";
% cases(7).pixelSize = 0.00375009375;
% cases(7).dt        = 1/102247;

%% ---------------- RUN ALL CASES ----------------
summaryRows = table();

allLoc = struct();
allLoc.caseName = strings(0,1);
allLoc.Re       = nan(0,1);
allLoc.kDh      = nan(0,1);
allLoc.pixelSize = nan(0,1);
allLoc.inception2x_xy = cell(0,1);  % [x y] where first growth jump >= factor

trackFigOutDir = fullfile(figDir, "TrackDiagnostics");
if ~isfolder(trackFigOutDir), mkdir(trackFigOutDir); end

gifOutDir = fullfile(resultsDir, "diagnostic gifs");
if plotOpts.saveDiagnosticGifs && ~isfolder(gifOutDir), mkdir(gifOutDir); end

allSize = struct();
allSize.caseName = strings(0,1);
allSize.Re       = nan(0,1);
allSize.kDh      = nan(0,1);
allSize.size_eqd = cell(0,1);

gateSummaryRows = table();
cachePolicyTag = build_cache_policy_tag(parserOpts, qcOpts, flowOpts, activationOpts);

for i = 1:numel(cases)
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
    fprintf(['Gate summary: total=%d, injected=%d, activated=%d, ', ...
        'reject(short=%d nonfinite=%d nonmonoT=%d topo=%d flow=%d origin=%d noAct=%d wall=%d)\n'], ...
        g.nTracksTotal, g.nInjected, g.nActivated, ...
        g.nRejectedTooShort, g.nRejectedNonFinite, g.nRejectedNonMonotonicTime, ...
        g.nRejectedTopology, g.nRejectedFlow, g.nRejectedOrigin, ...
        g.nRejectedNoActivation, g.nRejectedWallBand);

    % Store per-case summary
    row = table( ...
        string(cases(i).name), cases(i).Re, cases(i).kDh, ...
        metrics.nTracksTotal, metrics.nInjected, metrics.nActivated, ...
        metrics.A_over_I, metrics.A_over_I_ci_low, metrics.A_over_I_ci_high, ...
        metrics.tau_mean, metrics.tau_std, numel(metrics.tau_values), ...
        'VariableNames', {'Case','Re','kDh','nTracksTotal','nInjected','nActivated', ...
        'A_over_I','A_over_I_ci_low','A_over_I_ci_high','tau_mean','tau_std','nTau'});

    summaryRows = [summaryRows; row]; %#ok<AGROW>

    % Accumulate inception (2x growth) locations for plotting
    allLoc.caseName(end+1,1) = string(cases(i).name);
    allLoc.Re(end+1,1)       = cases(i).Re;
    allLoc.kDh(end+1,1)      = cases(i).kDh;
    allLoc.pixelSize(end+1,1) = cases(i).pixelSize;
    allLoc.inception2x_xy{end+1,1} = metrics.inception2x_xy;

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

    if plotOpts.makeTrackDiagnostics
        plot_verification_tracks_for_case(cases(i), metrics, trackFigOutDir, plotOpts);
    end

    if plotOpts.saveDiagnosticGifs
        save_diagnostic_track_gifs(cases(i), out, gifOutDir, plotOpts, minNetDx_phys, minTrackSpots, areaJumpFactor);
    end

    % Close any accidental figures
    close all force;
end

% Sort by Re then k/D_h
summaryRows = sortrows(summaryRows, {'Re','kDh'});

% Save summary table
summaryCsv = fullfile(resultsDir, "summary_AI_tau_vs_kDh_Re.csv");
writetable(summaryRows, summaryCsv);
fprintf("\nSaved: %s\n", summaryCsv);

if ~isempty(gateSummaryRows)
    gateSummaryRows = sortrows(gateSummaryRows, {'Re','kDh'});
    gateSummaryCsv = fullfile(resultsDir, "track_gate_summary.csv");
    writetable(gateSummaryRows, gateSummaryCsv);
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

fprintf("\nAll done. Results in: %s\n", resultsDir);

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
