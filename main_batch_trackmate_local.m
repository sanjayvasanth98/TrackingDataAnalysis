%% main_batch_trackmate_local.m
% Batch TrackMate XML analysis (ARC + local)
% - Computes framewise Injection (left-moving availability) and Activation/Growth events
% - Supports one or multiple Reynolds-number sets for the same roughness ladder
% - Plots A/I vs k/d, inception (2x growth) locations, Tau vs k/d,
%   and upstream moving microbubble size distributions
%
% Output:
%   resultsDir/
%     Summary_tracks.csv
%     fit_AI_vs_kD_by_Re.txt
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

% Local quick-test limit mode: "events" caps left-moving tracks, "frames" caps GIF frames.
% Toggle between the two depending on what you want to inspect.
localLimitMode = "frames"; % <---edit: "events" or "frames"
localEventLimit = 500;     % <---edit: used when localLimitMode = "events"
localFrameLimit = 200;     % <---edit: used when localLimitMode = "frames"

if strcmpi(localLimitMode, "events")
    maxLeftMovingTracks = localEventLimit;
    localGifMaxFrames = Inf;
else
    maxLeftMovingTracks = Inf;
    localGifMaxFrames = localFrameLimit;
end
if isArc
    maxLeftMovingTracks = Inf;
    localGifMaxFrames = Inf;
end

% Case selection:
% - "all"            -> run every case below
% - 1                -> run one case by index
% - [1 3 6]          -> run multiple cases by index
% - "5um"            -> run one case by name
% - ["5um","30um"]   -> run multiple cases by name
caseSelection = "all"; % <---edit

% Cache parsed outputs (.mat) to avoid re-parsing XML on reruns
useMatCache = true; % <---edit
forceReparse = false; % <---edit

% Where to save all results
resultsDir = "E:\March Re 90,000 inception data\Processed images\results\results 29 local"; % <---edit
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
flowOpts.rightOriginFrac = 0.25; % rightmost 75% source band
flowOpts.netLeftMinConsecutiveNegSteps = 2; % <---edit: minimum consecutive counterflow steps for net-left membership
flowOpts.netLeftMinConsecutiveNetDx_mm = 0.00; % <---edit: minimum cumulative counterflow displacement over that run
flowOpts.netLeftBandRescueEnabled = true; % <---edit: special inclusion for short left-moving tracks in a target x-band
flowOpts.netLeftBandRescueX_mm = [0.5 1.5]; % <---edit: x-band for rescue inclusion
flowOpts.netLeftBandRescueMinCounterflowSteps = 1; % <---edit: minimum counterflow steps inside rescue band
flowOpts.netLeftBandRescueMinCounterflowDx_mm = 0.00; % <---edit: minimum cumulative counterflow dx inside rescue band
flowOpts.netLeftBandRescueRequireActivation = true; % <---edit: require track to activate to be rescued
flowOpts.netLeftBandRescueRequireActXInBand = false; % <---edit: if true, activation x must also lie inside rescue band
flowOpts.netLeftBandRescueBypassOriginGate = true; % <---edit: let rescued tracks pass even if they do not satisfy right-origin gate

qcOpts = struct();
qcOpts.minTrackSpots = 4;
qcOpts.maxTrackGaps = 1;
qcOpts.maxLeftMovingTracks = maxLeftMovingTracks;
qcOpts.rejectSplitMergeComplex = false; % <---edit: include breakup/split/merge tracks if they satisfy left-moving/counterflow gates
qcOpts.wallBandEnabled = false;
qcOpts.wallBandYLimits_mm = [];
qcOpts.excludeOriginBoxEnabled = true; % <---edit: exclude tracks that originate inside a forbidden x-y box
qcOpts.excludeOriginBoxX_mm = [0 0.5]; % <---edit: forbidden origin x-range (mm)
qcOpts.excludeOriginBoxY_mm = [0 1.2]; % <---edit: forbidden origin y-range (mm)
qcOpts.maxStepDy_mm = 0.1; % <---edit: max per-step y-displacement (mm), rejects mis-linked tracks

activationOpts = struct();
activationOpts.minPreJumpArea_px2 = 3; % <---edit: pre-jump area floor (px^2), eliminates sub-pixel noise
activationOpts.minPostJumpArea_px2 = 20; % <---edit: post-jump area floor (px^2), eliminates out-of-focus artifacts
activationOpts.preWindowFrames = 2; % number of frames used to compute pre-jump baseline A_pre
activationOpts.minPrePoints = 1; % minimum finite pre-jump points required to test activation
activationOpts.minInitialGrowthRatio = 2.0; % <---edit: sanity floor for initial growth (aNext/aPre)
activationOpts.sustainedWindowFrames = 3; % <---edit: frames after i+1 to check for sustained growth
activationOpts.largeAreaThreshold_px2 = 100; % <---edit: post-jump area above which shorter sustained window is used
activationOpts.largeAreaSustainedFrames = 2; % <---edit: sustained window frames for large-area activations
activationOpts.sustainedMinRatio = 1.05; % <---edit: per-step growth ratio threshold
activationOpts.sustainedMinDelta_px2 = 5; % <---edit: absolute growth alternative for small bubbles (px^2)
activationOpts.sustainedMinPassingPairs = 2; % <---edit: minimum consecutive pairs that must show growth
activationOpts.sustainedMinPairsAvailable = 2; % <---edit: below this, fall to burst fallback
activationOpts.enableBurstFallback = true; % allow burst detection when sustained window is too short
activationOpts.burstMinGrowthRatio = 3.0; % <---edit: burst fallback: aNext/aPre threshold
activationOpts.burstMinArea_px2 = 50; % <---edit: burst fallback: minimum post-jump area (px^2)
activationOpts.postActivationCheckFrames = 1; % <---edit: number of frames after jump to verify area remains large
activationOpts.postActivationMinArea_px2 = 50; % <---edit: minimum area (px^2) required in each post-activation check frame
activationOpts.rejectProximityMerge = true; % <---edit: reject activations where another spot is nearby at pre-jump frame
activationOpts.mergeProximityRadius_mm = 0; % <---edit: proximity radius (mm) for merge rejection (~13 px)
activationOpts.includeMicrobubbleActivationRescue = true; % <---edit: include activated microbubble-start tracks even when non-strict
activationOpts.microbubbleStartAreaRange_px2 = [1 120]; % <---edit: start-area range (px^2) for green track inclusion
activationOpts.microbubbleRequireOutsideStrictPrimary = true; % <---edit: only add green tracks that are outside strict-blue membership
activationOpts.microbubbleSeedAreaRange_px2 = activationOpts.microbubbleStartAreaRange_px2; % legacy alias
activationOpts.microbubbleRequireNonLeftMoving = activationOpts.microbubbleRequireOutsideStrictPrimary; % legacy alias

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
plotOpts.makeVideoOverlayGifs = true; % <---edit: overlay tracks on source AVI and export GIF
plotOpts.saveVideoOverlayGifs = true; % <---edit
plotOpts.diagnosticGifTrailLength = 10;
plotOpts.diagnosticGifDelayTime = 0.08;
plotOpts.diagnosticGifMaxFrames = localGifMaxFrames;
plotOpts.videoOverlayTrailLength = 10; % <---edit
plotOpts.videoOverlayDelayTime = 0.08; % <---edit
plotOpts.videoOverlayFadeHalfLifeFrames = 30; % <---edit
plotOpts.videoOverlayUseContiguousRange = true; % <---edit
plotOpts.videoOverlayMaxFrames = localGifMaxFrames;
plotOpts.videoOverlayMarkerSize = 15; % <---edit
plotOpts.videoOverlayXLim_mm = [0 4.8]; % <---edit
plotOpts.upstreamSizeXLim_um = []; % <---edit
plotOpts.inceptionImageSize_px = [1280 320]; % [width height]
plotOpts.inceptionXLim_mm = [0 4.8]; % <---edit
plotOpts.inceptionYLim_mm = [0 1.2]; % <---edit
plotOpts.maxActivationsPerCase = 250; % cap activation points shown per case in capped inception plot
plotOpts.makeCappedActivationInceptionPlots = true; % <---edit
plotOpts.makeAllActivationInceptionPlots = true; % <---edit
plotOpts.themes = enabled_plot_themes(plotOpts);

%% ---------------- DEFINE CASES (ONE OR MULTIPLE RE) ----------------
% Required fields per case:
%   name, Re, kD, xmlFile, pixelSize, dt
% If you include one Re only, everything still runs.
% Edit case definitions below for your runs. <---edit

cases = struct([]);

cases(1).name      = "5um";
cases(1).Re        = 95000;
cases(1).kD       = 0.030; % <-- set your k/d
cases(1).xmlFile   = "E:\March Re 90,000 inception data\Processed images\Smooth variation 2\test_nofilter_smoothvar2_48lit.xml"; % <---edit: set your XML path
cases(1).pixelSize = 0.00375009375;  % mm/px
cases(1).dt        = 1/102247;
cases(1).videoFile = "E:\March Re 90,000 inception data\Processed images\Smooth variation 2\Smooth variation 2000.avi"; % <---edit: optional AVI path for overlay GIF export

cases(2).name      = "12um";
cases(2).Re        = 95000;
cases(2).kD       = 0.080;
cases(2).xmlFile   = "E:\March Re 90,000 inception data\Processed images\P10S100\P10S100_48lit.xml";
cases(2).pixelSize = 0.00375009375;
cases(2).dt        = 1/102247;
cases(2).videoFile = "E:\March Re 90,000 inception data\Processed images\results\2000 frames\P10S100 2000.avi";

cases(3).name      = "20um";
cases(3).Re        = 95000;
cases(3).kD       = 0.141;
cases(3).xmlFile   = "E:\March Re 90,000 inception data\Processed images\P10S70\P10S70_48lit.xml";
cases(3).pixelSize = 0.00375009375;
cases(3).dt        = 1/102247;
cases(3).videoFile = "E:\March Re 90,000 inception data\Processed images\results\2000 frames\P10S70 2000.avi";

cases(4).name      = "30um";
cases(4).Re        = 95000;
cases(4).kD       = 0.267;
cases(4).xmlFile   = "E:\March Re 90,000 inception data\Processed images\P10S50\P10S50_48lit.xml";
cases(4).pixelSize = 0.00375009375;
cases(4).dt        = 1/102247;
cases(4).videoFile = "E:\March Re 90,000 inception data\Processed images\results\2000 frames\P10S50 2000.avi";

cases(5).name      = "53um";
cases(5).Re        = 95000;
cases(5).kD       = 0.444;
cases(5).xmlFile   = "E:\March Re 90,000 inception data\Processed images\P10S30\P10S30_48lit.xml";
cases(5).pixelSize = 0.00375009375;
cases(5).dt        = 1/102247;
cases(5).videoFile = "E:\March Re 90,000 inception data\Processed images\results\2000 frames\P10S30 2000.avi";

cases(6).name      = "80um";
cases(6).Re        = 95000;
cases(6).kD       = 0.720;
cases(6).xmlFile   = "E:\March Re 90,000 inception data\Processed images\P10S20\P10S20_48lit.xml";
cases(6).pixelSize = 0.00375009375;
cases(6).dt        = 1/102247;
cases(6).videoFile = "E:\March Re 90,000 inception data\Processed images\results\2000 frames\P10S20 2000.avi";

for ci = 1:numel(cases)
    cases(ci).diagnosticTrackIds = []; % empty = all parsed tracks in this case; otherwise list TRACK_IDs
    if ~isfield(cases, 'videoFile') || isempty(cases(ci).videoFile) || ...
            (~isstring(cases(ci).videoFile) && ~ischar(cases(ci).videoFile)) || ...
            strlength(string(cases(ci).videoFile)) < 1
        cases(ci).videoFile = ""; % <---edit: optional AVI path for overlay GIF export
    end
end

% Example:
% cases(2).diagnosticTrackIds = [123 456];
% cases(1).videoFile = "E:\path\to\case1.avi"; % <---edit

% Example second Reynolds set (uncomment/edit as needed)
% cases(7).name      = "5um";
% cases(7).Re        = 124000;
% cases(7).kD       = 0.0005;
% cases(7).xmlFile   = "<path_to_Re124k_case_5um.xml>";
% cases(7).pixelSize = 0.00375009375;
% cases(7).dt        = 1/102247;

%% ----------- ROI DATA (unwanted area + wall mask) ----------------
% Set to the ROI_throat.mat saved by Testing/ROI_and_throatloader.m.
% Leave as "" to disable all ROI-based filtering and overlays.
roiFile = "E:\March Re 90,000 inception data\Processed images\results\2000 frames\Smooth variation 2000\ROI_throat.mat";  % <---edit

if isstring(roiFile) && strlength(roiFile) > 0 && isfile(roiFile)
    R = load(roiFile);
    roiData = struct();
    roiData.wallMask            = R.ROI_throat.wallMask;
    roiData.unwantedTrackMask   = R.ROI_throat.unwantedTrackMask;
    roiData.throat_xy_px        = R.ROI_throat.throat_xy_px;
    roiData.maskPixelSize       = cases(1).pixelSize;  % same calibration for all cases
    qcOpts.unwantedAreaMask           = roiData.unwantedTrackMask;
    qcOpts.unwantedAreaMaskPixelSize  = roiData.maskPixelSize;
    plotOpts.roiData = roiData;
    fprintf('Loaded ROI data: unwanted mask=%d px, wall mask=%d px\n', ...
        sum(roiData.unwantedTrackMask(:)), sum(roiData.wallMask(:)));
elseif isstring(roiFile) && strlength(roiFile) > 0
    warning('ROI file not found: %s — mask filtering disabled.', roiFile);
end

%% ---- Collapse analysis options -----------------------------------------
% See analyze_collapse_events.m for full signal definition.
collapseOpts.minTrackSpots            = 3;    % <---edit: min spots for a valid track
collapseOpts.collapseMinPeakArea_px2  = 50;   % <---edit: min peak area to qualify (px^2)
collapseOpts.collapseTruncationFactor = 0.6;  % <---edit: final/peak area ratio above which track is flagged as truncated (not a collapse)
collapseOpts.fftNDomFreqs             = 5;    % number of dominant FFT peaks to report
if exist('roiData','var') && isstruct(roiData) && isfield(roiData,'unwantedTrackMask')
    collapseOpts.roiUnwantedMask = roiData.unwantedTrackMask;
    collapseOpts.roiPixelSize    = cases(1).pixelSize;
else
    collapseOpts.roiUnwantedMask = [];
    collapseOpts.roiPixelSize    = 0;
end

%% ---- Void fraction analysis options ------------------------------------
% See analyze_void_fraction.m for signal definition.
% The ROI mask subtracts wall / unwanted-track pixels from the FOV
% denominator so alpha = bubble area / fluid-accessible area (not total frame).
if exist('roiData','var') && isstruct(roiData) && isfield(roiData,'unwantedTrackMask')
    voidFracOpts.roiUnwantedMask  = roiData.unwantedTrackMask;
    voidFracOpts.roiPixelSize     = cases(1).pixelSize;
    voidFracOpts.cameraPixelSize  = cases(1).pixelSize;
else
    voidFracOpts.roiUnwantedMask  = [];
    voidFracOpts.roiPixelSize     = 0;
    voidFracOpts.cameraPixelSize  = cases(1).pixelSize;
end

[cases, selectedCaseIdx, totalCaseCount] = select_cases(cases, caseSelection);
selectedLabels = strings(numel(cases), 1);
for si = 1:numel(cases)
    selectedLabels(si) = sprintf('%d:%s', selectedCaseIdx(si), char(cases(si).name));
end
fprintf('Selected %d/%d case(s): %s\n', numel(cases), totalCaseCount, strjoin(cellstr(selectedLabels), ', '));

%% ---------------- RUN ALL CASES ----------------
summaryRows = table();

allLoc = struct();
allLoc.caseName   = strings(0,1);
allLoc.Re         = nan(0,1);
allLoc.kD         = nan(0,1);
allLoc.pixelSize  = nan(0,1);
allLoc.nActivated = nan(0,1);
allLoc.nInjected  = nan(0,1);
allLoc.inception2x_xy = cell(0,1);  % [x y] activation points on left-moving + microbubble-rescue tracks

trackFigOutDir = fullfile(figDir, "TrackDiagnostics");
if ~isfolder(trackFigOutDir), mkdir(trackFigOutDir); end

gifOutDir = fullfile(resultsDir, "diagnostic gifs");
if plotOpts.saveDiagnosticGifs && ~isfolder(gifOutDir), mkdir(gifOutDir); end

videoGifOutDir = fullfile(resultsDir, "video overlay gifs");
if plotOpts.makeVideoOverlayGifs && plotOpts.saveVideoOverlayGifs && ~isfolder(videoGifOutDir), mkdir(videoGifOutDir); end

allCollapse = struct();
allCollapse.caseName = {};
allCollapse.kD       = [];
allCollapse.Re       = [];
allCollapse.dt       = [];
allCollapse.data     = {};

allVoidFrac = struct();
allVoidFrac.caseName = {};
allVoidFrac.kD       = [];
allVoidFrac.Re       = [];
allVoidFrac.data     = {};

allSize = struct();
allSize.caseName = strings(0,1);
allSize.Re       = nan(0,1);
allSize.kD      = nan(0,1);
allSize.size_eqd = cell(0,1);

gateSummaryRows = table();
cachePolicyTag = build_cache_policy_tag(parserOpts, qcOpts, flowOpts, activationOpts);
runTimer = tic;

for i = 1:numel(cases)
    caseTimer = tic;
    fprintf("\n=== Case %d/%d: %s (Re=%g, k/d=%.6g) ===\n", ...
        i, numel(cases), cases(i).name, cases(i).Re, cases(i).kD);

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
    activationOpts.pixelSize = cases(i).pixelSize;
    metrics = trackmate_case_metrics(out, qcOpts, flowOpts, activationOpts);
    g = metrics.gateStats;
    fprintf(['Strict gate summary: total=%d, injected=%d, activated=%d, ', ...
        'reject(originBox=%d short=%d nonfinite=%d nonmonoT=%d yStep=%d topo=%d flow=%d origin=%d noAct=%d wall=%d)\n'], ...
        g.nTracksTotal, g.nInjected, g.nActivated, ...
        g.nRejectedOriginWindow, g.nRejectedTooShort, g.nRejectedNonFinite, g.nRejectedNonMonotonicTime, ...
        g.nRejectedExcessiveYStep, g.nRejectedTopology, g.nRejectedFlow, g.nRejectedOrigin, ...
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
    allLoc.caseName(end+1,1)   = string(cases(i).name);
    allLoc.Re(end+1,1)         = cases(i).Re;
    allLoc.kD(end+1,1)         = cases(i).kD;
    allLoc.pixelSize(end+1,1)  = cases(i).pixelSize;
    allLoc.nActivated(end+1,1) = g.nActivated;
    allLoc.nInjected(end+1,1)  = g.nInjected;
    allLoc.inception2x_xy{end+1,1} = choose_inception_activation_xy(metrics);

    % Collapse frequency analysis (all tracks, no direction filter)
    collapseResult = analyze_collapse_events(out, cases(i).pixelSize, cases(i).dt, collapseOpts);
    allCollapse.caseName{end+1} = char(cases(i).name);
    allCollapse.kD(end+1)       = cases(i).kD;
    allCollapse.Re(end+1)       = cases(i).Re;
    allCollapse.dt(end+1)       = cases(i).dt;
    allCollapse.data{end+1}     = collapseResult;

    % Void fraction analysis (all detected spots, no filter)
    voidFracOpts.cameraPixelSize = cases(i).pixelSize;
    vfResult = analyze_void_fraction(out, voidFracOpts);
    allVoidFrac.caseName{end+1} = char(cases(i).name);
    allVoidFrac.kD(end+1)       = cases(i).kD;
    allVoidFrac.Re(end+1)       = cases(i).Re;
    allVoidFrac.data{end+1}     = vfResult;

    % Accumulate upstream-size samples for distribution plot
    allSize.caseName(end+1,1) = string(cases(i).name);
    allSize.Re(end+1,1)       = cases(i).Re;
    allSize.kD(end+1,1)      = cases(i).kD;
    allSize.size_eqd{end+1,1} = metrics.upstreamSize_eqd;

    gateRow = table( ...
        string(cases(i).name), cases(i).Re, cases(i).kD, ...
        g.nTracksTotal, g.nInjected, g.nActivated, ...
        g.nRejectedTooShort, g.nRejectedNonFinite, g.nRejectedNonMonotonicTime, ...
        g.nRejectedTopology, g.nRejectedFlow, g.nRejectedOrigin, ...
        g.nRejectedNoActivation, g.nRejectedWallBand, ...
        g.originThreshold, g.xStartMin, g.xStartMax, ...
        'VariableNames', {'Case','Re','kD','nTracksTotal','nInjected','nActivated', ...
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

    % Microbubble (green track) activation events
    nMicroAE = 0;
    if isfield(metrics, 'microbubbleActivationEvent_frame_nonLeft') && ~isempty(metrics.microbubbleActivationEvent_frame_nonLeft)
        nMicroAE = numel(metrics.microbubbleActivationEvent_frame_nonLeft);
    end
    nTotalAE = metrics.strictActivationEventsTotal + nMicroAE;

    % Store per-case summary
    row = table( ...
        string(cases(i).name), cases(i).Re, cases(i).kD, ...
        metrics.nTracksTotal, nValidTracks, ...
        nStrictRecirculationTracks, nStrictActivatedTracks, ...
        metrics.strictTrackFrameExposure, metrics.strictActivationEventsTotal, ...
        nMicroAE, nTotalAE, ...
        strictRecirculationFrac_total, strictActivationFrac_valid, ...
        metrics.A_over_I, metrics.A_over_I_ci_low, metrics.A_over_I_ci_high, A_over_I_err_low, A_over_I_err_high, ...
        tau_mean_val, tau_std_val, tau_sem, nTauVals, elapsed_case_sec, ...
        'VariableNames', {'Case','Re','kD','nTracksTotal','nValidTracks', ...
        'nLeftMovingTracks','nLeftMovingActivated', ...
        'leftMovingFrameExposure','AE_leftMoving', ...
        'AE_microbubble','AE_total', ...
        'leftMovingFrac_total','activationFrac_valid','A_over_I','A_over_I_ci_low','A_over_I_ci_high','A_over_I_err_low','A_over_I_err_high', ...
        'tau_mean','tau_std','tau_sem','nTau','elapsed_case_sec'});

    summaryRows = [summaryRows; row]; %#ok<AGROW>

    % Close any accidental figures
    close all force;
end

% Sort by Re then k/d
summaryRows = sortrows(summaryRows, {'Re','kD'});

% Save summary table
summaryCsv = fullfile(resultsDir, "Summary_tracks.csv");
write_table_csv_compat(summaryRows, summaryCsv);
fprintf("\nSaved: %s\n", summaryCsv);

if ~isempty(gateSummaryRows)
    gateSummaryRows = sortrows(gateSummaryRows, {'Re','kD'});
    gateSummaryCsv = fullfile(resultsDir, "track_gate_summary.csv");
    write_table_csv_compat(gateSummaryRows, gateSummaryCsv);
    fprintf("Saved: %s\n", gateSummaryCsv);
end

%% ---------------- SAVE PLOT DATA (.mat) ----------------
matDir = fullfile(resultsDir, "plot_data_mat");
if ~isfolder(matDir), mkdir(matDir); end
save(fullfile(matDir, "activation_summary_by_case.mat"), 'summaryRows');
save(fullfile(matDir, "inception_locations_by_case.mat"), 'allLoc');
save(fullfile(matDir, "upstream_size_distribution_by_case.mat"), 'allSize');
save(fullfile(matDir, "collapse_analysis_by_case.mat"), 'allCollapse');
save(fullfile(matDir, "void_fraction_by_case.mat"), 'allVoidFrac');
normParams = struct('U_throat_ms', 13.32, 'H_throat_m', 10e-3, ...
    't_conv_s', 10e-3 / 13.32, 'throatHeight_mm', 10);
save(fullfile(matDir, "normalization_parameters.mat"), 'normParams');
fprintf("Saved plot data .mat files to: %s\n", matDir);

%% ---------------- PLOT 1: A/I vs k/d (per Re) ----------------
fitTxtFile = fullfile(resultsDir, "fit_AI_vs_kD_by_Re.txt");
plot_ai_vs_kdh_re(summaryRows, figDir, fitTxtFile, plotOpts);

%% ---------------- PLOT 2: Inception (2x growth) locations ----------------
locFigOutDir = fullfile(figDir, "InceptionLocations");
if ~isfolder(locFigOutDir), mkdir(locFigOutDir); end
plot_inception_locations_by_re(allLoc, locFigOutDir, plotOpts);

%% ---------------- PLOT 2b: A/I based on capped activations ----------------
plot_ai_vs_kd_capped(allLoc, figDir, plotOpts);

%% ---------------- PLOT 3: mean residence time Tau vs k/d ----------------
plot_tau_vs_kdh_re(summaryRows, figDir, plotOpts);

%% ---------------- PLOT 4: upstream moving microbubble size distributions ----------------
distFigOutDir = fullfile(figDir, "UpstreamSizeDistributions");
if ~isfolder(distFigOutDir), mkdir(distFigOutDir); end
plot_upstream_size_distribution_by_re(allSize, distFigOutDir, binSize_phys, plotOpts);

%% ---------------- PLOT 5 & 6: Collapse frequency analysis ----------------
collapseFigDir = fullfile(figDir, "CollapseAnalysis");
if ~isfolder(collapseFigDir), mkdir(collapseFigDir); end
plot_collapse_rate_vs_frame(allCollapse, collapseFigDir, plotOpts);
plot_collapse_power_spectrum(allCollapse, collapseFigDir, plotOpts);
write_collapse_analysis_csv(allCollapse, fullfile(resultsDir, "collapse_analysis.csv"));

%% ---------------- PLOT 7: Void fraction vs k/d ----------------
plot_void_fraction_vs_kd(allVoidFrac, figDir, plotOpts);

if useMatCache && cacheUpdated
    save(cacheFile, 'cacheDB', '-v7.3');
    fprintf("Saved cache DB: %s (%d entries)\n", cacheFile, numel(cacheDB.key));
end

totalElapsedSec = toc(runTimer);
fprintf("Total elapsed time (selected cases): %.2f s (%s)\n", totalElapsedSec, format_elapsed_hms(totalElapsedSec));
fprintf("\nAll done. Results in: %s\n", resultsDir);

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
outCsv = fullfile(resultsDir, sprintf('framewise_counts_%s_Re_%g_kD_%g.csv', caseToken, caseDef.Re, caseDef.kD));
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
if isfield(metrics, 'strictActivationEvent_xy') && ~isempty(metrics.strictActivationEvent_xy)
    xy = metrics.strictActivationEvent_xy;
elseif isfield(metrics, 'activationEvent_xy') && ~isempty(metrics.activationEvent_xy)
    xy = metrics.activationEvent_xy;
elseif isfield(metrics, 'activationEvent_xy_netLeftLegacy') && ~isempty(metrics.activationEvent_xy_netLeftLegacy)
    xy = metrics.activationEvent_xy_netLeftLegacy;
elseif isfield(metrics, 'inception2x_xy') && ~isempty(metrics.inception2x_xy)
    xy = metrics.inception2x_xy;
end
if isfield(metrics, 'microbubbleActivationEvent_xy_nonLeft') && ~isempty(metrics.microbubbleActivationEvent_xy_nonLeft)
    xy = [xy; metrics.microbubbleActivationEvent_xy_nonLeft];
end
if ~isempty(xy)
    xy = unique(xy, 'rows', 'stable');
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
    'pv=%d|pts=%d|pft=%d|bulk=%s|dx=%.6g|neg=%.4g|pos=%.4g|origin=%d|originFrac=%.4g|netRun=%d|netRunDx=%.6g|', ...
    'bandOn=%d|bandX0=%.6g|bandX1=%.6g|bandSteps=%d|bandDx=%.6g|bandAct=%d|bandActX=%d|bandBypassOrigin=%d|', ...
    'spots=%d|gaps=%.4g|rsmc=%d|wall=%d|maxDy=%.6g|oxOn=%d|ox0=%.4g|ox1=%.4g|oy0=%.4g|oy1=%.4g|', ...
    'preA=%.4g|postA=%.4g|pre=%d|minpre=%d|initGR=%.4g|sWin=%d|lgA=%.4g|lgSW=%d|sRat=%.4g|sDel=%.4g|sPass=%d|sAvail=%d|', ...
    'burst=%d|bGR=%.4g|bA=%.4g|postChk=%d|postMinA=%.4g|proxMerge=%d|proxR=%.6g|mbAct=%d|mbLo=%.4g|mbHi=%.4g|mbNonLeft=%d'], ...
    parserOpts.parserVersion, double(logical(parserOpts.parseTrackedSpotsOnly)), double(logical(parserOpts.parseFilteredTracksOnly)), ...
    char(string(flowOpts.bulkDirection)), flowOpts.minNetDxCounterflow_mm, ...
    flowOpts.minNegativeStepFraction, flowOpts.maxPositiveStepFraction, ...
    double(logical(flowOpts.requireRightOrigin)), flowOpts.rightOriginFrac, ...
    flowOpts.netLeftMinConsecutiveNegSteps, flowOpts.netLeftMinConsecutiveNetDx_mm, ...
    double(logical(flowOpts.netLeftBandRescueEnabled)), flowOpts.netLeftBandRescueX_mm(1), flowOpts.netLeftBandRescueX_mm(2), ...
    flowOpts.netLeftBandRescueMinCounterflowSteps, flowOpts.netLeftBandRescueMinCounterflowDx_mm, ...
    double(logical(flowOpts.netLeftBandRescueRequireActivation)), double(logical(flowOpts.netLeftBandRescueRequireActXInBand)), ...
    double(logical(flowOpts.netLeftBandRescueBypassOriginGate)), ...
    qcOpts.minTrackSpots, qcOpts.maxTrackGaps, double(logical(qcOpts.rejectSplitMergeComplex)), ...
    double(logical(qcOpts.wallBandEnabled)), qcOpts.maxStepDy_mm, ...
    double(logical(qcOpts.excludeOriginBoxEnabled)), qcOpts.excludeOriginBoxX_mm(1), qcOpts.excludeOriginBoxX_mm(2), ...
    qcOpts.excludeOriginBoxY_mm(1), qcOpts.excludeOriginBoxY_mm(2), ...
    activationOpts.minPreJumpArea_px2, activationOpts.minPostJumpArea_px2, ...
    activationOpts.preWindowFrames, activationOpts.minPrePoints, ...
    activationOpts.minInitialGrowthRatio, ...
    activationOpts.sustainedWindowFrames, activationOpts.largeAreaThreshold_px2, ...
    activationOpts.largeAreaSustainedFrames, activationOpts.sustainedMinRatio, ...
    activationOpts.sustainedMinDelta_px2, activationOpts.sustainedMinPassingPairs, ...
    activationOpts.sustainedMinPairsAvailable, ...
    double(logical(activationOpts.enableBurstFallback)), ...
    activationOpts.burstMinGrowthRatio, activationOpts.burstMinArea_px2, ...
    activationOpts.postActivationCheckFrames, activationOpts.postActivationMinArea_px2, ...
    double(logical(activationOpts.rejectProximityMerge)), activationOpts.mergeProximityRadius_mm, ...
    double(logical(activationOpts.includeMicrobubbleActivationRescue)), ...
    activationOpts.microbubbleSeedAreaRange_px2(1), activationOpts.microbubbleSeedAreaRange_px2(2), ...
    double(logical(activationOpts.microbubbleRequireNonLeftMoving)));
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
