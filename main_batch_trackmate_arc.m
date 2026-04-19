%% main_batch_trackmate_arc.m
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
% runMode = "local"; % <---edit
runMode = "arc"; % <---edit

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
% On ARC the toggle works the same way — edit localLimitMode above.

% Case selection:
% - "all"            -> run every case below
% - 1                -> run one case by index
% - [1 3 6]          -> run multiple cases by index
% - "5um"            -> run one case by name
% - ["5um","30um"]   -> run multiple cases by name
caseSelection = "all"; % <---edit

% Multi-XML independent-sample auto-detection:
% Keep cases(i).xmlFile pointed at the base XML path or folder path. If it
% points to a folder, the folder name is used as the XML stem. If sibling
% files named <stem>_1.xml, <stem>_2.xml, ... exist, those numbered XMLs
% are used as independent random samples for that case. If no numbered XMLs
% exist, the unnumbered base XML is used as the single combined input.
% Combined results stay in resultsDir; per-case convergence results are
% written under resultsDir/results individual/<caseName>/.
convergenceFrameStep = 2000; % <---edit: cumulative frame interval for convergence CSVs

% Cache parsed outputs (.mat) to avoid re-parsing XML on reruns

useMatCache = true; % <---edit
forceReparse = false; % <---edit

% ----------- ROI DATA (unwanted area + wall mask) ----------------
% Set to the ROI_throat.mat saved by Testing/ROI_and_throatloader.m.
% Upload this file to ARC alongside the XML files. Leave as "" to disable.

roiFile = "/home/kbsanjayvasanth/Tracking dataanlaysis/roi/ROI_throat.mat";  % <---edit: ARC path

% Where to save all results
resultsDir = "/home/kbsanjayvasanth/Tracking dataanlaysis/results/results7"; % <---edit
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
parserOpts.parserVersion = 4;
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
plotOpts.makeTrackDiagnostics = false; % <---edit
plotOpts.saveDiagnosticGifs = false; % <---edit
plotOpts.makeVideoOverlayGifs = false; % <---edit: overlay tracks on source AVI and export GIF
plotOpts.saveVideoOverlayGifs = false; % <---edit
plotOpts.diagnosticGifTrailLength = 10;
plotOpts.diagnosticGifDelayTime = 0.08;
plotOpts.diagnosticGifMaxFrames = localGifMaxFrames;
plotOpts.videoOverlayTrailLength = 10; % <---edit
plotOpts.videoOverlayDelayTime = 0.08; % <---edit
plotOpts.videoOverlayFadeHalfLifeFrames = 30; % <---edit
plotOpts.videoOverlayUseContiguousRange = true; % <---edit
plotOpts.videoOverlayMaxFrames = localGifMaxFrames;
plotOpts.videoOverlayMarkerSize = 26; % <---edit
plotOpts.videoOverlayXLim_mm = [0 4.8]; % <---edit
plotOpts.upstreamSizeXLim_um = []; % <---edit
plotOpts.inceptionImageSize_px = [1280 320]; % [width height]
plotOpts.inceptionXLim_mm = [0 5]; % <---edit
plotOpts.inceptionYLim_mm = [0 1.2]; % <---edit
plotOpts.maxActivationsPerCase = 250; % cap activation points shown per case in capped inception plot
plotOpts.makeCappedActivationInceptionPlots = true; % <---edit
plotOpts.makeAllActivationInceptionPlots = true; % <---edit
plotOpts.breakupDRatioXScale = "log"; % <---edit: x scale for gamma vs d_child/d_parent
plotOpts.breakupDRatioXLim = [0.11, 1.4]; % <---edit: fixed x-window for gamma vs d_child/d_parent
plotOpts.breakupDRatioClipLowPercentile = 0.5; % <---edit: fallback if breakupDRatioXLim is empty
plotOpts.breakupDRatioClipPercentile = 99.5; % <---edit: fallback if breakupDRatioXLim is empty
plotOpts.breakupGammaYLim = [-0.5, 0.5]; % <---edit: fixed y-window for gamma vs d_child/d_parent
plotOpts.breakupGammaYClipPercentile = [1, 99]; % <---edit: robust auto y-window percentile
plotOpts.breakupDRatioMarkerSize = 30; % <---edit: gamma vs d_child/d_parent marker size
plotOpts.breakupDRatioMarkerAlpha = 0.40; % <---edit: gamma vs d_child/d_parent marker transparency
plotOpts.breakupDRatioTrendMaxBins = 12; % <---edit: binned mean trend bins
plotOpts.breakupDRatioTrendMinCount = 5; % <---edit: skip sparse bins in mean trend
plotOpts.breakupARMarkerSize = 30; % <---edit: gamma vs parent AR marker size
plotOpts.breakupARMarkerAlpha = 0.40; % <---edit: gamma vs parent AR marker transparency
plotOpts.breakupARXLim = [0, 4]; % <---edit: fixed x-window for gamma vs parent AR
plotOpts.breakupARGammaYLim = [-0.5, 0.5]; % <---edit: fixed y-window for gamma vs parent AR
plotOpts.breakupARTrendMaxBins = 12; % <---edit: binned mean trend bins for gamma vs parent AR
plotOpts.breakupARTrendMinCount = 5; % <---edit: skip sparse bins in AR mean trend
plotOpts.themes = enabled_plot_themes(plotOpts);

%% ---- Lagrangian acceleration / pressure-gradient proxy options ----------
lagAccelOpts = struct();
lagAccelOpts.sgWindowFrames = 7; % <---edit: Savitzky-Golay local polynomial window
lagAccelOpts.sgPolyOrder = 3; % <---edit: local polynomial order for smoothing/derivatives
lagAccelOpts.triggerWindowFrames = 5; % <---edit: frames immediately before activation
lagAccelOpts.minTriggerSamples = 5; % <---edit: require complete trigger/random windows
lagAccelOpts.minTrackFrames = 12; % <---edit: skip short tracks before differentiating
lagAccelOpts.requireConsecutiveFrames = true; % <---edit: avoid frame-gap derivative artifacts
lagAccelOpts.maxAllowedFrameGap = 1;
lagAccelOpts.excludeSmoothingEdgeFrames = true; % edge derivatives are less reliable
lagAccelOpts.trackPopulation = "strictPrimary"; % strict left-moving/injected tracks
lagAccelOpts.bulkDirection = flowOpts.bulkDirection;
lagAccelOpts.throatHeight_mm = 10;
lagAccelOpts.imageSize_px = plotOpts.inceptionImageSize_px;
lagAccelOpts.xLimNorm = plotOpts.inceptionXLim_mm ./ lagAccelOpts.throatHeight_mm;
lagAccelOpts.xLimNorm(2) = min(lagAccelOpts.xLimNorm(2), 0.5);
lagAccelOpts.yLimNorm = plotOpts.inceptionYLim_mm ./ lagAccelOpts.throatHeight_mm;
lagAccelOpts.fallbackURef_m_s = 13.32;
lagAccelOpts.fallbackDiameter_m = 100e-6;
lagAccelOpts.heatmapGridSize = [25 25];
lagAccelOpts.heatmapStats = ["median", "p90"];
lagAccelOpts.heatmapMinSamplesPerBin = 10;
lagAccelOpts.activationOverlayPerCase = 75;
lagAccelOpts.randomSeed = 42;
lagAccelOpts.minSpearmanN = 10;
lagAccelOpts.makeSanityPlots = true;
lagAccelOpts.nSanityTracks = 5;
lagAccelOpts.makeStationarySanityCheck = true;
lagAccelOpts.stationaryMaxNetDisplacement_px = 2;
lagAccelOpts.stationaryMaxPathLength_px = 5;

%% ---- Proximity activation / neighbor-response options ------------------
proximityActivationOpts = struct();
proximityActivationOpts.maxGamma = 20; % main interaction-zone limit, d/Rmax
proximityActivationOpts.extendedMaxGamma = 40; % saved/plotted supplemental range
proximityActivationOpts.gammaBins = [0 2 5 10 20 40];
proximityActivationOpts.similarSizeRatioRange = [0.5 2.0]; % neighbor R0 / primary R0
proximityActivationOpts.baselineFrames = 5;
proximityActivationOpts.postCollapseFrames = 5;
proximityActivationOpts.secondaryStartLagFrames = 1;
proximityActivationOpts.secondaryPostCollapseFrames = 5;
proximityActivationOpts.referenceFrameTolerance = 2;
proximityActivationOpts.minNeighborFramesInWindow = 3;
proximityActivationOpts.minCorrelationFrames = 5;
proximityActivationOpts.smoothingWindowFrames = 5;
proximityActivationOpts.smoothingPolyOrder = 2;
proximityActivationOpts.maxCrossCorrelationLagFrames = 5;
proximityActivationOpts.useOnlyActivatedPrimaryCollapses = true;
proximityActivationOpts.requireBasicValidNeighbors = true;
proximityActivationOpts.requireNeighborBaselineBeforePeak = true;
proximityActivationOpts.excludeAlreadyActivatedNeighbors = true;
proximityActivationOpts.randomSeed = 42;

%% ---------------- DEFINE CASES (ONE OR MULTIPLE RE) ----------------
% Required fields per case:
%   name, Re, kD, xmlFile, pixelSize, dt
%   xmlFile may be a base XML path or a folder path. If it points to a
%   folder, the folder name is used as the XML stem for numbered sample detection.
% Optional:
%   videoFile = base AVI path or folder path. If it points to a folder,
%   the folder name is used as the video stem. If sibling files named
%   <stem>_1.avi, <stem>_2.avi, ... exist, overlay GIFs use those per XML sample.
% If you include one Re only, everything still runs.
% Edit case definitions below for your runs. <---edit

cases = struct([]);

cases(1).name      = "5um";
cases(1).Re        = 95000;
cases(1).kD       = 0.030; % <-- set your k/d
cases(1).xmlFile   = "/home/kbsanjayvasanth/Tracking dataanlaysis/batch_data/Results/Smooth"; % <-- set XML path
cases(1).pixelSize = 0.00375009375;  % mm/px
cases(1).dt        = 1/102247;
cases(1).videoFile = "/home/kbsanjayvasanth/Tracking dataanlaysis/videos/Smooth variation 2000.avi";

cases(2).name      = "12um";
cases(2).Re        = 95000;
cases(2).kD       = 0.080;
cases(2).xmlFile   = "/home/kbsanjayvasanth/Tracking dataanlaysis/batch_data/Results/P10S100";
cases(2).pixelSize = 0.00375009375;
cases(2).dt        = 1/102247;
cases(2).videoFile = "/home/kbsanjayvasanth/Tracking dataanlaysis/videos/P10S100 2000.avi";

cases(3).name      = "20um";
cases(3).Re        = 95000;
cases(3).kD       = 0.141;
cases(3).xmlFile   = "/home/kbsanjayvasanth/Tracking dataanlaysis/batch_data/Results/P10S70";
cases(3).pixelSize = 0.00375009375;
cases(3).dt        = 1/102247;
cases(3).videoFile = "/home/kbsanjayvasanth/Tracking dataanlaysis/videos/P10S70 2000.avi";

cases(4).name      = "30um";
cases(4).Re        = 95000;
cases(4).kD       = 0.267;
cases(4).xmlFile   = "/home/kbsanjayvasanth/Tracking dataanlaysis/batch_data/Results/P10S50";
cases(4).pixelSize = 0.00375009375;
cases(4).dt        = 1/102247;
cases(4).videoFile = "/home/kbsanjayvasanth/Tracking dataanlaysis/videos/P10S50 2000.avi";

cases(5).name      = "53um";
cases(5).Re        = 95000;
cases(5).kD       = 0.444;
cases(5).xmlFile   = "/home/kbsanjayvasanth/Tracking dataanlaysis/batch_data/Results/P10S30";
cases(5).pixelSize = 0.00375009375;
cases(5).dt        = 1/102247;
cases(5).videoFile = "/home/kbsanjayvasanth/Tracking dataanlaysis/videos/P10S30 2000.avi";

cases(6).name      = "80um";
cases(6).Re        = 95000;
cases(6).kD       = 0.720;
cases(6).xmlFile   = "/home/kbsanjayvasanth/Tracking dataanlaysis/batch_data/Results/P10S20";
cases(6).pixelSize = 0.00375009375;
cases(6).dt        = 1/102247;
cases(6).videoFile = "/home/kbsanjayvasanth/Tracking dataanlaysis/videos/P10S20 2000.avi";

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
% cases(1).videoFile = "/home/.../case1.avi"; % <---edit

% Example second Reynolds set (uncomment/edit as needed)
% cases(7).name      = "5um";
% cases(7).Re        = 124000;
% cases(7).kD       = 0.0005;
% cases(7).xmlFile   = "<path_to_Re124k_case_5um.xml>";
% cases(7).pixelSize = 0.00375009375;
% cases(7).dt        = 1/102247;

%% ----------- ROI DATA (unwanted area + wall mask) ----------------
% roiFile is defined near the top of this script — edit it there.

if isstring(roiFile) && strlength(roiFile) > 0 && isfile(roiFile)
    R = load(roiFile);
    roiData = struct();
    roiData.wallMask            = R.ROI_throat.wallMask;
    roiData.unwantedTrackMask   = R.ROI_throat.unwantedTrackMask;
    roiData.throat_xy_px        = R.ROI_throat.throat_xy_px;
    roiData.maskPixelSize       = cases(1).pixelSize;
    qcOpts.unwantedAreaMask           = roiData.unwantedTrackMask;
    qcOpts.unwantedAreaMaskPixelSize  = roiData.maskPixelSize;
    plotOpts.roiData = roiData;
    fprintf('Loaded ROI data: unwanted mask=%d px, wall mask=%d px\n', ...
        sum(roiData.unwantedTrackMask(:)), sum(roiData.wallMask(:)));
elseif isstring(roiFile) && strlength(roiFile) > 0
    warning('ROI file not found: %s — mask filtering disabled.', roiFile);
end

%% ---- Collapse analysis options -----------------------------------------
collapseOpts.minTrackSpots            = 3;    % <---edit
collapseOpts.collapseMinPeakArea_px2  = 50;   % <---edit
collapseOpts.collapseTruncationFactor = 0.6;  % <---edit
collapseOpts.fftNDomFreqs             = 5;
if exist('roiData','var') && isstruct(roiData) && isfield(roiData,'unwantedTrackMask')
    collapseOpts.roiUnwantedMask = roiData.unwantedTrackMask;
    collapseOpts.roiPixelSize    = cases(1).pixelSize;
else
    collapseOpts.roiUnwantedMask = [];
    collapseOpts.roiPixelSize    = 0;
end

%% ---- Collapse recirculation attribution options ------------------------
% A collapse-generated microbubble is a small-start track born near a
% qualified collapse location within this short frame window.
collapseRecirculationOpts.minFrameLag = 0; % <---edit: first allowed frame after collapse
collapseRecirculationOpts.maxFrameLag = 3; % <---edit: last allowed frame after collapse
collapseRecirculationOpts.maxDistance_mm = 0.15; % <---edit: spatial match radius around collapse
collapseRecirculationOpts.microbubbleStartAreaRange_px2 = activationOpts.microbubbleStartAreaRange_px2;
collapseRecirculationOpts.requireBasicValid = true;

%% ---- Void fraction analysis options ------------------------------------
if exist('roiData','var') && isstruct(roiData) && isfield(roiData,'unwantedTrackMask')
    voidFracOpts.roiUnwantedMask  = roiData.unwantedTrackMask;
    voidFracOpts.roiPixelSize     = cases(1).pixelSize;
    voidFracOpts.cameraPixelSize  = cases(1).pixelSize;
else
    voidFracOpts.roiUnwantedMask  = [];
    voidFracOpts.roiPixelSize     = 0;
    voidFracOpts.cameraPixelSize  = cases(1).pixelSize;
end

%% ---- Breakup analysis options ------------------------------------------
breakupOpts.aspectRatioMin       = 1.5;    % <---edit: parent ellipse aspect ratio threshold (lowest AR of interest; post-filtered per threshold below)
breakupOpts.arThresholds         = [1.5, 2.0, 3.0, 4.0]; % <---edit: AR thresholds for breakup k/d and AR plots & .mat files
breakupOpts.childAreaMin_px2     = 100.0;  % <---edit: min child spot area (px^2)
breakupOpts.dRoughnessSpacing_mm = 0.384;  % <---edit: roughness spacing d (mm) for gamma normalisation

%% ---- Growth/collapse rate analysis options ----------------------------
growthCollapseOpts.U_m_s = 13.32; % <---edit: normalizing U for Re = 95,000
growthCollapseOpts.aspectRatioEdges = [0, 2, 5, Inf]; % AR groups: <2, 2-5, >=5
growthCollapseOpts.growthTrackMode = "strictActivated"; % activated upstream-moving tracks
growthCollapseOpts.minTrackSpots = collapseOpts.minTrackSpots;
growthCollapseOpts.collapseMinPeakArea_px2 = collapseOpts.collapseMinPeakArea_px2;
growthCollapseOpts.collapseTruncationFactor = collapseOpts.collapseTruncationFactor;
growthCollapseOpts.roiUnwantedMask = collapseOpts.roiUnwantedMask;
growthCollapseOpts.roiPixelSize = collapseOpts.roiPixelSize;

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
allLoc.leftMovingActivation_xy = cell(0,1);  % [x y] strict left-moving activation points only

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
allCollapse.pixelSize = [];
allCollapse.data     = {};

allCollapseRecirculation = repmat(struct('caseName',"", 'Re',NaN, 'kD',NaN, 'data',[]), 0, 1);

allGrowthCollapse = struct();
allGrowthCollapse.caseName = strings(0,1);
allGrowthCollapse.kD = nan(0,1);
allGrowthCollapse.Re = nan(0,1);
allGrowthCollapse.U_m_s = nan(0,1);
allGrowthCollapse.data = cell(0,1);

allLagAccel = [];
allProximityActivation = [];

allVoidFrac = struct();
allVoidFrac.caseName = {};
allVoidFrac.kD       = [];
allVoidFrac.Re       = [];
allVoidFrac.data     = {};

allBreakup = repmat(struct('caseName', "", 'Re', NaN, 'kD', 0, 'events', []), 0, 1);

allSize = struct();
allSize.caseName = strings(0,1);
allSize.Re       = nan(0,1);
allSize.kD       = nan(0,1);
allSize.size_eqd = cell(0,1);

gateSummaryRows = table();

individualInstanceSummaryRows = table();
individualConvergenceRows = table();

cachePolicyTag = build_cache_policy_tag(parserOpts, qcOpts, flowOpts, activationOpts);
runTimer = tic;

for i = 1:numel(cases)
    caseTimer = tic;
    fprintf("\n=== Case %d/%d: %s (Re=%g, k/d=%.6g) ===\n", ...
        i, numel(cases), cases(i).name, cases(i).Re, cases(i).kD);

    [xmlFiles_i, outChunks, cacheDB, cacheUpdated] = load_case_chunk_outputs( ...
        cases(i), maxTracksToParse, parserOpts, useMatCache, forceReparse, ...
        cacheDB, cacheUpdated, cachePolicyTag, ~isArc);
    nXmlChunks_i = numel(xmlFiles_i);
    hasMultipleXmlSamples = nXmlChunks_i > 1;
    if hasMultipleXmlSamples
        fprintf("  Found %d independent XML sample(s) for pooled analysis.\n", nXmlChunks_i);
    end

    % Pool independent XML samples without linking tracks across files. The
    % merge only offsets IDs/frame axes so counts can be evaluated together.
    out = merge_parsed_outputs(outChunks, xmlFiles_i);
    chunkInfo = out.meta.chunkInfo;
    chunkVideoFiles_i = {};
    if plotOpts.makeVideoOverlayGifs
        chunkVideoFiles_i = detect_chunk_video_files(cases(i).videoFile);
        if hasMultipleXmlSamples
            if isempty(chunkVideoFiles_i)
                fprintf("  Video overlay: no numbered AVI files found for this multi-XML case.\n");
            elseif numel(chunkVideoFiles_i) ~= nXmlChunks_i
                fprintf("  Video overlay: found %d AVI file(s) for %d XML sample(s); exporting first %d overlay GIF(s).\n", ...
                    numel(chunkVideoFiles_i), nXmlChunks_i, min(numel(chunkVideoFiles_i), nXmlChunks_i));
            end
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

    % Accumulate inception (2x growth) locations for plotting
    allLoc.caseName(end+1,1)   = string(cases(i).name);
    allLoc.Re(end+1,1)         = cases(i).Re;
    allLoc.kD(end+1,1)         = cases(i).kD;
    allLoc.pixelSize(end+1,1)  = cases(i).pixelSize;
    allLoc.nActivated(end+1,1) = g.nActivated;
    allLoc.nInjected(end+1,1)  = g.nInjected;
    allLoc.inception2x_xy{end+1,1} = choose_inception_activation_xy(metrics);
    allLoc.leftMovingActivation_xy{end+1,1} = choose_left_moving_activation_xy(metrics);

    % Lagrangian acceleration / pressure-gradient proxy analysis
    lagAccelResult = analyze_lagrangian_acceleration(out, metrics, cases(i), lagAccelOpts);
    allLagAccel(end+1,1) = lagAccelResult; %#ok<AGROW>

    % Collapse frequency analysis (all tracks, no direction filter)
    collapseResult = analyze_collapse_events(out, cases(i).pixelSize, cases(i).dt, collapseOpts);
    allCollapse.caseName{end+1} = char(cases(i).name);
    allCollapse.kD(end+1)       = cases(i).kD;
    allCollapse.Re(end+1)       = cases(i).Re;
    allCollapse.dt(end+1)       = cases(i).dt;
    allCollapse.pixelSize(end+1) = cases(i).pixelSize;
    allCollapse.data{end+1}     = collapseResult;

    collapseRecirculationResult = analyze_collapse_recirculation( ...
        out, metrics, collapseResult, cases(i), collapseRecirculationOpts);
    allCollapseRecirculation(end+1,1) = struct('caseName', string(cases(i).name), ...
        'Re', cases(i).Re, 'kD', cases(i).kD, 'data', collapseRecirculationResult);

    proximityActivationResult = analyze_proximity_activation( ...
        outChunks, xmlFiles_i, cases(i), qcOpts, flowOpts, activationOpts, collapseOpts, proximityActivationOpts);
    allProximityActivation(end+1,1) = proximityActivationResult; %#ok<AGROW>

    growthCollapseResult = analyze_growth_collapse_rates(out, metrics, cases(i), growthCollapseOpts);
    allGrowthCollapse.caseName(end+1,1) = string(cases(i).name);
    allGrowthCollapse.kD(end+1,1)       = cases(i).kD;
    allGrowthCollapse.Re(end+1,1)       = cases(i).Re;
    allGrowthCollapse.U_m_s(end+1,1)    = growthCollapseOpts.U_m_s;
    allGrowthCollapse.data{end+1,1}     = growthCollapseResult;

    % Void fraction analysis (all detected spots, no filter)
    voidFracOpts.cameraPixelSize = cases(i).pixelSize;
    vfResult = analyze_void_fraction(out, voidFracOpts);
    allVoidFrac.caseName{end+1} = char(cases(i).name);
    allVoidFrac.kD(end+1)       = cases(i).kD;
    allVoidFrac.Re(end+1)       = cases(i).Re;
    allVoidFrac.data{end+1}     = vfResult;

    % Breakup event analysis (elongated-parent split events)
    roiArg = [];
    if exist('roiData','var') && isstruct(roiData), roiArg = roiData; end
    breakupEventsByChunk = cell(nXmlChunks_i, 1);
    for c = 1:nXmlChunks_i
        breakupEventsByChunk{c} = analyze_breakup_events(xmlFiles_i{c}, cases(i).pixelSize, ...
            'roiData',                roiArg, ...
            'aspectRatioMin',         breakupOpts.aspectRatioMin, ...
            'childAreaMin_px2',       breakupOpts.childAreaMin_px2, ...
            'dRoughnessSpacing_mm',   breakupOpts.dRoughnessSpacing_mm);
    end
    bkEvents = combine_breakup_event_chunks(breakupEventsByChunk, chunkInfo);
    allBreakup(end+1,1) = struct('caseName', string(cases(i).name), ...
        'Re', cases(i).Re, 'kD', cases(i).kD, 'events', bkEvents);

    % Accumulate upstream-size samples for distribution plot
    allSize.caseName(end+1,1) = string(cases(i).name);
    allSize.Re(end+1,1)       = cases(i).Re;
    allSize.kD(end+1,1)       = cases(i).kD;
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

    if plotOpts.makeVideoOverlayGifs && ~hasMultipleXmlSamples
        save_video_overlay_gif_from_avi(cases(i), metrics, videoGifOutDir, plotOpts);
    end

    elapsed_case_sec = toc(caseTimer);
    fprintf('Case elapsed: %.2f s (%s)\n', elapsed_case_sec, format_elapsed_hms(elapsed_case_sec));

    % Store per-case summary
    row = make_activation_summary_row(cases(i), metrics, elapsed_case_sec);

    summaryRows = [summaryRows; row]; %#ok<AGROW>

    % ---- Independent XML-instance convergence for this case --------------
    [caseInstanceRows, caseConvergenceRows] = write_case_independent_convergence( ...
        resultsDir, cases(i), outChunks, xmlFiles_i, row, ...
        convergenceFrameStep, qcOpts, flowOpts, activationOpts);
    individualInstanceSummaryRows = append_table_compat(individualInstanceSummaryRows, caseInstanceRows);
    individualConvergenceRows = append_table_compat(individualConvergenceRows, caseConvergenceRows);

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

collapseRecirculationDir = fullfile(resultsDir, "collapse recirculation");
if ~isfolder(collapseRecirculationDir), mkdir(collapseRecirculationDir); end
write_collapse_recirculation_csv(allCollapseRecirculation, ...
    fullfile(collapseRecirculationDir, "collapse_recirculation_summary.csv"));
plot_collapse_recirculation_vs_kd(allCollapseRecirculation, collapseRecirculationDir, plotOpts);

%% ---------------- PLOT 1: A/I vs k/d (per Re) ----------------
%% ---------------- SAVE PLOT DATA (.mat) ----------------
matDir = fullfile(resultsDir, "plot_data_mat");
if ~isfolder(matDir), mkdir(matDir); end
save(fullfile(matDir, "activation_summary_by_case.mat"), 'summaryRows');
save(fullfile(matDir, "inception_locations_by_case.mat"), 'allLoc');
save(fullfile(matDir, "upstream_size_distribution_by_case.mat"), 'allSize');
save(fullfile(matDir, "collapse_analysis_by_case.mat"), 'allCollapse');
save(fullfile(matDir, "collapse_recirculation_by_case.mat"), 'allCollapseRecirculation');
save(fullfile(matDir, "growth_collapse_rate_by_case.mat"), 'allGrowthCollapse');
save(fullfile(matDir, "lagrangian_acceleration_by_case.mat"), 'allLagAccel', 'lagAccelOpts', '-v7.3');
save(fullfile(matDir, "proximity_activation_by_case.mat"), 'allProximityActivation', 'proximityActivationOpts', '-v7.3');
save(fullfile(matDir, "void_fraction_by_case.mat"), 'allVoidFrac');
% breakup_analysis_by_case.mat saved later with per-AR-threshold variants
normParams = struct('U_throat_ms', 13.32, 'H_throat_m', 10e-3, ...
    't_conv_s', 10e-3 / 13.32, 'throatHeight_mm', 10, ...
    'growthCollapse_U_ms', growthCollapseOpts.U_m_s, ...
    'growthCollapse_aspectRatioEdges', growthCollapseOpts.aspectRatioEdges);
save(fullfile(matDir, "normalization_parameters.mat"), 'normParams');
fprintf("Saved plot data .mat files to: %s\n", matDir);

%% ---------------- LAGRANGIAN ACCELERATION / PRESSURE PROXY --------------
lagAccelFigDir = fullfile(figDir, "lagrangian acceleration");
if ~isfolder(lagAccelFigDir), mkdir(lagAccelFigDir); end
lagAccelSummaryRows = lagrangian_acceleration_to_table(allLagAccel);
write_table_csv_compat(lagAccelSummaryRows, fullfile(lagAccelFigDir, "lagrangian_acceleration_summary.csv"));
save(fullfile(lagAccelFigDir, "lagrangian_acceleration_by_case.mat"), 'allLagAccel', 'lagAccelOpts', '-v7.3');
plot_lagrangian_acceleration_analysis(allLagAccel, lagAccelFigDir, plotOpts, lagAccelOpts);

%% ---------------- PROXIMITY ACTIVATION / NEIGHBOR RESPONSE --------------
proximityFigDir = fullfile(figDir, "Proximity activation");
if ~isfolder(proximityFigDir), mkdir(proximityFigDir); end
[proximityPairRows, proximityEventRows, proximityBinnedRows, proximitySummaryRows] = ...
    proximity_activation_to_tables(allProximityActivation);
if ~isempty(proximityPairRows) && height(proximityPairRows) > 0
    write_table_csv_compat(proximityPairRows, fullfile(proximityFigDir, "proximity_activation_pairs.csv"));
end
if ~isempty(proximityEventRows) && height(proximityEventRows) > 0
    write_table_csv_compat(proximityEventRows, fullfile(proximityFigDir, "proximity_activation_events.csv"));
end
if ~isempty(proximityBinnedRows) && height(proximityBinnedRows) > 0
    write_table_csv_compat(proximityBinnedRows, fullfile(proximityFigDir, "proximity_activation_binned_stats.csv"));
end
if ~isempty(proximitySummaryRows) && height(proximitySummaryRows) > 0
    write_table_csv_compat(proximitySummaryRows, fullfile(proximityFigDir, "proximity_activation_summary.csv"));
end
save(fullfile(proximityFigDir, "proximity_activation_by_case.mat"), ...
    'allProximityActivation', 'proximityActivationOpts', ...
    'proximityPairRows', 'proximityEventRows', 'proximityBinnedRows', 'proximitySummaryRows', '-v7.3');
plot_proximity_activation_analysis(allProximityActivation, proximityFigDir, plotOpts, proximityActivationOpts);

fitTxtFile = fullfile(resultsDir, "fit_AI_vs_kD_by_Re.txt");
plot_ai_vs_kdh_re(summaryRows, figDir, fitTxtFile, plotOpts);

%% ---------------- PLOT 2: Inception (2x growth) locations ----------------
locFigOutDir = fullfile(figDir, "InceptionLocations");
if ~isfolder(locFigOutDir), mkdir(locFigOutDir); end
plot_inception_locations_by_re(allLoc, locFigOutDir, plotOpts);
leftMovingLocPlotOpts = plotOpts;
leftMovingLocPlotOpts.inceptionLocationField = 'leftMovingActivation_xy';
leftMovingLocPlotOpts.inceptionLocationOutputStem = 'LeftMovingActivation_locations';
leftMovingLocPlotOpts.inceptionLocationWarningLabel = 'left-moving activation points';
leftMovingLocPlotOpts.inceptionLocationPlotDimensional = false;
leftMovingLocPlotOpts.inceptionLocationPlotNormalized = true;
plot_inception_locations_by_re(allLoc, locFigOutDir, leftMovingLocPlotOpts);

%% ---------------- PLOT 2b: A/I based on capped activations ----------------
plot_ai_vs_kd_capped(allLoc, figDir, plotOpts);

%% ---------------- PLOT 3: mean residence time Tau vs k/d ----------------
plot_tau_vs_kdh_re(summaryRows, figDir, plotOpts);
plot_tau_velocity_colormap(summaryRows, figDir, plotOpts, matDir);

%% ---------------- PLOT 4: upstream moving microbubble size distributions ----------------
distFigOutDir = fullfile(figDir, "UpstreamSizeDistributions");
if ~isfolder(distFigOutDir), mkdir(distFigOutDir); end
plot_upstream_size_distribution_by_re(allSize, distFigOutDir, binSize_phys, plotOpts);

%% ---------------- PLOT 5 & 6: Collapse frequency analysis ----------------
collapseFigDir = fullfile(figDir, "CollapseAnalysis");
if ~isfolder(collapseFigDir), mkdir(collapseFigDir); end
plot_collapse_rate_vs_frame(allCollapse, collapseFigDir, plotOpts);
plot_collapse_rate_vs_kd(allCollapse, collapseFigDir, plotOpts);
plot_collapse_size_distribution(allCollapse, collapseFigDir, plotOpts);
write_collapse_analysis_csv(allCollapse, fullfile(resultsDir, "collapse_analysis.csv"));

%% ---------------- PLOT 6b: Growth/collapse axial rate PDFs ----------------
growthCollapseFigDir = fullfile(figDir, "Growth and collapse rate");
if ~isfolder(growthCollapseFigDir), mkdir(growthCollapseFigDir); end
plot_growth_collapse_rate_pdf(allGrowthCollapse, growthCollapseFigDir, plotOpts, matDir);

%% ---------------- PLOT 7: Void fraction vs k/d ----------------
plot_void_fraction_vs_kd(allVoidFrac, figDir, plotOpts);

%% ---------------- PLOT 8: Breakup gamma vs d_child/d_parent ----------------
breakupFigDir = fullfile(figDir, "BreakupAnalysis");
if ~isfolder(breakupFigDir), mkdir(breakupFigDir); end

% Save the full breakup dataset. Plots are split by Re when multiple
% Reynolds numbers are present, so identical k/d values do not pool.
save(fullfile(matDir, "breakup_analysis_by_case.mat"), 'allBreakup');
write_breakup_analysis_csv(allBreakup, fullfile(resultsDir, "breakup_events.csv"));
plot_breakup_analysis_by_re(allBreakup, breakupFigDir, plotOpts, matDir, breakupOpts.arThresholds);

%% ================ INDEPENDENT-SAMPLE CONVERGENCE RESULTS ================
if ~isempty(individualConvergenceRows) && height(individualConvergenceRows) > 0
    individualRootDir = fullfile(resultsDir, "results individual");
    if ~isfolder(individualRootDir), mkdir(individualRootDir); end
    write_table_csv_compat(individualConvergenceRows, ...
        fullfile(individualRootDir, "convergence_summary_all_cases.csv"));
end
if ~isempty(individualInstanceSummaryRows) && height(individualInstanceSummaryRows) > 0
    individualRootDir = fullfile(resultsDir, "results individual");
    if ~isfolder(individualRootDir), mkdir(individualRootDir); end
    write_table_csv_compat(individualInstanceSummaryRows, ...
        fullfile(individualRootDir, "xml_instance_summary_all_cases.csv"));
end

if useMatCache && cacheUpdated
    save(cacheFile, 'cacheDB', '-v7.3');
    fprintf("Saved cache DB: %s (%d entries)\n", cacheFile, numel(cacheDB.key));
end

totalElapsedSec = toc(runTimer);
fprintf("Total elapsed time (selected cases): %.2f s (%s)\n", totalElapsedSec, format_elapsed_hms(totalElapsedSec));
fprintf("\nAll done. Results in: %s\n", resultsDir);

function frameTbl = write_framewise_case_csv(resultsDir, caseDef, metrics)
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

function summaryRow = make_activation_summary_row(caseDef, metrics, elapsedCaseSec)
nValidTracks = metrics.nBasicValidTracks;
nStrictRecirculationTracks = metrics.nStrictPrimaryTracks;
nStrictActivatedTracks = metrics.nStrictActivatedTracks;

leftMovingActivated_pct = NaN;
if nStrictRecirculationTracks > 0
    leftMovingActivated_pct = 100 * nStrictActivatedTracks / nStrictRecirculationTracks;
end

nMicroAE = 0;
if isfield(metrics, 'microbubbleActivationEvent_frame_nonLeft') && ~isempty(metrics.microbubbleActivationEvent_frame_nonLeft)
    nMicroAE = numel(metrics.microbubbleActivationEvent_frame_nonLeft);
end
nTotalAE = metrics.strictActivationEventsTotal + nMicroAE;

AE_leftMoving_pct = NaN;
if nTotalAE > 0
    AE_leftMoving_pct = 100 * metrics.strictActivationEventsTotal / nTotalAE;
end

strictRecirculationFrac_total = nStrictRecirculationTracks / max(metrics.nTracksTotal, 1);
strictActivationFrac_valid = nStrictActivatedTracks / max(nValidTracks, 1);
nTauVals = numel(metrics.tau_values);
tau_sem = NaN;
if nTauVals > 1
    tau_sem = metrics.tau_std / sqrt(nTauVals);
end

A_over_I_err_low = metrics.A_over_I - metrics.A_over_I_ci_low;
A_over_I_err_high = metrics.A_over_I_ci_high - metrics.A_over_I;

summaryRow = table( ...
    string(caseDef.name), caseDef.Re, caseDef.kD, ...
    metrics.nTracksTotal, nValidTracks, ...
    nStrictRecirculationTracks, nStrictActivatedTracks, leftMovingActivated_pct, ...
    metrics.strictTrackFrameExposure, metrics.strictActivationEventsTotal, ...
    nMicroAE, nTotalAE, AE_leftMoving_pct, ...
    strictRecirculationFrac_total, strictActivationFrac_valid, ...
    metrics.A_over_I, metrics.A_over_I_ci_low, metrics.A_over_I_ci_high, A_over_I_err_low, A_over_I_err_high, ...
    metrics.tau_mean, metrics.tau_std, tau_sem, nTauVals, ...
    metrics.activatedUpstreamVelocity_mean_m_s, metrics.activatedUpstreamVelocity_std_m_s, ...
    metrics.activatedUpstreamVelocity_sem_m_s, metrics.activatedUpstreamVelocity_n, ...
    elapsedCaseSec, ...
    'VariableNames', {'Case','Re','kD','nTracksTotal','nValidTracks', ...
    'nLeftMovingTracks','nLeftMovingActivated','leftMovingActivated_pct', ...
    'leftMovingFrameExposure','AE_leftMoving', ...
    'AE_microbubble','AE_total','AE_leftMoving_pct', ...
    'leftMovingFrac_total','activationFrac_valid','A_over_I','A_over_I_ci_low','A_over_I_ci_high','A_over_I_err_low','A_over_I_err_high', ...
    'tau_mean','tau_std','tau_sem','nTau', ...
    'activatedVelocity_mean_m_s','activatedVelocity_std_m_s','activatedVelocity_sem_m_s','activatedVelocity_n', ...
    'elapsed_case_sec'});
end

function [instanceRows, convergenceRows] = write_case_independent_convergence( ...
    resultsDir, caseDef, outChunks, xmlFiles, finalSummaryRow, frameStep, ...
    qcOpts, flowOpts, activationOpts)

instanceRows = table();
convergenceRows = table();
if isempty(outChunks)
    return;
end

caseToken = sanitize_case_token(caseDef.name);
caseDir = fullfile(resultsDir, "results individual", caseToken);
if ~isfolder(caseDir), mkdir(caseDir); end

[instanceRows, convergenceRows, totalFrames] = build_case_independent_convergence_tables( ...
    caseDef, outChunks, xmlFiles, frameStep, qcOpts, flowOpts, activationOpts);

finalMeta = table(frameStep, totalFrames, numel(outChunks), ...
    'VariableNames', {'ConvergenceFrameStep','TotalSourceFrames','XMLInstanceCount'});
write_table_csv_compat([finalMeta, finalSummaryRow], fullfile(caseDir, "converged_summary.csv"));

if ~isempty(instanceRows) && height(instanceRows) > 0
    write_table_csv_compat(instanceRows, fullfile(caseDir, "xml_instance_summary.csv"));
end
if ~isempty(convergenceRows) && height(convergenceRows) > 0
    write_table_csv_compat(convergenceRows, fullfile(caseDir, "convergence_summary.csv"));
end

fprintf('Independent sample convergence for %s: %d XML file(s), %d source frame(s), step=%d -> %s\n', ...
    char(string(caseDef.name)), numel(outChunks), totalFrames, frameStep, caseDir);
end

function [instanceRows, convergenceRows, totalFrames] = build_case_independent_convergence_tables( ...
    caseDef, outChunks, xmlFiles, frameStep, qcOpts, flowOpts, activationOpts)

instanceRows = table();
convergenceRows = table();
totalFrames = estimate_total_source_frames(outChunks);

activationOpts.pixelSize = caseDef.pixelSize;

for c = 1:numel(outChunks)
    outOne = outChunks{c};
    [frameMin, frameMax, frameCount] = parsed_output_frame_span(outOne);
    if isempty(outOne) || ~isstruct(outOne) || ~isfield(outOne, 'trajectories') || isempty(outOne.trajectories)
        continue;
    end

    metricsOne = trackmate_case_metrics(outOne, qcOpts, flowOpts, activationOpts);
    summaryOne = make_activation_summary_row(caseDef, metricsOne, 0);
    metaOne = table(c, string(xmlFiles{c}), frameMin, frameMax, frameCount, ...
        'VariableNames', {'XMLInstance','XMLFile','SourceFrameMin','SourceFrameMax','SourceFrameCount'});
    instanceRows = append_table_compat(instanceRows, [metaOne, summaryOne]);
end

if ~(isfinite(frameStep) && frameStep > 0 && totalFrames > 0)
    return;
end

targets = (frameStep:frameStep:totalFrames).';
if isempty(targets) || targets(end) ~= totalFrames
    targets(end+1,1) = totalFrames;
end

previousSummary = table();
for bi = 1:numel(targets)
    targetFrames = targets(bi);
    [outCum, usedFrames, instancesUsed, lastInstance, lastInstanceFrames] = ...
        cumulative_output_by_source_frames(outChunks, xmlFiles, targetFrames);
    if isempty(outCum) || usedFrames <= 0
        continue;
    end

    metricsCum = trackmate_case_metrics(outCum, qcOpts, flowOpts, activationOpts);
    summaryCum = make_activation_summary_row(caseDef, metricsCum, 0);
    changeTbl = summary_change_table(summaryCum, previousSummary);
    metaCum = table(bi, targetFrames, usedFrames, instancesUsed, lastInstance, lastInstanceFrames, ...
        'VariableNames', {'FrameBin','CumulativeFrameTarget','CumulativeFramesUsed', ...
        'XMLInstancesUsed','LastXMLInstance','LastXMLFramesUsed'});

    convergenceRows = append_table_compat(convergenceRows, [metaCum, summaryCum, changeTbl]);
    previousSummary = summaryCum;
end
end

function T = summary_change_table(currentSummary, previousSummary)
varNames = currentSummary.Properties.VariableNames;
skipNames = ["Case", "Re", "kD", "elapsed_case_sec"];
changeVals = [];
changeNames = strings(0,1);

for vi = 1:numel(varNames)
    name = string(varNames{vi});
    if any(strcmp(name, skipNames))
        continue;
    end

    curRaw = currentSummary.(varNames{vi});
    if ~(isnumeric(curRaw) || islogical(curRaw))
        continue;
    end
    curVal = double(curRaw(1));
    prevVal = NaN;
    if ~isempty(previousSummary) && height(previousSummary) > 0 && ismember(varNames{vi}, previousSummary.Properties.VariableNames)
        prevRaw = previousSummary.(varNames{vi});
        if isnumeric(prevRaw) || islogical(prevRaw)
            prevVal = double(prevRaw(1));
        end
    end

    deltaVal = curVal - prevVal;
    pctVal = NaN;
    if isfinite(prevVal) && abs(prevVal) > eps
        pctVal = 100 * deltaVal / abs(prevVal);
    end

    changeVals = [changeVals, deltaVal, pctVal]; %#ok<AGROW>
    changeNames(end+1,1) = "delta_" + name; %#ok<AGROW>
    changeNames(end+1,1) = "pctChange_" + name; %#ok<AGROW>
end

if isempty(changeVals)
    T = table();
else
    T = array2table(changeVals, 'VariableNames', cellstr(changeNames));
end
end

function [outCum, usedFrames, instancesUsed, lastInstance, lastInstanceFrames] = ...
    cumulative_output_by_source_frames(outChunks, xmlFiles, targetFrames)

outCum = [];
usedFrames = 0;
instancesUsed = 0;
lastInstance = 0;
lastInstanceFrames = 0;
remainingFrames = max(0, round(targetFrames));
parts = {};
partXmlFiles = {};

for c = 1:numel(outChunks)
    if remainingFrames <= 0
        break;
    end

    outOne = outChunks{c};
    [frameMin, ~, frameCount] = parsed_output_frame_span(outOne);
    if frameCount <= 0 || ~isfinite(frameMin)
        continue;
    end

    framesThisInstance = min(remainingFrames, frameCount);
    frameStop = frameMin + framesThisInstance - 1;
    outPart = subset_parsed_output_by_frame_window(outOne, frameMin, frameStop);
    if parsed_output_has_content(outPart)
        parts{end+1,1} = outPart; %#ok<AGROW>
        partXmlFiles{end+1,1} = xmlFiles{c}; %#ok<AGROW>
    end

    usedFrames = usedFrames + framesThisInstance;
    instancesUsed = instancesUsed + 1;
    lastInstance = c;
    lastInstanceFrames = framesThisInstance;
    remainingFrames = remainingFrames - framesThisInstance;
end

if ~isempty(parts)
    outCum = merge_parsed_outputs(parts, partXmlFiles);
end
end

function tf = parsed_output_has_content(out)
tf = isstruct(out) && isfield(out, 'trajectories') && ~isempty(out.trajectories) && ...
    isfield(out, 'spots') && istable(out.spots) && height(out.spots) > 0;
end

function totalFrames = estimate_total_source_frames(outChunks)
totalFrames = 0;
for c = 1:numel(outChunks)
    [~, ~, frameCount] = parsed_output_frame_span(outChunks{c});
    if isfinite(frameCount) && frameCount > 0
        totalFrames = totalFrames + frameCount;
    end
end
end

function [frameMin, frameMax, frameCount] = parsed_output_frame_span(out)
frameMin = NaN;
frameMax = NaN;
frameCount = 0;
frames = nan(0,1);

if isstruct(out) && isfield(out, 'spots') && istable(out.spots) && ...
        height(out.spots) > 0 && ismember('FRAME', out.spots.Properties.VariableNames)
    frames = [frames; double(out.spots.FRAME(:))]; %#ok<AGROW>
end
if isstruct(out) && isfield(out, 'trajectories') && ~isempty(out.trajectories) && ...
        isstruct(out.trajectories) && isfield(out.trajectories, 'frame')
    for k = 1:numel(out.trajectories)
        frames = [frames; double(out.trajectories(k).frame(:))]; %#ok<AGROW>
    end
end

frames = frames(isfinite(frames));
if isstruct(out) && isfield(out, 'meta') && isstruct(out.meta) && ...
        isfield(out.meta, 'imageNFrames') && isfinite(out.meta.imageNFrames) && out.meta.imageNFrames > 0
    nImageFrames = round(out.meta.imageNFrames);
    if isempty(frames)
        frameMin = 0;
    else
        observedMin = min(frames);
        if observedMin > nImageFrames
            frameMin = round(observedMin);
        else
            frameMin = 0;
        end
    end
    frameMax = frameMin + nImageFrames - 1;
    frameCount = nImageFrames;
    return;
end

if isempty(frames)
    return;
end

frameMin = min(frames);
frameMax = max(frames);
frameCount = max(0, round(frameMax) - round(frameMin) + 1);
end

function outSub = subset_parsed_output_by_frame_window(out, frameStart, frameStop)
outSub = out;
if isempty(out) || ~isstruct(out)
    return;
end

keptSpotIds = nan(0,1);
if isfield(out, 'spots') && istable(out.spots) && height(out.spots) > 0 && ...
        ismember('FRAME', out.spots.Properties.VariableNames)
    spotFrame = double(out.spots.FRAME(:));
    keepSpots = isfinite(spotFrame) & spotFrame >= frameStart & spotFrame <= frameStop;
    outSub.spots = out.spots(keepSpots, :);
    if ismember('ID', outSub.spots.Properties.VariableNames)
        keptSpotIds = double(outSub.spots.ID(:));
    end
end

if isfield(out, 'trajectories')
    trajSub = out.trajectories([]);
else
    trajSub = struct([]);
end
if isfield(out, 'trajectories') && ~isempty(out.trajectories)
    for k = 1:numel(out.trajectories)
        tr = out.trajectories(k);
        if ~isfield(tr, 'frame') || isempty(tr.frame)
            continue;
        end
        trFrame = double(tr.frame(:));
        keepPoints = isfinite(trFrame) & trFrame >= frameStart & trFrame <= frameStop;
        if ~any(keepPoints)
            continue;
        end
        tr = subset_trajectory_points(tr, keepPoints);
        trajSub(end+1,1) = tr; %#ok<AGROW>
    end
end
outSub.trajectories = trajSub;

keptTrackIds = collect_traj_scalar_field_local(trajSub, 'TRACK_ID');
if isfield(out, 'tracks') && istable(out.tracks) && height(out.tracks) > 0 && ...
        ismember('TRACK_ID', out.tracks.Properties.VariableNames)
    keepTracks = ismember(out.tracks.TRACK_ID, keptTrackIds);
    outSub.tracks = out.tracks(keepTracks, :);
else
    if isfield(out, 'tracks')
        outSub.tracks = out.tracks;
    end
end

if isfield(out, 'edges') && istable(out.edges) && height(out.edges) > 0
    keepEdges = true(height(out.edges), 1);
    if ismember('TRACK_ID', out.edges.Properties.VariableNames)
        keepEdges = keepEdges & ismember(out.edges.TRACK_ID, keptTrackIds);
    end
    if ~isempty(keptSpotIds)
        if ismember('SPOT_SOURCE_ID', out.edges.Properties.VariableNames)
            keepEdges = keepEdges & ismember(out.edges.SPOT_SOURCE_ID, keptSpotIds);
        end
        if ismember('SPOT_TARGET_ID', out.edges.Properties.VariableNames)
            keepEdges = keepEdges & ismember(out.edges.SPOT_TARGET_ID, keptSpotIds);
        end
    end
    outSub.edges = out.edges(keepEdges, :);
else
    if isfield(out, 'edges')
        outSub.edges = out.edges;
    end
end

if isfield(outSub, 'meta') && isstruct(outSub.meta)
    outSub.meta.convergenceFrameStart = frameStart;
    outSub.meta.convergenceFrameStop = frameStop;
    if isfield(outSub, 'spots') && istable(outSub.spots)
        outSub.meta.nSpotsParsed = height(outSub.spots);
    end
    if isfield(outSub, 'tracks') && istable(outSub.tracks)
        outSub.meta.nTracksParsed = height(outSub.tracks);
    end
    if isfield(outSub, 'edges') && istable(outSub.edges)
        outSub.meta.nEdgesParsed = height(outSub.edges);
    end
end
end

function tr = subset_trajectory_points(tr, keepPoints)
fields = fieldnames(tr);
nPoints = numel(keepPoints);
segmentKeep = keepPoints(1:end-1) & keepPoints(2:end);

for fi = 1:numel(fields)
    name = fields{fi};
    val = tr.(name);
    if ~(isnumeric(val) || islogical(val))
        continue;
    end
    if isvector(val) && numel(val) == nPoints
        tr.(name) = val(keepPoints);
    elseif isvector(val) && numel(val) == max(nPoints - 1, 0)
        tr.(name) = val(segmentKeep);
    end
end
end

function vals = collect_traj_scalar_field_local(traj, fieldName)
vals = nan(0,1);
if isempty(traj) || ~isstruct(traj) || ~isfield(traj, fieldName)
    return;
end
vals = nan(numel(traj), 1);
for k = 1:numel(traj)
    raw = traj(k).(fieldName);
    if isempty(raw)
        continue;
    end
    vals(k) = raw(1);
end
vals = vals(isfinite(vals));
end

function outTbl = append_table_compat(outTbl, newRows)
if isempty(newRows) || height(newRows) == 0
    return;
end
if isempty(outTbl) || width(outTbl) == 0
    outTbl = newRows;
else
    outTbl = [outTbl; newRows]; %#ok<AGROW>
end
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

function xy = choose_left_moving_activation_xy(metrics)
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

function [xmlFiles, outChunks, cacheDB, cacheUpdated] = load_case_chunk_outputs( ...
    caseDef, maxTracksToParse, parserOpts, useMatCache, forceReparse, ...
    cacheDB, cacheUpdated, cachePolicyTag, verboseFlag)
xmlFiles = detect_chunk_xml_files(caseDef.xmlFile);
nXml = numel(xmlFiles);
outChunks = cell(nXml, 1);

for c = 1:nXml
    chunkCaseDef = caseDef;
    chunkCaseDef.xmlFile = xmlFiles{c};

    caseKey = build_case_cache_key(chunkCaseDef, maxTracksToParse, cachePolicyTag);
    idxCache = find(cacheDB.key == caseKey, 1, 'first');
    useCacheEntry = false;

    if useMatCache && ~forceReparse && ~isempty(idxCache)
        cachedOut = cacheDB.out{idxCache};
        if is_cache_entry_compatible(cachedOut, parserOpts)
            outChunks{c} = cachedOut;
            useCacheEntry = true;
            fprintf("  Using cache entry %d/%d for XML sample %d/%d\n", ...
                idxCache, numel(cacheDB.key), c, nXml);
        else
            fprintf("  Cache entry %d/%d incompatible for XML sample %d/%d. Reparsing.\n", ...
                idxCache, numel(cacheDB.key), c, nXml);
        end
    end

    if ~useCacheEntry
        fprintf("  Parsing XML sample %d/%d: %s\n", c, nXml, xmlFiles{c});
        outChunks{c} = analyze_trackmate_xml_arc( ...
            xmlFiles{c}, ...
            pixelSize = caseDef.pixelSize, ...
            dt        = caseDef.dt, ...
            maxTracks = maxTracksToParse, ...
            parseROI  = false, ...
            verbose   = verboseFlag, ...
            makePlots = false, ...
            parseTrackedSpotsOnly = parserOpts.parseTrackedSpotsOnly, ...
            parseFilteredTracksOnly = parserOpts.parseFilteredTracksOnly);

        if useMatCache
            if isempty(idxCache)
                cacheDB.key(end+1,1) = caseKey;
                cacheDB.out{end+1,1} = outChunks{c};
            else
                cacheDB.out{idxCache,1} = outChunks{c};
            end
            cacheUpdated = true;
        end
    end
end
end

function eventsOut = combine_breakup_event_chunks(eventsByChunk, chunkInfo)
if isempty(eventsByChunk)
    eventsOut = struct([]);
    return;
end

eventsOut = eventsByChunk{1}([]);
for c = 1:numel(eventsByChunk)
    ev = eventsByChunk{c};
    if isempty(ev)
        continue;
    end

    frameOffset = 0;
    if nargin >= 2 && numel(chunkInfo) >= c && isstruct(chunkInfo(c)) && ...
            isfield(chunkInfo(c), 'frameOffset') && isfinite(chunkInfo(c).frameOffset)
        frameOffset = chunkInfo(c).frameOffset;
    end

    ev = shift_breakup_event_frames(ev, frameOffset);
    eventsOut = [eventsOut; ev(:)]; %#ok<AGROW>
end
end

function events = shift_breakup_event_frames(events, frameOffset)
if isempty(events) || ~isfinite(frameOffset) || frameOffset == 0
    return;
end

for k = 1:numel(events)
    if isfield(events, 'parentFrame') && isfinite(events(k).parentFrame)
        events(k).parentFrame = events(k).parentFrame + frameOffset;
    end
    if isfield(events, 'childFrame') && isfinite(events(k).childFrame)
        events(k).childFrame = events(k).childFrame + frameOffset;
    end
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
