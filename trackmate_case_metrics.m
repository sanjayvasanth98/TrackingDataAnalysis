function metrics = trackmate_case_metrics(out, varargin)
%TRACKMATE_CASE_METRICS  Unified upstream/activation metrics for TrackMate trajectories.
% New call style:
%   metrics = trackmate_case_metrics(out, qcOpts, flowOpts, activationOpts)
%
% Backward-compatible call style:
%   metrics = trackmate_case_metrics(out, minNetDx, minTrackSpots, areaJumpFactor)

[qcOpts, flowOpts, activationOpts] = resolve_metric_options(varargin{:});

traj = out.trajectories;
spots = out.spots;
tracks = table();
if isfield(out, 'tracks') && istable(out.tracks)
    tracks = out.tracks;
end

nTotal = numel(traj);

% Net-left framewise set (retained as legacy outputs).
netLeftInception_xy = zeros(0,2);
netLeftActivationEvent_xy = zeros(0,2);
netLeftTrack_xy = cell(0,1);
netLeftUpstreamSize_eqd = nan(0,1);
netLeftTau_values = nan(0,1);
netLeftActivationEventFrames = nan(0,1);
netLeftActivationEventTrackIds = nan(0,1);
netLeftFrameCells = cell(0,1);
netLeftTrackIds = nan(0,1);

% Strict recirculation set (authoritative primary outputs).
strictInjInception_xy = zeros(0,2);
strictActLocation_xy = zeros(0,2);
strictTrack_xy = cell(0,1);
strictUpstreamSize_eqd = nan(0,1);
strictTau_values = nan(0,1);
strictActivationEventFrames = nan(0,1);
strictActivationEventTrackIds = nan(0,1);
strictFrameCells = cell(0,1);
strictPrimaryTrackIds = nan(0,1);
strictFrameAxis = nan(0,1);
strictFrameVisible = nan(0,1);
strictFrameActEvents = nan(0,1);
strictFrameCumExposure = nan(0,1);
strictFrameCumActEvents = nan(0,1);
netLeftFrameAxis = nan(0,1);
netLeftFrameVisible = nan(0,1);
netLeftFrameActEvents = nan(0,1);
netLeftFrameCumExposure = nan(0,1);
netLeftFrameCumActEvents = nan(0,1);

gateStats = struct();
gateStats.nTracksTotal = nTotal;
gateStats.nRejectedTooShort = 0;
gateStats.nRejectedNonFinite = 0;
gateStats.nRejectedNonMonotonicTime = 0;
gateStats.nRejectedTopology = 0;
gateStats.nRejectedFlow = 0;
gateStats.nRejectedOrigin = 0;
gateStats.nRejectedNoActivation = 0;
gateStats.nRejectedWallBand = 0;
gateStats.nInjected = 0;
gateStats.nActivated = 0;
gateStats.originThreshold = NaN;
gateStats.xStartMin = NaN;
gateStats.xStartMax = NaN;

trackCatalog = repmat(make_catalog_template(), nTotal, 1);

if isempty(traj) || isempty(spots) || ~ismember('ID', spots.Properties.VariableNames)
    metrics = pack_metrics();
    return;
end

% Map SpotID -> row.
[spotIdSorted, sortIdx] = sort(spots.ID);
spotRowById = containers.Map('KeyType', 'double', 'ValueType', 'double');
for k = 1:numel(spotIdSorted)
    spotRowById(spotIdSorted(k)) = sortIdx(k);
end

% Map TRACK_ID -> row in tracks table (for strict legacy topology QC).
trackRowById = containers.Map('KeyType', 'double', 'ValueType', 'double');
if ~isempty(tracks) && ismember('TRACK_ID', tracks.Properties.VariableNames)
    for k = 1:height(tracks)
        tid = tracks.TRACK_ID(k);
        if isfinite(tid) && ~isKey(trackRowById, tid)
            trackRowById(tid) = k;
        end
    end
end

% Precompute strict-origin threshold once.
xStarts = nan(0,1);
for k = 1:nTotal
    xk = get_struct_column(traj(k), 'x_phys');
    if ~isempty(xk) && isfinite(xk(1))
        xStarts(end+1,1) = xk(1); %#ok<AGROW>
    end
end
[originThreshold, xMinAll, xMaxAll] = compute_origin_threshold(xStarts, flowOpts);
gateStats.originThreshold = originThreshold;
gateStats.xStartMin = xMinAll;
gateStats.xStartMax = xMaxAll;

% Build per-track prepared records once (single source of truth).
prep = repmat(make_prep_template(), nTotal, 1);
for k = 1:nTotal
    prep(k) = build_track_prep(traj(k), spotRowById, spots, qcOpts, activationOpts);
    trackCatalog(k) = build_track_catalog_entry(prep(k));

    % Net-left framewise set (legacy comparison).
    if ~(prep(k).isBasicValid && prep(k).isLeftMoving)
        continue;
    end
    if isfinite(qcOpts.maxLeftMovingTracks) && qcOpts.maxLeftMovingTracks >= 0 ...
            && numel(netLeftTrackIds) >= qcOpts.maxLeftMovingTracks
        continue;
    end

    netLeftTrackIds(end+1,1) = prep(k).TRACK_ID; %#ok<AGROW>
    netLeftTrack_xy{end+1,1} = [prep(k).x, prep(k).y]; %#ok<AGROW>
    netLeftInception_xy(end+1,:) = [prep(k).x(1), prep(k).y(1)]; %#ok<AGROW>

    frameUnique = unique(prep(k).frame(isfinite(prep(k).frame)));
    if isempty(frameUnique)
        frameUnique = (1:numel(prep(k).x)).';
    end
    netLeftFrameCells{end+1,1} = frameUnique(:); %#ok<AGROW>

    nUse = numel(prep(k).areaVals);
    if prep(k).hasActivation
        nUse = prep(k).idxJump;
    end

    if nUse >= 1
        areaUse = prep(k).areaVals(1:nUse);
        eqd = sqrt(4 .* areaUse ./ pi);
        eqd = eqd(isfinite(eqd) & eqd > 0);
        if ~isempty(eqd)
            netLeftUpstreamSize_eqd = [netLeftUpstreamSize_eqd; eqd(:)]; %#ok<AGROW>
        end
    end

    if prep(k).hasActivation
        netLeftActivationEvent_xy(end+1,:) = [prep(k).actX, prep(k).actY]; %#ok<AGROW>
        netLeftActivationEventFrames(end+1,1) = prep(k).actFrame; %#ok<AGROW>
        netLeftActivationEventTrackIds(end+1,1) = prep(k).TRACK_ID; %#ok<AGROW>
        netLeftTau_values(end+1,1) = max(prep(k).t(prep(k).actIdx) - prep(k).t(1), 0); %#ok<AGROW>
    end
end

% Strict gate accounting + authoritative strict primary set.
for k = 1:nTotal
    if prep(k).isShort
        gateStats.nRejectedTooShort = gateStats.nRejectedTooShort + 1;
        continue;
    end
    if prep(k).isNonFinite
        gateStats.nRejectedNonFinite = gateStats.nRejectedNonFinite + 1;
        continue;
    end
    if prep(k).isNonMonotonicTime
        gateStats.nRejectedNonMonotonicTime = gateStats.nRejectedNonMonotonicTime + 1;
        continue;
    end

    if ~passes_topology_gate(traj(k), trackRowById, tracks, qcOpts)
        gateStats.nRejectedTopology = gateStats.nRejectedTopology + 1;
        continue;
    end

    [isCounterflow, ~, ~] = passes_counterflow_gate(prep(k).x, flowOpts);
    if ~isCounterflow
        gateStats.nRejectedFlow = gateStats.nRejectedFlow + 1;
        continue;
    end

    if flowOpts.requireRightOrigin && isfinite(originThreshold)
        if ~passes_origin_gate(prep(k).x(1), originThreshold, flowOpts)
            gateStats.nRejectedOrigin = gateStats.nRejectedOrigin + 1;
            continue;
        end
    end

    gateStats.nInjected = gateStats.nInjected + 1;
    trackCatalog(k).isStrictPrimary = true;

    strictPrimaryTrackIds(end+1,1) = prep(k).TRACK_ID; %#ok<AGROW>
    strictTrack_xy{end+1,1} = [prep(k).x, prep(k).y]; %#ok<AGROW>
    strictInjInception_xy(end+1,:) = [prep(k).x(1), prep(k).y(1)]; %#ok<AGROW>

    frameUnique = unique(prep(k).frame(isfinite(prep(k).frame)));
    if isempty(frameUnique)
        frameUnique = (1:numel(prep(k).x)).';
    end
    strictFrameCells{end+1,1} = frameUnique(:); %#ok<AGROW>

    nUse = numel(prep(k).areaVals);
    if prep(k).hasActivation
        nUse = prep(k).idxJump;
    end
    if nUse >= 1
        areaUse = prep(k).areaVals(1:nUse);
        eqd = sqrt(4 .* areaUse ./ pi);
        eqd = eqd(isfinite(eqd) & eqd > 0);
        if ~isempty(eqd)
            strictUpstreamSize_eqd = [strictUpstreamSize_eqd; eqd(:)]; %#ok<AGROW>
        end
    end

    if ~prep(k).hasActivation
        gateStats.nRejectedNoActivation = gateStats.nRejectedNoActivation + 1;
        continue;
    end

    if qcOpts.wallBandEnabled && ~passes_wall_band(prep(k).actY, qcOpts.wallBandYLimits_mm)
        gateStats.nRejectedWallBand = gateStats.nRejectedWallBand + 1;
        continue;
    end

    gateStats.nActivated = gateStats.nActivated + 1;
    trackCatalog(k).isStrictActivated = true;
    trackCatalog(k).strictActivationFrame = prep(k).actFrame;
    trackCatalog(k).strictActivationX = prep(k).actX;
    trackCatalog(k).strictActivationY = prep(k).actY;
    trackCatalog(k).strictActivationIndex = prep(k).actIdx;
    strictActLocation_xy(end+1,:) = [prep(k).actX, prep(k).actY]; %#ok<AGROW>
    strictActivationEventFrames(end+1,1) = prep(k).actFrame; %#ok<AGROW>
    strictActivationEventTrackIds(end+1,1) = prep(k).TRACK_ID; %#ok<AGROW>
    strictTau_values(end+1,1) = max(prep(k).t(prep(k).actIdx) - prep(k).t(1), 0); %#ok<AGROW>
end

% Framewise counts.
[strictFrameAxis, strictFrameVisible, strictFrameActEvents, strictFrameCumExposure, strictFrameCumActEvents] = ...
    build_framewise_counts(strictFrameCells, strictActivationEventFrames);
[netLeftFrameAxis, netLeftFrameVisible, netLeftFrameActEvents, netLeftFrameCumExposure, netLeftFrameCumActEvents] = ...
    build_framewise_counts(netLeftFrameCells, netLeftActivationEventFrames);

metrics = pack_metrics();

    function outMetrics = pack_metrics()
        nStrictPrimaryTracks = numel(strictPrimaryTrackIds);
        strictActTrackIdsFinite = strictActivationEventTrackIds(isfinite(strictActivationEventTrackIds));
        nStrictActivatedTracks = numel(unique(strictActTrackIdsFinite));
        if nStrictActivatedTracks == 0
            nStrictActivatedTracks = size(strictActLocation_xy, 1);
        end

        strictTrackFrameExposure = sum(strictFrameVisible);
        strictActivationEventsTotal = sum(strictFrameActEvents);
        A_over_I = strictActivationEventsTotal / max(strictTrackFrameExposure, 1);
        [A_over_I_ci_low, A_over_I_ci_high] = wilson_ci(strictActivationEventsTotal, strictTrackFrameExposure, 0.95);

        nNetLeftTracks = numel(netLeftTrackIds);
        netLeftActTrackIdsFinite = netLeftActivationEventTrackIds(isfinite(netLeftActivationEventTrackIds));
        nNetLeftActivatedTracks = numel(unique(netLeftActTrackIdsFinite));
        if nNetLeftActivatedTracks == 0
            nNetLeftActivatedTracks = size(netLeftActivationEvent_xy, 1);
        end
        netLeftTrackFrameExposure = sum(netLeftFrameVisible);
        netLeftActivationEventsTotal = sum(netLeftFrameActEvents);
        A_over_I_netLeftLegacy = netLeftActivationEventsTotal / max(netLeftTrackFrameExposure, 1);
        [A_over_I_ci_low_netLeftLegacy, A_over_I_ci_high_netLeftLegacy] = ...
            wilson_ci(netLeftActivationEventsTotal, netLeftTrackFrameExposure, 0.95);

        nInjectedStrict = gateStats.nInjected;
        nActivatedStrict = gateStats.nActivated;
        A_over_I_strictLegacy = nActivatedStrict / max(nInjectedStrict, 1);
        [A_over_I_ci_low_strictLegacy, A_over_I_ci_high_strictLegacy] = ...
            wilson_ci(nActivatedStrict, nInjectedStrict, 0.95);

        if isempty(strictFrameVisible)
            strictMeanVisiblePerFrame = 0;
            strictPeakVisiblePerFrame = 0;
        else
            strictMeanVisiblePerFrame = mean(strictFrameVisible);
            strictPeakVisiblePerFrame = max(strictFrameVisible);
        end

        if isempty(netLeftFrameVisible)
            netLeftMeanVisiblePerFrame = 0;
            netLeftPeakVisiblePerFrame = 0;
        else
            netLeftMeanVisiblePerFrame = mean(netLeftFrameVisible);
            netLeftPeakVisiblePerFrame = max(netLeftFrameVisible);
        end

        nBasicValidTracks = sum([prep.isBasicValid]);

        outMetrics = struct();
        outMetrics.nTracksTotal = nTotal;
        outMetrics.nBasicValidTracks = nBasicValidTracks;

        % Strict-primary metrics (authoritative framewise workflow).
        outMetrics.nStrictPrimaryTracks = nStrictPrimaryTracks;
        outMetrics.nStrictActivatedTracks = nStrictActivatedTracks;
        outMetrics.strictTrackFrameExposure = strictTrackFrameExposure;
        outMetrics.strictActivationEventsTotal = strictActivationEventsTotal;
        outMetrics.strictMeanVisiblePerFrame = strictMeanVisiblePerFrame;
        outMetrics.strictPeakVisiblePerFrame = strictPeakVisiblePerFrame;
        outMetrics.A_over_I = A_over_I;
        outMetrics.A_over_I_ci_low = A_over_I_ci_low;
        outMetrics.A_over_I_ci_high = A_over_I_ci_high;

        % Strict-primary authoritative arrays.
        outMetrics.strictPrimaryTrackIds = strictPrimaryTrackIds;
        outMetrics.strictActivatedTrackIds = unique(strictActTrackIdsFinite);
        outMetrics.strictTrack_xy = strictTrack_xy;
        outMetrics.strictInception_xy = strictInjInception_xy;
        outMetrics.strictActivationEvent_xy = strictActLocation_xy;
        outMetrics.strictActivationEvent_frame = strictActivationEventFrames;
        outMetrics.strictActivationEvent_trackId = strictActivationEventTrackIds;
        outMetrics.strict_tau_values = strictTau_values;
        outMetrics.strict_upstreamSize_eqd = strictUpstreamSize_eqd;
        outMetrics.strict_frame_axis = strictFrameAxis;
        outMetrics.strict_frame_nVisible = strictFrameVisible;
        outMetrics.strict_frame_nActivationEvents = strictFrameActEvents;
        outMetrics.strict_frame_cumExposure = strictFrameCumExposure;
        outMetrics.strict_frame_cumActivationEvents = strictFrameCumActEvents;

        % Compatibility aliases (primary now strict recirculation).
        outMetrics.nLeftMovingTracks = nStrictPrimaryTracks;
        outMetrics.nActivatedLeftMovingTracks = nStrictActivatedTracks;
        outMetrics.leftMovingTrackFrameExposure = strictTrackFrameExposure;
        outMetrics.activationEventsTotal = strictActivationEventsTotal;
        outMetrics.meanLeftMovingPerFrame = strictMeanVisiblePerFrame;
        outMetrics.peakLeftMovingPerFrame = strictPeakVisiblePerFrame;
        outMetrics.nInjected = nStrictPrimaryTracks;
        outMetrics.nActivated = nStrictActivatedTracks;
        outMetrics.injInception_xy = strictInjInception_xy;
        outMetrics.actLocation_xy = strictActLocation_xy;
        outMetrics.inception2x_xy = strictActLocation_xy;
        outMetrics.upstreamTrack_xy = strictTrack_xy;
        outMetrics.leftMovingTrack_xy = strictTrack_xy;
        outMetrics.activationEvent_xy = strictActLocation_xy;
        outMetrics.activationEvent_frame = strictActivationEventFrames;
        outMetrics.activationEvent_trackId = strictActivationEventTrackIds;
        outMetrics.tau_values = strictTau_values;
        outMetrics.tau_mean = mean(strictTau_values, 'omitnan');
        outMetrics.tau_std = std(strictTau_values, 0, 'omitnan');
        outMetrics.upstreamSize_eqd = strictUpstreamSize_eqd;
        outMetrics.frame_axis = strictFrameAxis;
        outMetrics.frame_nLeftMovingVisible = strictFrameVisible;
        outMetrics.frame_nActivationEvents = strictFrameActEvents;
        outMetrics.frame_cumExposure = strictFrameCumExposure;
        outMetrics.frame_cumActivationEvents = strictFrameCumActEvents;

        % Net-left legacy framewise outputs for audit continuity.
        outMetrics.nLeftMovingTracks_netLeftLegacy = nNetLeftTracks;
        outMetrics.nActivatedLeftMovingTracks_netLeftLegacy = nNetLeftActivatedTracks;
        outMetrics.leftMovingTrackFrameExposure_netLeftLegacy = netLeftTrackFrameExposure;
        outMetrics.activationEventsTotal_netLeftLegacy = netLeftActivationEventsTotal;
        outMetrics.meanLeftMovingPerFrame_netLeftLegacy = netLeftMeanVisiblePerFrame;
        outMetrics.peakLeftMovingPerFrame_netLeftLegacy = netLeftPeakVisiblePerFrame;
        outMetrics.leftMovingTrackIds_netLeftLegacy = unique(netLeftTrackIds(isfinite(netLeftTrackIds)), 'stable');
        outMetrics.maxLeftMovingTracks_netLeftLegacy = qcOpts.maxLeftMovingTracks;
        outMetrics.A_over_I_netLeftLegacy = A_over_I_netLeftLegacy;
        outMetrics.A_over_I_ci_low_netLeftLegacy = A_over_I_ci_low_netLeftLegacy;
        outMetrics.A_over_I_ci_high_netLeftLegacy = A_over_I_ci_high_netLeftLegacy;
        outMetrics.injInception_xy_netLeftLegacy = netLeftInception_xy;
        outMetrics.activationEvent_xy_netLeftLegacy = netLeftActivationEvent_xy;
        outMetrics.activationEvent_frame_netLeftLegacy = netLeftActivationEventFrames;
        outMetrics.activationEvent_trackId_netLeftLegacy = netLeftActivationEventTrackIds;
        outMetrics.leftMovingTrack_xy_netLeftLegacy = netLeftTrack_xy;
        outMetrics.tau_values_netLeftLegacy = netLeftTau_values;
        outMetrics.upstreamSize_eqd_netLeftLegacy = netLeftUpstreamSize_eqd;
        outMetrics.frame_axis_netLeftLegacy = netLeftFrameAxis;
        outMetrics.frame_nLeftMovingVisible_netLeftLegacy = netLeftFrameVisible;
        outMetrics.frame_nActivationEvents_netLeftLegacy = netLeftFrameActEvents;
        outMetrics.frame_cumExposure_netLeftLegacy = netLeftFrameCumExposure;
        outMetrics.frame_cumActivationEvents_netLeftLegacy = netLeftFrameCumActEvents;

        % Secondary strict legacy metrics.
        outMetrics.nInjected_strictLegacy = nInjectedStrict;
        outMetrics.nActivated_strictLegacy = nActivatedStrict;
        outMetrics.A_over_I_strictLegacy = A_over_I_strictLegacy;
        outMetrics.A_over_I_ci_low_strictLegacy = A_over_I_ci_low_strictLegacy;
        outMetrics.A_over_I_ci_high_strictLegacy = A_over_I_ci_high_strictLegacy;
        outMetrics.injInception_xy_strictLegacy = strictInjInception_xy;
        outMetrics.actLocation_xy_strictLegacy = strictActLocation_xy;

        outMetrics.trackCatalog = trackCatalog;
        outMetrics.gateStats = gateStats;
        outMetrics.qcOpts = qcOpts;
        outMetrics.flowOpts = flowOpts;
        outMetrics.activationOpts = activationOpts;
    end
end

function tpl = make_prep_template()
tpl = struct( ...
    'TRACK_ID', NaN, ...
    'spotIds', nan(0,1), ...
    'nTrackSpots', 0, ...
    'x', nan(0,1), ...
    'y', nan(0,1), ...
    't', nan(0,1), ...
    'frame', nan(0,1), ...
    'areaVals', nan(0,1), ...
    'isShort', false, ...
    'isNonFinite', false, ...
    'isNonMonotonicTime', false, ...
    'isBasicValid', false, ...
    'isLeftMoving', false, ...
    'idxJump', NaN, ...
    'hasActivation', false, ...
    'actIdx', NaN, ...
    'actX', NaN, ...
    'actY', NaN, ...
    'actFrame', NaN);
end

function tpl = make_catalog_template()
tpl = struct( ...
    'TRACK_ID', NaN, ...
    'frame', nan(0,1), ...
    'x', nan(0,1), ...
    'y', nan(0,1), ...
    'isBasicValid', false, ...
    'isLeftMoving', false, ...
    'isStrictPrimary', false, ...
    'isStrictActivated', false, ...
    'isActivated', false, ...
    'activationFrame', NaN, ...
    'activationX', NaN, ...
    'activationY', NaN, ...
    'activationIndex', NaN, ...
    'strictActivationFrame', NaN, ...
    'strictActivationX', NaN, ...
    'strictActivationY', NaN, ...
    'strictActivationIndex', NaN);
end

function prep = build_track_prep(thisTraj, spotRowById, spots, qcOpts, activationOpts)
prep = make_prep_template();

if isfield(thisTraj, 'TRACK_ID')
    prep.TRACK_ID = thisTraj.TRACK_ID;
end
if isfield(thisTraj, 'spotIds')
    prep.spotIds = thisTraj.spotIds(:);
end
prep.nTrackSpots = numel(prep.spotIds);
prep.isShort = (prep.nTrackSpots < qcOpts.minTrackSpots);

if isfield(thisTraj, 'x_phys')
    prep.x = thisTraj.x_phys(:);
end
if isfield(thisTraj, 'y_phys')
    prep.y = thisTraj.y_phys(:);
end
if isfield(thisTraj, 't')
    prep.t = thisTraj.t(:);
end

n = numel(prep.x);
rawFrame = nan(0,1);
if isfield(thisTraj, 'frame')
    rawFrame = thisTraj.frame(:);
end
prep.frame = normalize_frame_values(rawFrame, n);

if n < 2 || numel(prep.y) ~= n || numel(prep.t) ~= n || any(~isfinite([prep.x; prep.y; prep.t]))
    prep.isNonFinite = true;
    return;
end

dtTrack = diff(prep.t);
if isempty(dtTrack) || any(~isfinite(dtTrack)) || any(dtTrack <= 0)
    prep.isNonMonotonicTime = true;
    return;
end

if prep.isShort
    return;
end

prep.isBasicValid = true;
prep.isLeftMoving = (prep.x(end) < prep.x(1));

prep.areaVals = extract_area_values(prep.spotIds, spotRowById, spots);
idxJump = find_sustained_growth_activation(prep.areaVals, activationOpts);
if ~isempty(idxJump)
    actIdx = idxJump + 1;
    if actIdx >= 1 && actIdx <= n
        prep.idxJump = idxJump;
        prep.hasActivation = true;
        prep.actIdx = actIdx;
        prep.actX = prep.x(actIdx);
        prep.actY = prep.y(actIdx);
        prep.actFrame = prep.frame(actIdx);
    end
end
end

function catalog = build_track_catalog_entry(prep)
catalog = make_catalog_template();
catalog.TRACK_ID = prep.TRACK_ID;
catalog.frame = prep.frame;
catalog.x = prep.x;
catalog.y = prep.y;
catalog.isBasicValid = prep.isBasicValid;
catalog.isLeftMoving = prep.isLeftMoving;
catalog.isStrictPrimary = false;
catalog.isStrictActivated = false;
catalog.isActivated = prep.hasActivation;
catalog.activationFrame = prep.actFrame;
catalog.activationX = prep.actX;
catalog.activationY = prep.actY;
catalog.activationIndex = prep.actIdx;
catalog.strictActivationFrame = NaN;
catalog.strictActivationX = NaN;
catalog.strictActivationY = NaN;
catalog.strictActivationIndex = NaN;
end

function areaVals = extract_area_values(spotIds, spotRowById, spots)
areaVals = nan(numel(spotIds), 1);
if isempty(spotIds) || ~ismember('AREA', spots.Properties.VariableNames)
    return;
end

for ii = 1:numel(spotIds)
    sid = spotIds(ii);
    if isKey(spotRowById, sid)
        row = spotRowById(sid);
        areaVals(ii) = spots.AREA(row);
    end
end
end

function frameVals = normalize_frame_values(rawFrame, n)
if nargin < 2 || ~isfinite(n) || n < 0
    n = numel(rawFrame);
end

if n == 0
    frameVals = nan(0,1);
    return;
end

if isempty(rawFrame) || numel(rawFrame) ~= n || any(~isfinite(rawFrame(:)))
    frameVals = (1:n).';
    return;
end

frameVals = rawFrame(:);
r = round(frameVals);
if max(abs(frameVals - r)) < 1e-9
    frameVals = r;
end
end

function [frameAxis, nLeftVisible, nActEvents, cumExposure, cumActEvents] = ...
    build_framewise_counts(leftMovingFrameCells, activationFrames)

allFrames = nan(0,1);
for i = 1:numel(leftMovingFrameCells)
    f = leftMovingFrameCells{i};
    if isempty(f)
        continue;
    end
    f = unique(f(isfinite(f)));
    if ~isempty(f)
        allFrames = [allFrames; f(:)]; %#ok<AGROW>
    end
end

activationFrames = activationFrames(:);
activationFrames = activationFrames(isfinite(activationFrames));

frameAxis = unique([allFrames; activationFrames]);
if isempty(frameAxis)
    nLeftVisible = nan(0,1);
    nActEvents = nan(0,1);
    cumExposure = nan(0,1);
    cumActEvents = nan(0,1);
    return;
end

nFrames = numel(frameAxis);
nLeftVisible = zeros(nFrames,1);
nActEvents = zeros(nFrames,1);

if ~isempty(allFrames)
    [~, locLeft] = ismember(allFrames, frameAxis);
    locLeft = locLeft(locLeft > 0);
    if ~isempty(locLeft)
        nLeftVisible = accumarray(locLeft, 1, [nFrames, 1]);
    end
end

if ~isempty(activationFrames)
    [~, locAct] = ismember(activationFrames, frameAxis);
    locAct = locAct(locAct > 0);
    if ~isempty(locAct)
        nActEvents = accumarray(locAct, 1, [nFrames, 1]);
    end
end

cumExposure = cumsum(nLeftVisible);
cumActEvents = cumsum(nActEvents);
end

function x = get_struct_column(s, fieldName)
x = nan(0,1);
if isfield(s, fieldName)
    x = s.(fieldName)(:);
end
end

function pass = passes_topology_gate(thisTraj, trackRowById, tracks, qcOpts)
pass = true;

if isempty(tracks) || ~isKey(trackRowById, thisTraj.TRACK_ID)
    return;
end

row = trackRowById(thisTraj.TRACK_ID);
nGaps = get_track_field(tracks, row, 'NUMBER_GAPS', 0);
nSplits = get_track_field(tracks, row, 'NUMBER_SPLITS', 0);
nMerges = get_track_field(tracks, row, 'NUMBER_MERGES', 0);
nComplex = get_track_field(tracks, row, 'NUMBER_COMPLEX', 0);

if isfinite(qcOpts.maxTrackGaps) && isfinite(nGaps) && (nGaps > qcOpts.maxTrackGaps)
    pass = false;
    return;
end

if qcOpts.rejectSplitMergeComplex
    if (nSplits > 0) || (nMerges > 0) || (nComplex > 0)
        pass = false;
    end
end
end

function v = get_track_field(tracks, row, fieldName, defaultValue)
v = defaultValue;
if isempty(tracks) || row < 1 || row > height(tracks)
    return;
end
if ~ismember(fieldName, tracks.Properties.VariableNames)
    return;
end
raw = tracks{row, fieldName};
if isempty(raw)
    return;
end
v = raw(1);
end

function [pass, negFrac, posFrac] = passes_counterflow_gate(x, flowOpts)
pass = false;
negFrac = NaN;
posFrac = NaN;

if numel(x) < 2
    return;
end

dx = diff(x);
dx = dx(isfinite(dx));
if isempty(dx)
    return;
end

negFrac = mean(dx < 0);
posFrac = mean(dx > 0);
netDx = x(end) - x(1);

if strcmpi(flowOpts.bulkDirection, 'left_to_right')
    pass = (netDx <= -flowOpts.minNetDxCounterflow_mm) && ...
        (negFrac >= flowOpts.minNegativeStepFraction) && ...
        (posFrac <= flowOpts.maxPositiveStepFraction);
else
    pass = (netDx >= flowOpts.minNetDxCounterflow_mm) && ...
        (posFrac >= flowOpts.minNegativeStepFraction) && ...
        (negFrac <= flowOpts.maxPositiveStepFraction);
end
end

function pass = passes_origin_gate(xStart, originThreshold, flowOpts)
if ~flowOpts.requireRightOrigin || ~isfinite(originThreshold) || ~isfinite(xStart)
    pass = true;
    return;
end

if strcmpi(flowOpts.bulkDirection, 'left_to_right')
    pass = (xStart >= originThreshold);
else
    pass = (xStart <= originThreshold);
end
end

function [originThreshold, xMin, xMax] = compute_origin_threshold(xStarts, flowOpts)
originThreshold = NaN;
xMin = NaN;
xMax = NaN;

if isempty(xStarts)
    return;
end

xStarts = xStarts(isfinite(xStarts));
if isempty(xStarts)
    return;
end

xMin = min(xStarts);
xMax = max(xStarts);

if ~isfinite(xMin) || ~isfinite(xMax)
    return;
end

originFrac = min(1, max(0, flowOpts.rightOriginFrac));
rangeX = xMax - xMin;

if strcmpi(flowOpts.bulkDirection, 'left_to_right')
    originThreshold = xMin + originFrac * rangeX;
else
    originThreshold = xMin + (1 - originFrac) * rangeX;
end
end

function pass = passes_wall_band(yVal, yLimits)
pass = true;
if isempty(yLimits) || numel(yLimits) ~= 2 || any(~isfinite(yLimits))
    return;
end
pass = (yVal >= min(yLimits)) && (yVal <= max(yLimits));
end

function idxJump = find_sustained_growth_activation(areaVals, opts)
idxJump = [];
n = numel(areaVals);

if n < 2
    return;
end

for i = 1:(n - 1)
    aCurr = areaVals(i);
    aNext = areaVals(i + 1);

    if ~(isfinite(aCurr) && isfinite(aNext) && aCurr > 0 && aNext > 0)
        continue;
    end

    preStart = max(1, i - opts.preWindowFrames + 1);
    preVals = areaVals(preStart:i);
    preVals = preVals(isfinite(preVals) & preVals > 0);
    if numel(preVals) < opts.minPrePoints
        continue;
    end

    aPre = median(preVals);
    if ~(isfinite(aPre) && aPre > 0)
        continue;
    end

    if (aNext / aCurr) < opts.areaJumpFactor
        continue;
    end
    if (aNext / aPre) < opts.areaJumpFactor
        continue;
    end

    postStart = i + 2;
    postEnd = min(n, i + 1 + opts.requiredPostFrames);
    postVals = nan(0,1);
    if postStart <= postEnd
        postVals = areaVals(postStart:postEnd);
    end
    postVals = postVals(isfinite(postVals) & postVals > 0);

    if numel(postVals) >= opts.requiredPostFrames
        postRatio = postVals ./ aPre;
        if median(postRatio) >= opts.postMedianFactor && max(postRatio) >= opts.postMaxFactor
            idxJump = i;
            return;
        end
    elseif numel(postVals) == 1 && opts.enableBurstFallback
        if (aNext / aPre) >= opts.burstJumpFactor && (postVals(1) / aPre) >= opts.burstPostFactor
            idxJump = i;
            return;
        end
    else
        % No post-jump evidence -> reject this jump candidate.
    end
end
end

function [qcOpts, flowOpts, activationOpts] = resolve_metric_options(varargin)
qcOpts = default_qc_opts();
flowOpts = default_flow_opts();
activationOpts = default_activation_opts();

if isempty(varargin)
    return;
end

if numel(varargin) >= 3 && isnumeric(varargin{1}) && isnumeric(varargin{2}) && isnumeric(varargin{3})
    % Backward compatibility mode.
    flowOpts.minNetDxCounterflow_mm = max(0, double(varargin{1}));
    qcOpts.minTrackSpots = max(2, round(double(varargin{2})));
    activationOpts.areaJumpFactor = max(1, double(varargin{3}));
    return;
end

if numel(varargin) >= 1 && isstruct(varargin{1}) && ~isempty(varargin{1})
    qcOpts = merge_struct(qcOpts, varargin{1});
end
if numel(varargin) >= 2 && isstruct(varargin{2}) && ~isempty(varargin{2})
    flowOpts = merge_struct(flowOpts, varargin{2});
end
if numel(varargin) >= 3 && isstruct(varargin{3}) && ~isempty(varargin{3})
    activationOpts = merge_struct(activationOpts, varargin{3});
end

qcOpts.minTrackSpots = max(2, round(qcOpts.minTrackSpots));
qcOpts.maxTrackGaps = max(0, round(qcOpts.maxTrackGaps));
if ~isfield(qcOpts, 'maxLeftMovingTracks') || isempty(qcOpts.maxLeftMovingTracks) || ~isfinite(qcOpts.maxLeftMovingTracks)
    qcOpts.maxLeftMovingTracks = Inf;
else
    qcOpts.maxLeftMovingTracks = max(0, round(double(qcOpts.maxLeftMovingTracks)));
end
flowOpts.minNetDxCounterflow_mm = max(0, flowOpts.minNetDxCounterflow_mm);
flowOpts.minNegativeStepFraction = min(1, max(0, flowOpts.minNegativeStepFraction));
flowOpts.maxPositiveStepFraction = min(1, max(0, flowOpts.maxPositiveStepFraction));
flowOpts.rightOriginFrac = min(1, max(0, flowOpts.rightOriginFrac));
activationOpts.areaJumpFactor = max(1, activationOpts.areaJumpFactor);
activationOpts.preWindowFrames = max(1, round(activationOpts.preWindowFrames));
activationOpts.minPrePoints = max(1, round(activationOpts.minPrePoints));
activationOpts.requiredPostFrames = max(1, round(activationOpts.requiredPostFrames));
end

function out = merge_struct(base, override)
out = base;
fn = fieldnames(override);
for i = 1:numel(fn)
    out.(fn{i}) = override.(fn{i});
end
end

function qcOpts = default_qc_opts()
qcOpts = struct();
qcOpts.minTrackSpots = 5;
qcOpts.maxTrackGaps = 1;
qcOpts.maxLeftMovingTracks = Inf;
qcOpts.rejectSplitMergeComplex = true;
qcOpts.wallBandEnabled = false;
qcOpts.wallBandYLimits_mm = [];
end

function flowOpts = default_flow_opts()
flowOpts = struct();
flowOpts.bulkDirection = "left_to_right";
flowOpts.minNetDxCounterflow_mm = 0.08;
flowOpts.minNegativeStepFraction = 0.65;
flowOpts.maxPositiveStepFraction = 0.30;
flowOpts.requireRightOrigin = true;
flowOpts.rightOriginFrac = 0.60;
end

function activationOpts = default_activation_opts()
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
end

function [ciLow, ciHigh] = wilson_ci(k, n, confidenceLevel)
if nargin < 3 || ~isfinite(confidenceLevel)
    confidenceLevel = 0.95;
end

if ~(isfinite(n) && n > 0 && isfinite(k) && k >= 0)
    ciLow = NaN;
    ciHigh = NaN;
    return;
end

z = -sqrt(2) * erfcinv(confidenceLevel + 1);
pHat = k / n;
z2 = z^2;
denom = 1 + z2 / n;
center = (pHat + z2 / (2 * n)) / denom;
radius = z * sqrt((pHat * (1 - pHat) / n) + (z2 / (4 * n^2))) / denom;

ciLow = max(0, center - radius);
ciHigh = min(1, center + radius);
end
