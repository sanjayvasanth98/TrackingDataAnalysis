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
netLeftMemberMask = false(nTotal,1);

% Non-strict microbubble-start set (diagnostic inclusion for activated tracks).
microRescueTrack_xy = cell(0,1);
microRescueTrackIds = nan(0,1);
microRescueActivationEvent_xy = zeros(0,2);
microRescueActivationEventFrames = nan(0,1);
microRescueActivationEventTrackIds = nan(0,1);
microRescueActivationSeedArea_px2 = nan(0,1);

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
gateStats.nRejectedExcessiveYStep = 0;
gateStats.nRejectedOriginWindow = 0;
gateStats.nRejectedTopology = 0;
gateStats.nRejectedFlow = 0;
gateStats.nRejectedOrigin = 0;
gateStats.nRejectedNoActivation = 0;
gateStats.nRejectedWallBand = 0;
gateStats.nRejectedProximityMerge = 0;
gateStats.nRejectedUnwantedArea = 0;
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
    prep(k) = build_track_prep(traj(k), spotRowById, spots, qcOpts, flowOpts, activationOpts);
end

% Post-processing: reject activations caused by proximity merge (overlapping blobs).
[prep, nProxRejected] = reject_proximity_merge_activations(prep, spots, activationOpts);
gateStats.nRejectedProximityMerge = nProxRejected;

% Build catalogs and net-left accounting.
for k = 1:nTotal
    trackCatalog(k) = build_track_catalog_entry(prep(k));

    % Net-left framewise set (legacy comparison).
    isNetLeftMember = prep(k).isBasicValid && prep(k).isLeftMoving;
    usedBandRescue = false;
    if ~isNetLeftMember
        usedBandRescue = prep(k).isBasicValid && passes_netleft_band_rescue(prep(k), flowOpts);
        isNetLeftMember = usedBandRescue;
    end
    if isNetLeftMember && flowOpts.requireRightOrigin && isfinite(originThreshold)
        if ~(usedBandRescue && flowOpts.netLeftBandRescueBypassOriginGate)
            isNetLeftMember = passes_origin_gate(prep(k).x(1), originThreshold, flowOpts);
        end
    end
    netLeftMemberMask(k) = isNetLeftMember;
    if ~isNetLeftMember
        continue;
    end
    if prep(k).isUnwantedArea
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
        eqd = sqrt(4 .* areaUse ./ pi) * activationOpts.pixelSize; % px -> mm
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
    if prep(k).isOriginExcludedBox
        gateStats.nRejectedOriginWindow = gateStats.nRejectedOriginWindow + 1;
        continue;
    end
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
    if prep(k).isExcessiveYStep
        gateStats.nRejectedExcessiveYStep = gateStats.nRejectedExcessiveYStep + 1;
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

    if prep(k).isUnwantedArea
        gateStats.nRejectedUnwantedArea = gateStats.nRejectedUnwantedArea + 1;
        continue;
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
        eqd = sqrt(4 .* areaUse ./ pi) * activationOpts.pixelSize; % px -> mm
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

% Include microbubble-start tracks (configured px^2 range) outside strict set.
if activationOpts.includeMicrobubbleActivationRescue
    startRange = sort(double(activationOpts.microbubbleStartAreaRange_px2(1:2)));
    for k = 1:nTotal
        if ~prep(k).isBasicValid
            continue;
        end
        if activationOpts.microbubbleRequireOutsideStrictPrimary && is_true_field(trackCatalog(k), 'isStrictPrimary')
            continue;
        end
        if prep(k).isUnwantedArea
            continue;
        end

        startArea = prep(k).startArea;
        if ~(isfinite(startArea) && startArea >= startRange(1) && startArea <= startRange(2))
            continue;
        end

        trackCatalog(k).isMicrobubbleRescueNonLeft = true;
        trackCatalog(k).microbubbleSeedArea_px2 = startArea;

        microRescueTrackIds(end+1,1) = prep(k).TRACK_ID; %#ok<AGROW>
        microRescueTrack_xy{end+1,1} = [prep(k).x, prep(k).y]; %#ok<AGROW>

        if prep(k).hasActivation
            microRescueActivationEvent_xy(end+1,:) = [prep(k).actX, prep(k).actY]; %#ok<AGROW>
            microRescueActivationEventFrames(end+1,1) = prep(k).actFrame; %#ok<AGROW>
            microRescueActivationEventTrackIds(end+1,1) = prep(k).TRACK_ID; %#ok<AGROW>
            microRescueActivationSeedArea_px2(end+1,1) = startArea; %#ok<AGROW>
        end
    end
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
        A_over_I = strictActivationEventsTotal / max(nStrictPrimaryTracks, 1);
        [A_over_I_ci_low, A_over_I_ci_high] = wilson_ci(strictActivationEventsTotal, nStrictPrimaryTracks, 0.95);

        nNetLeftTracks = numel(netLeftTrackIds);
        netLeftActTrackIdsFinite = netLeftActivationEventTrackIds(isfinite(netLeftActivationEventTrackIds));
        nNetLeftActivatedTracks = numel(unique(netLeftActTrackIdsFinite));
        if nNetLeftActivatedTracks == 0
            nNetLeftActivatedTracks = size(netLeftActivationEvent_xy, 1);
        end
        netLeftTrackFrameExposure = sum(netLeftFrameVisible);
        netLeftActivationEventsTotal = sum(netLeftFrameActEvents);
        A_over_I_netLeftLegacy = netLeftActivationEventsTotal / max(nNetLeftTracks, 1);
        [A_over_I_ci_low_netLeftLegacy, A_over_I_ci_high_netLeftLegacy] = ...
            wilson_ci(netLeftActivationEventsTotal, nNetLeftTracks, 0.95);

        microTrackIdsFinite = microRescueTrackIds(isfinite(microRescueTrackIds));
        nMicroRescueTracks = numel(unique(microTrackIdsFinite));
        nMicroRescueActivationEvents = size(microRescueActivationEvent_xy, 1);

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

        % Non-strict microbubble-start outputs (for diagnostics/visualization).
        outMetrics.nMicrobubbleRescueTracks_nonLeft = nMicroRescueTracks;
        outMetrics.nMicrobubbleRescueActivationEvents_nonLeft = nMicroRescueActivationEvents;
        outMetrics.microbubbleRescueTrackIds_nonLeft = unique(microTrackIdsFinite, 'stable');
        outMetrics.microbubbleRescueTrack_xy_nonLeft = microRescueTrack_xy;
        outMetrics.microbubbleActivationEvent_xy_nonLeft = microRescueActivationEvent_xy;
        outMetrics.microbubbleActivationEvent_frame_nonLeft = microRescueActivationEventFrames;
        outMetrics.microbubbleActivationEvent_trackId_nonLeft = microRescueActivationEventTrackIds;
        outMetrics.microbubbleActivationSeedArea_px2_nonLeft = microRescueActivationSeedArea_px2;
        outMetrics.microbubbleActivationAreaRange_px2 = activationOpts.microbubbleStartAreaRange_px2(:).';

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
    'startArea', NaN, ...
    'isOriginExcludedBox', false, ...
    'isShort', false, ...
    'isNonFinite', false, ...
    'isNonMonotonicTime', false, ...
    'isExcessiveYStep', false, ...
    'isBasicValid', false, ...
    'isLeftMoving', false, ...
    'idxJump', NaN, ...
    'hasActivation', false, ...
    'actIdx', NaN, ...
    'actX', NaN, ...
    'actY', NaN, ...
    'actFrame', NaN, ...
    'actSeedArea', NaN, ...
    'isProximityRejected', false, ...
    'isUnwantedArea', false);
end

function tpl = make_catalog_template()
tpl = struct( ...
    'TRACK_ID', NaN, ...
    'frame', nan(0,1), ...
    'x', nan(0,1), ...
    'y', nan(0,1), ...
    'isBasicValid', false, ...
    'isLeftMoving', false, ...
    'isOriginExcludedBox', false, ...
    'isStrictPrimary', false, ...
    'isStrictActivated', false, ...
    'isMicrobubbleRescueNonLeft', false, ...
    'isActivated', false, ...
    'activationFrame', NaN, ...
    'activationX', NaN, ...
    'activationY', NaN, ...
    'activationIndex', NaN, ...
    'microbubbleSeedArea_px2', NaN, ...
    'strictActivationFrame', NaN, ...
    'strictActivationX', NaN, ...
    'strictActivationY', NaN, ...
    'strictActivationIndex', NaN);
end

function prep = build_track_prep(thisTraj, spotRowById, spots, qcOpts, flowOpts, activationOpts)
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

if is_origin_excluded_start(prep.x(1), prep.y(1), qcOpts)
    prep.isOriginExcludedBox = true;
    return;
end

if prep.isShort
    return;
end

% Reject tracks with excessive per-step y-displacement (mis-linked spots)
if qcOpts.maxStepDy_mm > 0
    dySteps = abs(diff(prep.y));
    if any(dySteps > qcOpts.maxStepDy_mm)
        prep.isExcessiveYStep = true;
        return;
    end
end

prep.isBasicValid = true;
prep.isLeftMoving = passes_netleft_motion_gate(prep.x, flowOpts);

prep.areaVals = extract_area_values(prep.spotIds, spotRowById, spots);
prep.startArea = estimate_track_start_area(prep.areaVals);
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
        prep.actSeedArea = estimate_activation_seed_area(prep.areaVals, idxJump, activationOpts);
    end
end

% Check if any spot falls inside the unwanted track area mask.
if ~isempty(qcOpts.unwantedAreaMask) && qcOpts.unwantedAreaMaskPixelSize > 0
    prep.isUnwantedArea = any_spot_in_mask(prep.x, prep.y, ...
        qcOpts.unwantedAreaMask, qcOpts.unwantedAreaMaskPixelSize);
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
catalog.isOriginExcludedBox = prep.isOriginExcludedBox;
catalog.isStrictPrimary = false;
catalog.isStrictActivated = false;
catalog.isMicrobubbleRescueNonLeft = false;
catalog.isActivated = prep.hasActivation;
catalog.activationFrame = prep.actFrame;
catalog.activationX = prep.actX;
catalog.activationY = prep.actY;
catalog.activationIndex = prep.actIdx;
catalog.microbubbleSeedArea_px2 = prep.startArea;
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

function startArea = estimate_track_start_area(areaVals)
startArea = NaN;
if isempty(areaVals)
    return;
end
idx = find(isfinite(areaVals) & areaVals > 0, 1, 'first');
if isempty(idx)
    return;
end
startArea = areaVals(idx);
end

function seedArea = estimate_activation_seed_area(areaVals, idxJump, opts)
seedArea = NaN;

if isempty(areaVals) || ~isfinite(idxJump)
    return;
end
n = numel(areaVals);
i = round(idxJump);
if i < 1 || i > n
    return;
end

preStart = max(1, i - opts.preWindowFrames + 1);
preVals = areaVals(preStart:i);
preVals = preVals(isfinite(preVals) & preVals > 0);
if numel(preVals) >= opts.minPrePoints
    seedArea = median(preVals);
    return;
end

aCurr = areaVals(i);
if isfinite(aCurr) && aCurr > 0
    seedArea = aCurr;
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

function tf = is_true_field(s, fieldName)
tf = false;
if ~isstruct(s) || ~isfield(s, fieldName)
    return;
end
v = s.(fieldName);
if islogical(v)
    tf = any(v(:));
elseif isnumeric(v)
    tf = any(v(:) ~= 0);
end
end

function tf = is_origin_excluded_start(xStart, yStart, qcOpts)
tf = false;
if ~isfield(qcOpts, 'excludeOriginBoxEnabled') || ~logical(qcOpts.excludeOriginBoxEnabled)
    return;
end
if ~(isfinite(xStart) && isfinite(yStart))
    return;
end
if ~isfield(qcOpts, 'excludeOriginBoxX_mm') || numel(qcOpts.excludeOriginBoxX_mm) < 2
    return;
end
if ~isfield(qcOpts, 'excludeOriginBoxY_mm') || numel(qcOpts.excludeOriginBoxY_mm) < 2
    return;
end

xLim = sort(double(qcOpts.excludeOriginBoxX_mm(1:2)));
yLim = sort(double(qcOpts.excludeOriginBoxY_mm(1:2)));
if any(~isfinite(xLim)) || any(~isfinite(yLim))
    return;
end

tf = (xStart >= xLim(1)) && (xStart <= xLim(2)) && ...
     (yStart >= yLim(1)) && (yStart <= yLim(2));
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

function pass = passes_netleft_motion_gate(x, flowOpts)
pass = false;

if numel(x) < 2
    return;
end

dx = diff(x);
dx = dx(isfinite(dx));
if isempty(dx)
    return;
end

minRunSteps = max(1, round(flowOpts.netLeftMinConsecutiveNegSteps));
minRunNetDx = max(0, double(flowOpts.netLeftMinConsecutiveNetDx_mm));

if strcmpi(flowOpts.bulkDirection, 'left_to_right')
    isCounterflowStep = (dx < 0);
    counterflowStepMag = -dx;
else
    isCounterflowStep = (dx > 0);
    counterflowStepMag = dx;
end

if ~any(isCounterflowStep)
    return;
end

runEdge = diff([false; isCounterflowStep(:); false]);
runStart = find(runEdge == 1);
runStop = find(runEdge == -1) - 1;

for r = 1:numel(runStart)
    s = runStart(r);
    e = runStop(r);
    if (e - s + 1) < minRunSteps
        continue;
    end

    runNetDx = sum(counterflowStepMag(s:e), 'omitnan');
    if isfinite(runNetDx) && runNetDx >= minRunNetDx
        pass = true;
        return;
    end
end
end

function pass = passes_netleft_band_rescue(prep, flowOpts)
pass = false;
if ~isfield(flowOpts, 'netLeftBandRescueEnabled') || ~flowOpts.netLeftBandRescueEnabled
    return;
end

xBand = nan(1,2);
if isfield(flowOpts, 'netLeftBandRescueX_mm') && numel(flowOpts.netLeftBandRescueX_mm) >= 2
    xBand = double(flowOpts.netLeftBandRescueX_mm(1:2));
end
if any(~isfinite(xBand))
    return;
end
xBand = sort(xBand);

if numel(prep.x) < 2
    return;
end

if flowOpts.netLeftBandRescueRequireActivation
    if ~prep.hasActivation
        return;
    end
    if flowOpts.netLeftBandRescueRequireActXInBand
        if ~(isfinite(prep.actX) && prep.actX >= xBand(1) && prep.actX <= xBand(2))
            return;
        end
    end
end

dx = diff(prep.x);
x0 = prep.x(1:end-1);
x1 = prep.x(2:end);
finiteStep = isfinite(dx) & isfinite(x0) & isfinite(x1);
if ~any(finiteStep)
    return;
end

if strcmpi(flowOpts.bulkDirection, 'left_to_right')
    isCounterflowStep = (dx < 0);
    counterflowStepMag = -dx;
else
    isCounterflowStep = (dx > 0);
    counterflowStepMag = dx;
end

inBandStep = (min(x0, x1) <= xBand(2)) & (max(x0, x1) >= xBand(1));
cand = finiteStep & isCounterflowStep & inBandStep;

if ~any(cand)
    return;
end

minSteps = max(1, round(flowOpts.netLeftBandRescueMinCounterflowSteps));
minDx = max(0, double(flowOpts.netLeftBandRescueMinCounterflowDx_mm));

if sum(cand) < minSteps
    return;
end

sumDx = sum(counterflowStepMag(cand), 'omitnan');
if ~(isfinite(sumDx) && sumDx >= minDx)
    return;
end

pass = true;
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

function [prep, nRejected] = reject_proximity_merge_activations(prep, spots, activationOpts)
%REJECT_PROXIMITY_MERGE_ACTIVATIONS  Reject activations caused by overlapping blobs.
%   At the pre-jump frame, check if any OTHER spot is within
%   mergeProximityRadius_mm. If so, the area jump is likely a thresholding
%   artifact from two bubbles overlapping in the binary image.
nRejected = 0;

if ~activationOpts.rejectProximityMerge || activationOpts.mergeProximityRadius_mm <= 0
    return;
end

% Build per-frame spot position index for efficient lookup.
hasFrameCol = ismember('FRAME', spots.Properties.VariableNames);
hasXCol = ismember('X', spots.Properties.VariableNames);
hasYCol = ismember('Y', spots.Properties.VariableNames);
hasIDCol = ismember('ID', spots.Properties.VariableNames);
if ~hasFrameCol || ~hasXCol || ~hasYCol || ~hasIDCol
    return;
end

% Get pixelSize to convert spot pixel coords to mm (prep.x/y are in mm).
if isfield(activationOpts, 'pixelSize') && isfinite(activationOpts.pixelSize) && activationOpts.pixelSize > 0
    pxSz = activationOpts.pixelSize;
else
    pxSz = 1;  % fallback: assume spots are already in physical units
end

spotFrames = spots.FRAME(:);
spotX = spots.X(:) * pxSz;  % convert pixels -> mm
spotY = spots.Y(:) * pxSz;  % convert pixels -> mm
spotIDs = spots.ID(:);
uniqueFrames = unique(spotFrames(isfinite(spotFrames)));
frameSpotIdx = containers.Map('KeyType', 'double', 'ValueType', 'any');
for fi = 1:numel(uniqueFrames)
    fr = uniqueFrames(fi);
    frameSpotIdx(fr) = find(spotFrames == fr);
end

radius = activationOpts.mergeProximityRadius_mm;
radiusSq = radius * radius;

for k = 1:numel(prep)
    if ~prep(k).hasActivation || ~isfinite(prep(k).idxJump)
        continue;
    end

    preIdx = prep(k).idxJump;
    if preIdx < 1 || preIdx > numel(prep(k).frame) || preIdx > numel(prep(k).spotIds)
        continue;
    end

    preFrame = prep(k).frame(preIdx);
    preSpotId = prep(k).spotIds(preIdx);
    preX = prep(k).x(preIdx);
    preY = prep(k).y(preIdx);

    if ~isfinite(preFrame) || ~isfinite(preX) || ~isfinite(preY)
        continue;
    end

    if ~isKey(frameSpotIdx, preFrame)
        continue;
    end

    rowsAtFrame = frameSpotIdx(preFrame);
    foundNearby = false;
    for ri = 1:numel(rowsAtFrame)
        row = rowsAtFrame(ri);
        if spotIDs(row) == preSpotId
            continue;
        end
        dx = spotX(row) - preX;
        dy = spotY(row) - preY;
        if (dx * dx + dy * dy) < radiusSq
            foundNearby = true;
            break;
        end
    end

    if foundNearby
        prep(k).hasActivation = false;
        prep(k).isProximityRejected = true;
        prep(k).actX = NaN;
        prep(k).actY = NaN;
        prep(k).actFrame = NaN;
        prep(k).actIdx = NaN;
        prep(k).actSeedArea = NaN;
        prep(k).idxJump = NaN;
        nRejected = nRejected + 1;
    end
end

if nRejected > 0
    fprintf('  Proximity merge rejection: %d activation(s) rejected (radius=%.4f mm)\n', nRejected, radius);
end
end

function idxJump = find_sustained_growth_activation(areaVals, opts)
%FIND_SUSTAINED_GROWTH_ACTIVATION  Detect first cavitation activation event.
%   Uses absolute area floors + sustained monotonic growth instead of a
%   fixed jump-ratio threshold.
idxJump = [];
n = numel(areaVals);

if n < 2
    return;
end

for i = 1:(n - 1)
    aCurr = areaVals(i);
    aNext = areaVals(i + 1);

    % Stage 0: basic validity
    if ~(isfinite(aCurr) && isfinite(aNext) && aCurr > 0 && aNext > 0)
        continue;
    end

    % Stage 1: absolute area floors (eliminates out-of-focus artifacts)
    if aCurr < opts.minPreJumpArea_px2
        continue;
    end
    if aNext < opts.minPostJumpArea_px2
        continue;
    end

    % Stage 2: initial growth sanity check — aNext vs aCurr
    aPre = aCurr;
    if (aNext / aPre) < opts.minInitialGrowthRatio
        continue;
    end

    % Stage 3: sustained monotonic growth window
    % Large post-jump area (>= threshold) uses a shorter sustained window
    if aNext >= opts.largeAreaThreshold_px2
        swFrames = opts.largeAreaSustainedFrames;
    else
        swFrames = opts.sustainedWindowFrames;
    end
    winEnd = min(n, i + 1 + swFrames);
    winVals = areaVals((i + 1):winEnd);
    winVals_valid = winVals(isfinite(winVals) & winVals > 0);

    if numel(winVals_valid) >= 2
        nPairs = numel(winVals_valid) - 1;
        nPassing = 0;
        for j = 1:nPairs
            ratio_j = winVals_valid(j + 1) / winVals_valid(j);
            delta_j = winVals_valid(j + 1) - winVals_valid(j);
            if ratio_j >= opts.sustainedMinRatio || delta_j >= opts.sustainedMinDelta_px2
                nPassing = nPassing + 1;
            end
        end

        if nPairs >= opts.sustainedMinPairsAvailable && nPassing >= opts.sustainedMinPassingPairs
            if passes_post_activation_area_check(areaVals, i, n, opts)
                idxJump = i;
                return;
            else
                continue;
            end
        elseif nPairs < opts.sustainedMinPairsAvailable
            % Not enough pairs — fall through to burst fallback (Stage 4)
        else
            % Enough pairs but too few passed — reject this candidate
            continue;
        end
    end

    % Stage 4: burst fallback (track/video ends before sustained window)
    if opts.enableBurstFallback
        if (aNext / aPre) >= opts.burstMinGrowthRatio && aNext >= opts.burstMinArea_px2
            if passes_post_activation_area_check(areaVals, i, n, opts)
                idxJump = i;
                return;
            else
                continue;
            end
        end
    end
end
end

function ok = passes_post_activation_area_check(areaVals, i, n, opts)
%PASSES_POST_ACTIVATION_AREA_CHECK  Verify that the next postActivationCheckFrames
%   frames after the jump each have area >= postActivationMinArea_px2.
%   If not enough frames remain, check whatever is available.
ok = true;
nCheck = opts.postActivationCheckFrames;
minA = opts.postActivationMinArea_px2;
for j = 1:nCheck
    idx = i + 1 + j;  % i+1 is the jump frame; check i+2, i+3, ...
    if idx > n
        break;  % track ends — accept with whatever we've seen
    end
    a = areaVals(idx);
    if ~isfinite(a) || a < minA
        ok = false;
        return;
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
    activationOpts.minInitialGrowthRatio = max(1, double(varargin{3}));
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
if ~isfield(qcOpts, 'excludeOriginBoxEnabled') || isempty(qcOpts.excludeOriginBoxEnabled)
    qcOpts.excludeOriginBoxEnabled = false;
end
if ~isfield(qcOpts, 'excludeOriginBoxX_mm') || numel(qcOpts.excludeOriginBoxX_mm) < 2
    qcOpts.excludeOriginBoxX_mm = [0 0.5];
end
if ~isfield(qcOpts, 'excludeOriginBoxY_mm') || numel(qcOpts.excludeOriginBoxY_mm) < 2
    qcOpts.excludeOriginBoxY_mm = [0 1.2];
end
qcOpts.excludeOriginBoxX_mm = sort(double(qcOpts.excludeOriginBoxX_mm(1:2)));
qcOpts.excludeOriginBoxY_mm = sort(double(qcOpts.excludeOriginBoxY_mm(1:2)));
if any(~isfinite(qcOpts.excludeOriginBoxX_mm))
    qcOpts.excludeOriginBoxX_mm = [0 0.5];
end
if any(~isfinite(qcOpts.excludeOriginBoxY_mm))
    qcOpts.excludeOriginBoxY_mm = [0 1.2];
end
flowOpts.minNetDxCounterflow_mm = max(0, flowOpts.minNetDxCounterflow_mm);
flowOpts.minNegativeStepFraction = min(1, max(0, flowOpts.minNegativeStepFraction));
flowOpts.maxPositiveStepFraction = min(1, max(0, flowOpts.maxPositiveStepFraction));
flowOpts.rightOriginFrac = min(1, max(0, flowOpts.rightOriginFrac));
flowOpts.netLeftMinConsecutiveNegSteps = max(1, round(flowOpts.netLeftMinConsecutiveNegSteps));
flowOpts.netLeftMinConsecutiveNetDx_mm = max(0, flowOpts.netLeftMinConsecutiveNetDx_mm);
if ~isfield(flowOpts, 'netLeftBandRescueEnabled') || isempty(flowOpts.netLeftBandRescueEnabled)
    flowOpts.netLeftBandRescueEnabled = false;
end
if ~isfield(flowOpts, 'netLeftBandRescueX_mm') || numel(flowOpts.netLeftBandRescueX_mm) < 2
    flowOpts.netLeftBandRescueX_mm = [0.5 1.5];
end
if ~isfield(flowOpts, 'netLeftBandRescueMinCounterflowSteps') || isempty(flowOpts.netLeftBandRescueMinCounterflowSteps)
    flowOpts.netLeftBandRescueMinCounterflowSteps = 1;
end
if ~isfield(flowOpts, 'netLeftBandRescueMinCounterflowDx_mm') || isempty(flowOpts.netLeftBandRescueMinCounterflowDx_mm)
    flowOpts.netLeftBandRescueMinCounterflowDx_mm = 0.0;
end
if ~isfield(flowOpts, 'netLeftBandRescueRequireActivation') || isempty(flowOpts.netLeftBandRescueRequireActivation)
    flowOpts.netLeftBandRescueRequireActivation = true;
end
if ~isfield(flowOpts, 'netLeftBandRescueRequireActXInBand') || isempty(flowOpts.netLeftBandRescueRequireActXInBand)
    flowOpts.netLeftBandRescueRequireActXInBand = false;
end
if ~isfield(flowOpts, 'netLeftBandRescueBypassOriginGate') || isempty(flowOpts.netLeftBandRescueBypassOriginGate)
    flowOpts.netLeftBandRescueBypassOriginGate = true;
end
flowOpts.netLeftBandRescueMinCounterflowSteps = max(1, round(flowOpts.netLeftBandRescueMinCounterflowSteps));
flowOpts.netLeftBandRescueMinCounterflowDx_mm = max(0, flowOpts.netLeftBandRescueMinCounterflowDx_mm);
activationOpts.minPreJumpArea_px2 = max(0, activationOpts.minPreJumpArea_px2);
activationOpts.minPostJumpArea_px2 = max(0, activationOpts.minPostJumpArea_px2);
activationOpts.preWindowFrames = max(1, round(activationOpts.preWindowFrames));
activationOpts.minPrePoints = max(1, round(activationOpts.minPrePoints));
activationOpts.minInitialGrowthRatio = max(1, activationOpts.minInitialGrowthRatio);
activationOpts.sustainedWindowFrames = max(1, round(activationOpts.sustainedWindowFrames));
activationOpts.sustainedMinPassingPairs = max(1, round(activationOpts.sustainedMinPassingPairs));
activationOpts.sustainedMinPairsAvailable = max(1, round(activationOpts.sustainedMinPairsAvailable));
if ~isfield(activationOpts, 'includeMicrobubbleActivationRescue') || isempty(activationOpts.includeMicrobubbleActivationRescue)
    activationOpts.includeMicrobubbleActivationRescue = false;
end
if ~isfield(activationOpts, 'microbubbleStartAreaRange_px2') || numel(activationOpts.microbubbleStartAreaRange_px2) < 2
    if isfield(activationOpts, 'microbubbleSeedAreaRange_px2') && numel(activationOpts.microbubbleSeedAreaRange_px2) >= 2
        activationOpts.microbubbleStartAreaRange_px2 = activationOpts.microbubbleSeedAreaRange_px2(1:2);
    else
        activationOpts.microbubbleStartAreaRange_px2 = [1 120];
    end
end
startRange = double(activationOpts.microbubbleStartAreaRange_px2(1:2));
if any(~isfinite(startRange))
    startRange = [1 120];
end
startRange = sort(startRange);
startRange(startRange < 0) = 0;
activationOpts.microbubbleStartAreaRange_px2 = startRange;

if ~isfield(activationOpts, 'microbubbleRequireOutsideStrictPrimary') || isempty(activationOpts.microbubbleRequireOutsideStrictPrimary)
    if isfield(activationOpts, 'microbubbleRequireNonLeftMoving') && ~isempty(activationOpts.microbubbleRequireNonLeftMoving)
        activationOpts.microbubbleRequireOutsideStrictPrimary = logical(activationOpts.microbubbleRequireNonLeftMoving);
    else
        activationOpts.microbubbleRequireOutsideStrictPrimary = true;
    end
end

if ~isfield(activationOpts, 'rejectProximityMerge') || isempty(activationOpts.rejectProximityMerge)
    activationOpts.rejectProximityMerge = true;
end
if ~isfield(activationOpts, 'mergeProximityRadius_mm') || isempty(activationOpts.mergeProximityRadius_mm)
    activationOpts.mergeProximityRadius_mm = 0.05;
end
activationOpts.mergeProximityRadius_mm = max(0, activationOpts.mergeProximityRadius_mm);

% Keep legacy aliases populated for downstream compatibility.
activationOpts.microbubbleSeedAreaRange_px2 = activationOpts.microbubbleStartAreaRange_px2;
activationOpts.microbubbleRequireNonLeftMoving = activationOpts.microbubbleRequireOutsideStrictPrimary;
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
qcOpts.minTrackSpots = 4;
qcOpts.maxTrackGaps = 1;
qcOpts.maxLeftMovingTracks = Inf;
qcOpts.rejectSplitMergeComplex = true;
qcOpts.wallBandEnabled = false;
qcOpts.wallBandYLimits_mm = [];
qcOpts.excludeOriginBoxEnabled = false;
qcOpts.excludeOriginBoxX_mm = [0 0.5];
qcOpts.excludeOriginBoxY_mm = [0 1.2];
qcOpts.maxStepDy_mm = 0.1;
qcOpts.unwantedAreaMask = [];
qcOpts.unwantedAreaMaskPixelSize = 0;
end

function flowOpts = default_flow_opts()
flowOpts = struct();
flowOpts.bulkDirection = "left_to_right";
flowOpts.minNetDxCounterflow_mm = 0.08;
flowOpts.minNegativeStepFraction = 0.65;
flowOpts.maxPositiveStepFraction = 0.30;
flowOpts.requireRightOrigin = true;
flowOpts.rightOriginFrac = 0.60;
flowOpts.netLeftMinConsecutiveNegSteps = 3;
flowOpts.netLeftMinConsecutiveNetDx_mm = 0.0;
flowOpts.netLeftBandRescueEnabled = false;
flowOpts.netLeftBandRescueX_mm = [0.5 1.5];
flowOpts.netLeftBandRescueMinCounterflowSteps = 1;
flowOpts.netLeftBandRescueMinCounterflowDx_mm = 0.0;
flowOpts.netLeftBandRescueRequireActivation = true;
flowOpts.netLeftBandRescueRequireActXInBand = false;
flowOpts.netLeftBandRescueBypassOriginGate = true;
end

function activationOpts = default_activation_opts()
activationOpts = struct();
activationOpts.minPreJumpArea_px2 = 3;
activationOpts.minPostJumpArea_px2 = 20;
activationOpts.preWindowFrames = 2;
activationOpts.minPrePoints = 1;
activationOpts.minInitialGrowthRatio = 2.0;
activationOpts.sustainedWindowFrames = 3;
activationOpts.largeAreaThreshold_px2 = 150;
activationOpts.largeAreaSustainedFrames = 2;
activationOpts.sustainedMinRatio = 1.05;
activationOpts.sustainedMinDelta_px2 = 5;
activationOpts.sustainedMinPassingPairs = 2;
activationOpts.sustainedMinPairsAvailable = 2;
activationOpts.enableBurstFallback = true;
activationOpts.burstMinGrowthRatio = 3.0;
activationOpts.burstMinArea_px2 = 40;
activationOpts.postActivationCheckFrames = 2;
activationOpts.postActivationMinArea_px2 = 100;
activationOpts.rejectProximityMerge = true;
activationOpts.mergeProximityRadius_mm = 0.05;
activationOpts.pixelSize = 1;
activationOpts.includeMicrobubbleActivationRescue = false;
activationOpts.microbubbleStartAreaRange_px2 = [1 120];
activationOpts.microbubbleRequireOutsideStrictPrimary = true;
activationOpts.microbubbleSeedAreaRange_px2 = activationOpts.microbubbleStartAreaRange_px2;
activationOpts.microbubbleRequireNonLeftMoving = activationOpts.microbubbleRequireOutsideStrictPrimary;
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

function result = any_spot_in_mask(x_mm, y_mm, mask, pixelSize)
% Return true if ANY spot position (x,y in image mm coords) falls inside mask.
[nRows, nCols] = size(mask);
result = false;
for i = 1:numel(x_mm)
    c = round(x_mm(i) / pixelSize);
    r = round(y_mm(i) / pixelSize);
    if isfinite(c) && isfinite(r) && c >= 1 && c <= nCols && r >= 1 && r <= nRows
        if mask(r, c)
            result = true;
            return;
        end
    end
end
end
