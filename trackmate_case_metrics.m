function metrics = trackmate_case_metrics(out, varargin)
%TRACKMATE_CASE_METRICS  Robust injected/activated metrics for TrackMate trajectories.
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

injInception_xy = zeros(0,2);
actLocation_xy  = zeros(0,2);
inception2x_xy  = zeros(0,2);
upstreamTrack_xy = cell(0,1);

tau_values = nan(0,1);
upstreamSize_eqd = nan(0,1);

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

if isempty(traj) || isempty(spots) || ~ismember('ID', spots.Properties.VariableNames)
    metrics = pack_metrics();
    return;
end

% Map SpotID -> row
[spotIdSorted, sortIdx] = sort(spots.ID);
spotRowById = containers.Map('KeyType', 'double', 'ValueType', 'double');
for k = 1:numel(spotIdSorted)
    spotRowById(spotIdSorted(k)) = sortIdx(k);
end

% Map TRACK_ID -> row in tracks table (for topology QC).
trackRowById = containers.Map('KeyType', 'double', 'ValueType', 'double');
if ~isempty(tracks) && ismember('TRACK_ID', tracks.Properties.VariableNames)
    for k = 1:height(tracks)
        tid = tracks.TRACK_ID(k);
        if isfinite(tid) && ~isKey(trackRowById, tid)
            trackRowById(tid) = k;
        end
    end
end

xStarts = nan(0,1);
for k = 1:nTotal
    xk = traj(k).x_phys(:);
    if ~isempty(xk) && isfinite(xk(1))
        xStarts(end+1,1) = xk(1); %#ok<AGROW>
    end
end
[originThreshold, xMinAll, xMaxAll] = compute_origin_threshold(xStarts, flowOpts);
gateStats.originThreshold = originThreshold;
gateStats.xStartMin = xMinAll;
gateStats.xStartMax = xMaxAll;

for k = 1:nTotal
    nTrackSpots = numel(traj(k).spotIds);
    if nTrackSpots < qcOpts.minTrackSpots
        gateStats.nRejectedTooShort = gateStats.nRejectedTooShort + 1;
        continue;
    end

    x = traj(k).x_phys(:);
    y = traj(k).y_phys(:);
    t = traj(k).t(:);

    if numel(x) < 2 || numel(y) ~= numel(x) || numel(t) ~= numel(x) || any(~isfinite([x; y; t]))
        gateStats.nRejectedNonFinite = gateStats.nRejectedNonFinite + 1;
        continue;
    end

    dtTrack = diff(t);
    if isempty(dtTrack) || any(~isfinite(dtTrack)) || any(dtTrack <= 0)
        gateStats.nRejectedNonMonotonicTime = gateStats.nRejectedNonMonotonicTime + 1;
        continue;
    end

    if ~passes_topology_gate(traj(k), trackRowById, tracks, qcOpts)
        gateStats.nRejectedTopology = gateStats.nRejectedTopology + 1;
        continue;
    end

    [isCounterflow, ~, ~] = passes_counterflow_gate(x, flowOpts);
    if ~isCounterflow
        gateStats.nRejectedFlow = gateStats.nRejectedFlow + 1;
        continue;
    end

    if flowOpts.requireRightOrigin && isfinite(originThreshold)
        if ~passes_origin_gate(x(1), originThreshold, flowOpts)
            gateStats.nRejectedOrigin = gateStats.nRejectedOrigin + 1;
            continue;
        end
    end

    gateStats.nInjected = gateStats.nInjected + 1;

    % Inception location = first spot of this injected track.
    injInception_xy(end+1,:) = [x(1), y(1)]; %#ok<AGROW>
    upstreamTrack_xy{end+1,1} = [x, y]; %#ok<AGROW>

    areaVals = nan(numel(traj(k).spotIds),1);
    for ii = 1:numel(traj(k).spotIds)
        sid = traj(k).spotIds(ii);
        if isKey(spotRowById, sid)
            r = spotRowById(sid);
            if ismember('AREA', spots.Properties.VariableNames)
                areaVals(ii) = spots.AREA(r);
            end
        end
    end

    idxJump = find_sustained_growth_activation(areaVals, activationOpts);

    if isempty(idxJump)
        gateStats.nRejectedNoActivation = gateStats.nRejectedNoActivation + 1;
        nUse = numel(areaVals);
    else
        nUse = idxJump;
    end

    if nUse >= 1
        areaUse = areaVals(1:nUse);
        eqd = sqrt(4 .* areaUse ./ pi);
        eqd = eqd(isfinite(eqd) & eqd > 0);
        if ~isempty(eqd)
            upstreamSize_eqd = [upstreamSize_eqd; eqd(:)]; %#ok<AGROW>
        end
    end

    if isempty(idxJump)
        continue;
    end

    actX = x(idxJump+1);
    actY = y(idxJump+1);

    if qcOpts.wallBandEnabled && ~passes_wall_band(actY, qcOpts.wallBandYLimits_mm)
        gateStats.nRejectedWallBand = gateStats.nRejectedWallBand + 1;
        continue;
    end

    gateStats.nActivated = gateStats.nActivated + 1;
    actLocation_xy(end+1,:) = [actX, actY]; %#ok<AGROW>
    inception2x_xy(end+1,:) = [actX, actY]; %#ok<AGROW>

    if numel(t) >= (idxJump + 1) && isfinite(t(idxJump+1)) && isfinite(t(1))
        tau_values(end+1,1) = max(t(idxJump+1) - t(1), 0); %#ok<AGROW>
    end
end

metrics = pack_metrics();

    function outMetrics = pack_metrics()
        nInjected = gateStats.nInjected;
        nActivated = gateStats.nActivated;
        A_over_I = nActivated / max(nInjected, 1);
        [A_over_I_ci_low, A_over_I_ci_high] = wilson_ci(nActivated, nInjected, 0.95);

        outMetrics = struct();
        outMetrics.nTracksTotal     = nTotal;
        outMetrics.nInjected        = nInjected;
        outMetrics.nActivated       = nActivated;
        outMetrics.A_over_I         = A_over_I;
        outMetrics.A_over_I_ci_low  = A_over_I_ci_low;
        outMetrics.A_over_I_ci_high = A_over_I_ci_high;
        outMetrics.injInception_xy  = injInception_xy;
        outMetrics.actLocation_xy   = actLocation_xy;
        outMetrics.inception2x_xy   = inception2x_xy;
        outMetrics.tau_values       = tau_values;
        outMetrics.tau_mean         = mean(tau_values, 'omitnan');
        outMetrics.tau_std          = std(tau_values, 0, 'omitnan');
        outMetrics.upstreamSize_eqd = upstreamSize_eqd;
        outMetrics.upstreamTrack_xy = upstreamTrack_xy;
        outMetrics.gateStats        = gateStats;
        outMetrics.qcOpts           = qcOpts;
        outMetrics.flowOpts         = flowOpts;
        outMetrics.activationOpts   = activationOpts;
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
