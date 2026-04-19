function lagData = analyze_lagrangian_acceleration(out, metrics, caseDef, opts)
%ANALYZE_LAGRANGIAN_ACCELERATION  Smoothed track acceleration proxy analysis.
%
%   Computes raw and Savitzky-Golay-style smoothed Lagrangian acceleration
%   along TrackMate trajectories. Acceleration is normalized as
%       a* = |a| d_mean / U_ref^2
%   where d_mean is the case-mean equivalent bubble diameter and U_ref is
%   the case-mean upstream axial speed from the selected track population.

if nargin < 4 || isempty(opts)
    opts = struct();
end
opts = apply_defaults(opts);

lagData = make_case_result_template();
lagData.caseName = string(caseDef.name);
lagData.Re = caseDef.Re;
lagData.kD = caseDef.kD;
lagData.pixelSize_mm = caseDef.pixelSize;
lagData.opts = opts;

if isempty(out) || ~isstruct(out) || ~isfield(out, 'trajectories') || isempty(out.trajectories) || ...
        isempty(metrics) || ~isstruct(metrics) || ~isfield(metrics, 'trackCatalog') || isempty(metrics.trackCatalog)
    lagData.summary = summarize_case(lagData);
    return;
end

catalog = metrics.trackCatalog(:);
catMap = build_catalog_map(catalog);
areaMap = build_spot_area_map(out);

[dPool_m, uPool_m_s] = collect_normalization_pools(out.trajectories(:), catalog, catMap, areaMap, caseDef, opts);
lagData.dMean_m = finite_mean(dPool_m);
lagData.U_ref_m_s = finite_mean(uPool_m_s);
if ~(isfinite(lagData.dMean_m) && lagData.dMean_m > 0)
    lagData.dMean_m = opts.fallbackDiameter_m;
end
if ~(isfinite(lagData.U_ref_m_s) && lagData.U_ref_m_s > 0)
    lagData.U_ref_m_s = opts.fallbackURef_m_s;
end

yExtent_mm = opts.imageSize_px(2) * caseDef.pixelSize;
halfWindow = floor(opts.sgWindowFrames / 2);
caseSeed = opts.randomSeed + round(caseDef.Re) + round(caseDef.kD * 1e6);
rng(caseSeed, 'twister');

tracks = repmat(make_track_result_template(), 0, 1);
sanityTracks = repmat(make_sanity_track_template(), 0, 1);
skip = make_skip_counts();

for ti = 1:numel(out.trajectories)
    tr = out.trajectories(ti);
    [cat, hasCatalog] = lookup_catalog_entry(tr, catalog, catMap);
    if ~hasCatalog
        skip.noCatalog = skip.noCatalog + 1;
        continue;
    end

    isMainTrack = track_in_population(cat, opts.trackPopulation);
    raw = extract_track_arrays(tr, areaMap, caseDef.pixelSize);
    if raw.n < opts.minTrackFrames
        skip.tooShort = skip.tooShort + 1;
        continue;
    end
    if opts.requireConsecutiveFrames && has_excessive_frame_gap(raw.frame, opts.maxAllowedFrameGap)
        skip.frameGaps = skip.frameGaps + 1;
        continue;
    end
    if ~raw.isUsable
        skip.nonFinite = skip.nonFinite + 1;
        continue;
    end

    isStationaryCandidate = opts.makeStationarySanityCheck && ...
        is_stationary_candidate(cat, raw, opts, caseDef.pixelSize);

    if ~isMainTrack && ~isStationaryCandidate
        continue;
    end

    [smoothX_m, ~, ax_m_s2, smoothValidX] = local_poly_kinematics(raw.t_s, raw.x_m, opts.sgWindowFrames, opts.sgPolyOrder);
    [smoothY_m, ~, ay_m_s2, smoothValidY] = local_poly_kinematics(raw.t_s, raw.y_m, opts.sgWindowFrames, opts.sgPolyOrder);
    [~, axRaw_m_s2] = finite_difference_kinematics(raw.t_s, raw.x_m);
    [~, ayRaw_m_s2] = finite_difference_kinematics(raw.t_s, raw.y_m);

    aMag_m_s2 = hypot(ax_m_s2, ay_m_s2);
    aRawMag_m_s2 = hypot(axRaw_m_s2, ayRaw_m_s2);
    aStar = aMag_m_s2 .* lagData.dMean_m ./ (lagData.U_ref_m_s ^ 2);
    rawAStar = aRawMag_m_s2 .* lagData.dMean_m ./ (lagData.U_ref_m_s ^ 2);

    validAccel = smoothValidX & smoothValidY & isfinite(aStar) & aStar >= 0;
    if opts.excludeSmoothingEdgeFrames && numel(validAccel) > 2 * halfWindow
        validAccel(1:halfWindow) = false;
        validAccel((end-halfWindow+1):end) = false;
    end
    validRawAccel = isfinite(rawAStar) & rawAStar >= 0;

    if isStationaryCandidate
        lagData.stationaryAstar = [lagData.stationaryAstar; aStar(validAccel)]; %#ok<AGROW>
        lagData.nStationaryTracks = lagData.nStationaryTracks + 1;
    end

    if ~isMainTrack
        continue;
    end

    if ~any(validAccel)
        skip.noValidAcceleration = skip.noValidAcceleration + 1;
        continue;
    end

    trackResult = make_track_result_template();
    trackResult.TRACK_ID = get_scalar_field(cat, 'TRACK_ID', NaN);
    trackResult.isActivated = get_bool_field(cat, 'isStrictActivated');
    trackResult.isStrictPrimary = get_bool_field(cat, 'isStrictPrimary');
    trackResult.nPoints = raw.n;
    trackResult.nAccelSamples = sum(validAccel);
    trackResult.activationIndex = get_scalar_field(cat, 'strictActivationIndex', NaN);
    if ~isfinite(trackResult.activationIndex)
        trackResult.activationIndex = get_scalar_field(cat, 'activationIndex', NaN);
    end
    trackResult.activationFrame = get_scalar_field(cat, 'strictActivationFrame', NaN);
    if ~isfinite(trackResult.activationFrame)
        trackResult.activationFrame = get_scalar_field(cat, 'activationFrame', NaN);
    end
    trackResult.growthRatio = compute_growth_ratio(raw.area_px2);
    trackResult.allAstar = aStar(validAccel);
    trackResult.rawAllAstar = rawAStar(validRawAccel);
    trackResult.xNorm = smoothX_m(validAccel) * 1000 ./ opts.throatHeight_mm;
    trackResult.yNorm = (yExtent_mm - smoothY_m(validAccel) * 1000) ./ opts.throatHeight_mm;

    if trackResult.isActivated && isfinite(trackResult.activationIndex)
        actIdx = round(trackResult.activationIndex);
        [triggerIdx, triggerOk] = trigger_indices_for_activation(actIdx, raw.frame, validAccel, opts);
        if triggerOk
            trackResult.triggerAstar = aStar(triggerIdx);
            trackResult.peakTriggerAstar = max(trackResult.triggerAstar);
            trackResult.nTriggerSamples = numel(trackResult.triggerAstar);
        end
    else
        [randomIdx, randomOk] = random_window_indices(raw.frame, validAccel, opts);
        if randomOk
            trackResult.randomWindowAstar = aStar(randomIdx);
            trackResult.nRandomWindowSamples = numel(trackResult.randomWindowAstar);
        end
    end

    if trackResult.isActivated
        actIdx = round(trackResult.activationIndex);
        if isfinite(actIdx) && actIdx >= 1 && actIdx <= raw.n
            trackResult.activationXY_mm = [raw.x_m(actIdx), raw.y_m(actIdx)] * 1000;
            trackResult.activationXYNorm = [trackResult.activationXY_mm(1) / opts.throatHeight_mm, ...
                (yExtent_mm - trackResult.activationXY_mm(2)) / opts.throatHeight_mm];
        end
    end

    tracks(end+1,1) = trackResult; %#ok<AGROW>
    lagData.sampleAstar = [lagData.sampleAstar; trackResult.allAstar(:)]; %#ok<AGROW>
    lagData.sampleRawAstar = [lagData.sampleRawAstar; trackResult.rawAllAstar(:)]; %#ok<AGROW>
    lagData.sampleXNorm = [lagData.sampleXNorm; trackResult.xNorm(:)]; %#ok<AGROW>
    lagData.sampleYNorm = [lagData.sampleYNorm; trackResult.yNorm(:)]; %#ok<AGROW>
    lagData.sampleIsActivatedTrack = [lagData.sampleIsActivatedTrack; ...
        repmat(trackResult.isActivated, numel(trackResult.allAstar), 1)]; %#ok<AGROW>

    if trackResult.isActivated
        lagData.activatedAllAstar = [lagData.activatedAllAstar; trackResult.allAstar(:)]; %#ok<AGROW>
        lagData.activatedTriggerAstar = [lagData.activatedTriggerAstar; trackResult.triggerAstar(:)]; %#ok<AGROW>
        if isfinite(trackResult.peakTriggerAstar) && isfinite(trackResult.growthRatio)
            lagData.peakTriggerAstar = [lagData.peakTriggerAstar; trackResult.peakTriggerAstar]; %#ok<AGROW>
            lagData.growthRatio = [lagData.growthRatio; trackResult.growthRatio]; %#ok<AGROW>
        end
        if all(isfinite(trackResult.activationXY_mm))
            lagData.activationXY_mm = [lagData.activationXY_mm; trackResult.activationXY_mm]; %#ok<AGROW>
            lagData.activationXYNorm = [lagData.activationXYNorm; trackResult.activationXYNorm]; %#ok<AGROW>
        end
    else
        lagData.nonActivatedAllAstar = [lagData.nonActivatedAllAstar; trackResult.allAstar(:)]; %#ok<AGROW>
        lagData.nonActivatedRandomAstar = [lagData.nonActivatedRandomAstar; trackResult.randomWindowAstar(:)]; %#ok<AGROW>
    end

    if opts.makeSanityPlots && numel(sanityTracks) < opts.nSanityTracks && rand() < 0.20
        sanityTracks(end+1,1) = make_sanity_track(tr, raw, smoothX_m, smoothY_m, trackResult, yExtent_mm, opts); %#ok<AGROW>
    end
end

if opts.makeSanityPlots && numel(sanityTracks) < opts.nSanityTracks
    sanityTracks = top_up_sanity_tracks(sanityTracks, tracks, out.trajectories(:), catalog, catMap, areaMap, caseDef, opts, yExtent_mm);
end

lagData.tracks = tracks;
lagData.sanityTracks = sanityTracks;
lagData.skipCounts = skip;
lagData.summary = summarize_case(lagData);
end


% =========================================================================
function opts = apply_defaults(opts)
opts = default_field(opts, 'sgWindowFrames', 7);
opts.sgWindowFrames = max(3, round(opts.sgWindowFrames));
if mod(opts.sgWindowFrames, 2) == 0
    opts.sgWindowFrames = opts.sgWindowFrames + 1;
end
opts = default_field(opts, 'sgPolyOrder', 3);
opts.sgPolyOrder = max(1, min(round(opts.sgPolyOrder), opts.sgWindowFrames - 1));
opts = default_field(opts, 'triggerWindowFrames', 5);
opts.triggerWindowFrames = max(1, round(opts.triggerWindowFrames));
opts = default_field(opts, 'minTriggerSamples', opts.triggerWindowFrames);
opts = default_field(opts, 'minTrackFrames', max(opts.sgWindowFrames + 2, opts.triggerWindowFrames + opts.sgWindowFrames));
opts = default_field(opts, 'requireConsecutiveFrames', true);
opts = default_field(opts, 'maxAllowedFrameGap', 1);
opts = default_field(opts, 'excludeSmoothingEdgeFrames', true);
opts = default_field(opts, 'trackPopulation', "strictPrimary");
opts = default_field(opts, 'bulkDirection', "left_to_right");
opts = default_field(opts, 'randomSeed', 42);
opts = default_field(opts, 'throatHeight_mm', 10);
opts = default_field(opts, 'imageSize_px', [1280 320]);
opts = default_field(opts, 'fallbackURef_m_s', 13.32);
opts = default_field(opts, 'fallbackDiameter_m', 100e-6);
opts = default_field(opts, 'makeSanityPlots', true);
opts = default_field(opts, 'nSanityTracks', 5);
opts = default_field(opts, 'makeStationarySanityCheck', true);
opts = default_field(opts, 'stationaryMaxNetDisplacement_px', 2);
opts = default_field(opts, 'stationaryMaxPathLength_px', 5);
end


% =========================================================================
function s = default_field(s, name, value)
if ~isfield(s, name) || isempty(s.(name))
    s.(name) = value;
end
end


% =========================================================================
function lagData = make_case_result_template()
lagData = struct( ...
    'caseName', "", ...
    'Re', NaN, ...
    'kD', NaN, ...
    'pixelSize_mm', NaN, ...
    'dMean_m', NaN, ...
    'U_ref_m_s', NaN, ...
    'sampleAstar', nan(0,1), ...
    'sampleRawAstar', nan(0,1), ...
    'sampleXNorm', nan(0,1), ...
    'sampleYNorm', nan(0,1), ...
    'sampleIsActivatedTrack', false(0,1), ...
    'activatedAllAstar', nan(0,1), ...
    'nonActivatedAllAstar', nan(0,1), ...
    'activatedTriggerAstar', nan(0,1), ...
    'nonActivatedRandomAstar', nan(0,1), ...
    'peakTriggerAstar', nan(0,1), ...
    'growthRatio', nan(0,1), ...
    'activationXY_mm', zeros(0,2), ...
    'activationXYNorm', zeros(0,2), ...
    'stationaryAstar', nan(0,1), ...
    'nStationaryTracks', 0, ...
    'tracks', repmat(make_track_result_template(), 0, 1), ...
    'sanityTracks', repmat(make_sanity_track_template(), 0, 1), ...
    'skipCounts', make_skip_counts(), ...
    'summary', struct(), ...
    'opts', struct());
end


% =========================================================================
function tr = make_track_result_template()
tr = struct( ...
    'TRACK_ID', NaN, ...
    'isActivated', false, ...
    'isStrictPrimary', false, ...
    'activationIndex', NaN, ...
    'activationFrame', NaN, ...
    'activationXY_mm', [NaN NaN], ...
    'activationXYNorm', [NaN NaN], ...
    'nPoints', 0, ...
    'nAccelSamples', 0, ...
    'nTriggerSamples', 0, ...
    'nRandomWindowSamples', 0, ...
    'growthRatio', NaN, ...
    'peakTriggerAstar', NaN, ...
    'allAstar', nan(0,1), ...
    'rawAllAstar', nan(0,1), ...
    'triggerAstar', nan(0,1), ...
    'randomWindowAstar', nan(0,1), ...
    'xNorm', nan(0,1), ...
    'yNorm', nan(0,1));
end


% =========================================================================
function tr = make_sanity_track_template()
tr = struct( ...
    'TRACK_ID', NaN, ...
    'isActivated', false, ...
    'xRawNorm', nan(0,1), ...
    'yRawNorm', nan(0,1), ...
    'xSmoothNorm', nan(0,1), ...
    'ySmoothNorm', nan(0,1), ...
    'activationXYNorm', [NaN NaN]);
end


% =========================================================================
function skip = make_skip_counts()
skip = struct('noCatalog', 0, 'tooShort', 0, 'frameGaps', 0, ...
    'nonFinite', 0, 'noValidAcceleration', 0);
end


% =========================================================================
function [dPool_m, uPool_m_s] = collect_normalization_pools(trajectories, catalog, catMap, areaMap, caseDef, opts)
dPool_m = nan(0,1);
uPool_m_s = nan(0,1);

for ti = 1:numel(trajectories)
    [cat, hasCatalog] = lookup_catalog_entry(trajectories(ti), catalog, catMap);
    if ~hasCatalog || ~track_in_population(cat, opts.trackPopulation)
        continue;
    end
    raw = extract_track_arrays(trajectories(ti), areaMap, caseDef.pixelSize);
    if raw.n < opts.minTrackFrames || ~raw.isUsable
        continue;
    end
    if opts.requireConsecutiveFrames && has_excessive_frame_gap(raw.frame, opts.maxAllowedFrameGap)
        continue;
    end

    areaUse = raw.area_px2;
    actIdx = get_scalar_field(cat, 'strictActivationIndex', NaN);
    if isfinite(actIdx) && actIdx > 1
        areaUse = areaUse(1:min(round(actIdx)-1, numel(areaUse)));
    end
    d_m = sqrt(4 .* areaUse ./ pi) .* caseDef.pixelSize ./ 1000;
    dPool_m = [dPool_m; d_m(isfinite(d_m) & d_m > 0)]; %#ok<AGROW>

    dx_mm = diff(raw.x_m) * 1000;
    dt = diff(raw.t_s);
    if strcmpi(opts.bulkDirection, 'right_to_left')
        uStep = dx_mm ./ dt ./ 1000;
    else
        uStep = -dx_mm ./ dt ./ 1000;
    end
    uStep = uStep(isfinite(uStep) & uStep > 0);
    if isempty(uStep)
        speedStep = hypot(diff(raw.x_m), diff(raw.y_m)) ./ dt;
        uStep = speedStep(isfinite(speedStep) & speedStep > 0);
    end
    uPool_m_s = [uPool_m_s; uStep(:)]; %#ok<AGROW>
end
end


% =========================================================================
function raw = extract_track_arrays(tr, areaMap, pixelSize_mm)
raw = struct('TRACK_ID', NaN, 'spotIds', nan(0,1), 'frame', nan(0,1), ...
    't_s', nan(0,1), 'x_m', nan(0,1), 'y_m', nan(0,1), 'area_px2', nan(0,1), ...
    'n', 0, 'isUsable', false);

if isfield(tr, 'TRACK_ID'), raw.TRACK_ID = tr.TRACK_ID; end
if isfield(tr, 'spotIds'), raw.spotIds = tr.spotIds(:); end
if isfield(tr, 'frame'), raw.frame = double(tr.frame(:)); end
if isfield(tr, 't'), raw.t_s = double(tr.t(:)); end
if isfield(tr, 'x_phys') && ~isempty(tr.x_phys)
    raw.x_m = double(tr.x_phys(:)) ./ 1000;
elseif isfield(tr, 'x')
    raw.x_m = double(tr.x(:)) .* pixelSize_mm ./ 1000;
end
if isfield(tr, 'y_phys') && ~isempty(tr.y_phys)
    raw.y_m = double(tr.y_phys(:)) ./ 1000;
elseif isfield(tr, 'y')
    raw.y_m = double(tr.y(:)) .* pixelSize_mm ./ 1000;
end

n = min([numel(raw.x_m), numel(raw.y_m), numel(raw.t_s), numel(raw.frame)]);
if n < 1
    raw.n = 0;
    return;
end
raw.x_m = raw.x_m(1:n);
raw.y_m = raw.y_m(1:n);
raw.t_s = raw.t_s(1:n);
raw.frame = raw.frame(1:n);
if numel(raw.spotIds) >= n
    raw.spotIds = raw.spotIds(1:n);
else
    raw.spotIds(end+1:n,1) = NaN;
end
raw.area_px2 = lookup_area_values(raw.spotIds, areaMap);

valid = isfinite(raw.x_m) & isfinite(raw.y_m) & isfinite(raw.t_s) & isfinite(raw.frame);
raw.x_m = raw.x_m(valid);
raw.y_m = raw.y_m(valid);
raw.t_s = raw.t_s(valid);
raw.frame = raw.frame(valid);
raw.area_px2 = raw.area_px2(valid);
raw.spotIds = raw.spotIds(valid);

[raw.t_s, ord] = sort(raw.t_s);
raw.x_m = raw.x_m(ord);
raw.y_m = raw.y_m(ord);
raw.frame = raw.frame(ord);
raw.area_px2 = raw.area_px2(ord);
raw.spotIds = raw.spotIds(ord);
raw.n = numel(raw.t_s);
raw.isUsable = raw.n >= 3 && all(isfinite([raw.x_m; raw.y_m; raw.t_s])) && all(diff(raw.t_s) > 0);
end


% =========================================================================
function areaVals = lookup_area_values(spotIds, areaMap)
areaVals = nan(numel(spotIds), 1);
if isempty(spotIds) || isempty(areaMap)
    return;
end
for i = 1:numel(spotIds)
    sid = spotIds(i);
    if isfinite(sid) && isKey(areaMap, sid)
        areaVals(i) = areaMap(sid);
    end
end
end


% =========================================================================
function areaMap = build_spot_area_map(out)
areaMap = [];
if ~isfield(out, 'spots') || ~istable(out.spots) || height(out.spots) == 0 || ...
        ~ismember('ID', out.spots.Properties.VariableNames) || ~ismember('AREA', out.spots.Properties.VariableNames)
    return;
end
areaMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
for i = 1:height(out.spots)
    sid = out.spots.ID(i);
    if isfinite(sid)
        areaMap(sid) = out.spots.AREA(i);
    end
end
end


% =========================================================================
function catMap = build_catalog_map(catalog)
catMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
for i = 1:numel(catalog)
    tid = get_scalar_field(catalog(i), 'TRACK_ID', NaN);
    if isfinite(tid) && ~isKey(catMap, tid)
        catMap(tid) = i;
    end
end
end


% =========================================================================
function [cat, ok] = lookup_catalog_entry(tr, catalog, catMap)
cat = [];
ok = false;
if ~isfield(tr, 'TRACK_ID') || ~isfinite(tr.TRACK_ID) || isempty(catMap) || ~isKey(catMap, tr.TRACK_ID)
    return;
end
cat = catalog(catMap(tr.TRACK_ID));
ok = true;
end


% =========================================================================
function tf = track_in_population(cat, population)
population = lower(char(string(population)));
switch population
    case 'strictprimary'
        tf = get_bool_field(cat, 'isStrictPrimary');
    case 'basicvalid'
        tf = get_bool_field(cat, 'isBasicValid');
    case 'leftmoving'
        tf = get_bool_field(cat, 'isLeftMoving');
    otherwise
        tf = get_bool_field(cat, 'isStrictPrimary');
end
end


% =========================================================================
function tf = has_excessive_frame_gap(frame, maxAllowedGap)
frame = frame(:);
frame = frame(isfinite(frame));
if numel(frame) < 2
    tf = true;
    return;
end
df = diff(round(frame));
tf = any(df < 1 | df > maxAllowedGap);
end


% =========================================================================
function [xSmooth, velocity, acceleration, valid] = local_poly_kinematics(t, x, windowFrames, polyOrder)
t = t(:);
x = x(:);
n = numel(x);
xSmooth = nan(n,1);
velocity = nan(n,1);
acceleration = nan(n,1);
valid = false(n,1);
if n < polyOrder + 1
    return;
end

windowFrames = min(windowFrames, n);
if mod(windowFrames, 2) == 0
    windowFrames = windowFrames - 1;
end
windowFrames = max(windowFrames, polyOrder + 1);
halfWindow = floor(windowFrames / 2);

for i = 1:n
    i0 = max(1, i - halfWindow);
    i1 = min(n, i + halfWindow);
    while (i1 - i0 + 1) < windowFrames
        if i0 > 1
            i0 = i0 - 1;
        elseif i1 < n
            i1 = i1 + 1;
        else
            break;
        end
    end

    idx = i0:i1;
    good = isfinite(t(idx)) & isfinite(x(idx));
    idx = idx(good);
    if numel(idx) < polyOrder + 1
        continue;
    end

    tt = t(idx) - t(i);
    if numel(unique(tt)) < polyOrder + 1
        continue;
    end
    degree = min(polyOrder, numel(idx) - 1);
    p = polyfit(tt, x(idx), degree);
    xSmooth(i) = polyval(p, 0);
    if degree >= 1
        velocity(i) = p(end-1);
    end
    if degree >= 2
        acceleration(i) = 2 * p(end-2);
    end
    valid(i) = isfinite(xSmooth(i)) && isfinite(velocity(i)) && isfinite(acceleration(i));
end
end


% =========================================================================
function [velocity, acceleration] = finite_difference_kinematics(t, x)
t = t(:);
x = x(:);
n = numel(x);
velocity = nan(n,1);
acceleration = nan(n,1);
if n < 3
    return;
end
velocity = derivative_nonuniform(t, x);
acceleration = derivative_nonuniform(t, velocity);
end


% =========================================================================
function dxdt = derivative_nonuniform(t, x)
n = numel(x);
dxdt = nan(n,1);
if n < 2
    return;
end
for i = 1:n
    if i == 1
        i0 = 1; i1 = 2;
    elseif i == n
        i0 = n - 1; i1 = n;
    else
        i0 = i - 1; i1 = i + 1;
    end
    dt = t(i1) - t(i0);
    if isfinite(dt) && dt > 0 && isfinite(x(i0)) && isfinite(x(i1))
        dxdt(i) = (x(i1) - x(i0)) / dt;
    end
end
end


% =========================================================================
function [idx, ok] = trigger_indices_for_activation(actIdx, frame, validAccel, opts)
idx = nan(0,1);
ok = false;
if ~isfinite(actIdx)
    return;
end
actIdx = round(actIdx);
i0 = actIdx - opts.triggerWindowFrames;
i1 = actIdx - 1;
if i0 < 1 || i1 > numel(validAccel)
    return;
end
idx = (i0:i1).';
if any(~validAccel(idx)) || has_excessive_frame_gap(frame(idx), opts.maxAllowedFrameGap)
    idx = idx(validAccel(idx));
end
ok = numel(idx) >= opts.minTriggerSamples;
end


% =========================================================================
function [idx, ok] = random_window_indices(frame, validAccel, opts)
idx = nan(0,1);
ok = false;
n = numel(validAccel);
if n < opts.triggerWindowFrames
    return;
end

starts = [];
for s = 1:(n - opts.triggerWindowFrames + 1)
    cand = s:(s + opts.triggerWindowFrames - 1);
    if all(validAccel(cand)) && ~has_excessive_frame_gap(frame(cand), opts.maxAllowedFrameGap)
        starts(end+1,1) = s; %#ok<AGROW>
    end
end
if isempty(starts)
    return;
end
s = starts(randi(numel(starts)));
idx = (s:(s + opts.triggerWindowFrames - 1)).';
ok = true;
end


% =========================================================================
function ratio = compute_growth_ratio(areaVals)
ratio = NaN;
areaVals = areaVals(:);
areaVals = areaVals(isfinite(areaVals) & areaVals > 0);
if isempty(areaVals)
    return;
end
birthArea = areaVals(1);
maxArea = max(areaVals);
if isfinite(birthArea) && birthArea > 0 && isfinite(maxArea) && maxArea >= birthArea
    ratio = sqrt(maxArea / birthArea);
end
end


% =========================================================================
function tf = is_stationary_candidate(cat, raw, opts, pixelSize_mm)
tf = false;
if get_bool_field(cat, 'isActivated') || get_bool_field(cat, 'isStrictActivated')
    return;
end
if ~get_bool_field(cat, 'isBasicValid')
    return;
end
if raw.n < opts.minTrackFrames
    return;
end
netDisp_m = hypot(raw.x_m(end) - raw.x_m(1), raw.y_m(end) - raw.y_m(1));
pathSteps_m = hypot(diff(raw.x_m), diff(raw.y_m));
pathLength_m = sum(pathSteps_m(isfinite(pathSteps_m)));
netLimit_m = opts.stationaryMaxNetDisplacement_px * pixelSize_mm / 1000;
pathLimit_m = opts.stationaryMaxPathLength_px * pixelSize_mm / 1000;
tf = isfinite(netDisp_m) && isfinite(pathLength_m) && ...
    netDisp_m <= netLimit_m && pathLength_m <= pathLimit_m;
end


% =========================================================================
function sanity = make_sanity_track(tr, raw, smoothX_m, smoothY_m, trackResult, yExtent_mm, opts)
sanity = make_sanity_track_template();
if isfield(tr, 'TRACK_ID'), sanity.TRACK_ID = tr.TRACK_ID; end
sanity.isActivated = trackResult.isActivated;
sanity.xRawNorm = raw.x_m * 1000 ./ opts.throatHeight_mm;
sanity.yRawNorm = (yExtent_mm - raw.y_m * 1000) ./ opts.throatHeight_mm;
sanity.xSmoothNorm = smoothX_m * 1000 ./ opts.throatHeight_mm;
sanity.ySmoothNorm = (yExtent_mm - smoothY_m * 1000) ./ opts.throatHeight_mm;
sanity.activationXYNorm = trackResult.activationXYNorm;
end


% =========================================================================
function sanityTracks = top_up_sanity_tracks(sanityTracks, tracks, trajectories, catalog, catMap, areaMap, caseDef, opts, yExtent_mm)
if isempty(tracks)
    return;
end
needed = opts.nSanityTracks - numel(sanityTracks);
if needed <= 0
    return;
end
trackIds = [tracks.TRACK_ID];
trackIds = trackIds(isfinite(trackIds));
if isempty(trackIds)
    return;
end
trackIds = trackIds(randperm(numel(trackIds), min(numel(trackIds), needed * 3)));

for ii = 1:numel(trackIds)
    if numel(sanityTracks) >= opts.nSanityTracks
        break;
    end
    tid = trackIds(ii);
    trajIdx = find(arrayfun(@(tr) isfield(tr, 'TRACK_ID') && tr.TRACK_ID == tid, trajectories), 1, 'first');
    trkIdx = find([tracks.TRACK_ID] == tid, 1, 'first');
    if isempty(trajIdx) || isempty(trkIdx)
        continue;
    end
    [cat, hasCatalog] = lookup_catalog_entry(trajectories(trajIdx), catalog, catMap);
    if ~hasCatalog || ~track_in_population(cat, opts.trackPopulation)
        continue;
    end
    raw = extract_track_arrays(trajectories(trajIdx), areaMap, caseDef.pixelSize);
    if raw.n < opts.minTrackFrames || ~raw.isUsable
        continue;
    end
    [smoothX_m, ~, ~, ~] = local_poly_kinematics(raw.t_s, raw.x_m, opts.sgWindowFrames, opts.sgPolyOrder);
    [smoothY_m, ~, ~, ~] = local_poly_kinematics(raw.t_s, raw.y_m, opts.sgWindowFrames, opts.sgPolyOrder);
    sanityTracks(end+1,1) = make_sanity_track(trajectories(trajIdx), raw, smoothX_m, smoothY_m, tracks(trkIdx), yExtent_mm, opts); %#ok<AGROW>
end
end


% =========================================================================
function summary = summarize_case(lagData)
summary = struct();
summary.Case = lagData.caseName;
summary.Re = lagData.Re;
summary.kD = lagData.kD;
summary.dMean_um = lagData.dMean_m * 1e6;
summary.U_ref_m_s = lagData.U_ref_m_s;
summary.nTracksUsable = numel(lagData.tracks);
summary.nActivatedTracksUsable = sum(arrayfun(@(tr) tr.isActivated, lagData.tracks));
summary.nNonActivatedTracksUsable = summary.nTracksUsable - summary.nActivatedTracksUsable;
summary.nAllAccelSamples = numel(lagData.sampleAstar);
summary.nActivatedAllSamples = numel(lagData.activatedAllAstar);
summary.nNonActivatedAllSamples = numel(lagData.nonActivatedAllAstar);
summary.nActivatedTriggerSamples = numel(lagData.activatedTriggerAstar);
summary.nNonActivatedRandomSamples = numel(lagData.nonActivatedRandomAstar);
summary.nPeakGrowthPairs = numel(lagData.peakTriggerAstar);
summary.nStationaryTracks = lagData.nStationaryTracks;
summary.nStationarySamples = numel(lagData.stationaryAstar);
summary.activatedAllMedianAstar = finite_median(lagData.activatedAllAstar);
summary.activatedAllP90Astar = finite_prctile(lagData.activatedAllAstar, 90);
summary.nonActivatedAllMedianAstar = finite_median(lagData.nonActivatedAllAstar);
summary.nonActivatedAllP90Astar = finite_prctile(lagData.nonActivatedAllAstar, 90);
summary.activatedTriggerMedianAstar = finite_median(lagData.activatedTriggerAstar);
summary.activatedTriggerP90Astar = finite_prctile(lagData.activatedTriggerAstar, 90);
summary.nonActivatedRandomMedianAstar = finite_median(lagData.nonActivatedRandomAstar);
summary.nonActivatedRandomP90Astar = finite_prctile(lagData.nonActivatedRandomAstar, 90);
summary.stationaryP95Astar = finite_prctile(lagData.stationaryAstar, 95);
summary.sgWindowFrames = lagData.opts.sgWindowFrames;
summary.sgPolyOrder = lagData.opts.sgPolyOrder;
summary.triggerWindowFrames = lagData.opts.triggerWindowFrames;
summary.trackPopulation = string(lagData.opts.trackPopulation);
end


% =========================================================================
function val = get_scalar_field(s, fieldName, defaultVal)
val = defaultVal;
if isstruct(s) && isfield(s, fieldName) && ~isempty(s.(fieldName))
    raw = s.(fieldName);
    raw = raw(:);
    raw = raw(isfinite(raw));
    if ~isempty(raw)
        val = raw(1);
    end
end
end


% =========================================================================
function tf = get_bool_field(s, fieldName)
tf = false;
if ~isstruct(s) || ~isfield(s, fieldName) || isempty(s.(fieldName))
    return;
end
v = s.(fieldName);
if islogical(v)
    tf = any(v(:));
elseif isnumeric(v)
    tf = any(v(:) ~= 0);
end
end


% =========================================================================
function m = finite_mean(x)
x = x(:);
x = x(isfinite(x));
if isempty(x)
    m = NaN;
else
    m = mean(x);
end
end


% =========================================================================
function m = finite_median(x)
x = x(:);
x = x(isfinite(x));
if isempty(x)
    m = NaN;
else
    m = median(x);
end
end


% =========================================================================
function q = finite_prctile(x, p)
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
