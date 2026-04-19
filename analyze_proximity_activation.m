function proxData = analyze_proximity_activation(outChunks, xmlFiles, caseDef, qcOpts, flowOpts, activationOpts, collapseOpts, opts)
%ANALYZE_PROXIMITY_ACTIVATION  Neighbor response around cavitation-collapse events.
%
%   The analysis is intentionally run per XML sample.  The XML files inside
%   one case are independent random video segments, so neighbors are only
%   searched inside the same parsed TrackMate output.

if nargin < 2 || isempty(xmlFiles)
    xmlFiles = {};
end
if nargin < 8 || isempty(opts)
    opts = struct();
end
opts = apply_defaults(opts);

if ~iscell(outChunks)
    outChunks = {outChunks};
end
if isempty(xmlFiles)
    xmlFiles = repmat({''}, numel(outChunks), 1);
else
    xmlFiles = cellstr(string(xmlFiles(:)));
end
if numel(xmlFiles) < numel(outChunks)
    xmlFiles(end+1:numel(outChunks), 1) = {''};
end

proxData = make_result_template();
proxData.caseName = string(caseDef.name);
proxData.Re = caseDef.Re;
proxData.kD = caseDef.kD;
proxData.pixelSize_mm = caseDef.pixelSize;
proxData.opts = opts;

pairRows = {};
eventRows = {};
summarySkip = make_skip_counts();

rng(opts.randomSeed + round(caseDef.Re) + round(caseDef.kD * 1e6), 'twister');
for si = 1:numel(outChunks)
    out = outChunks{si};
    if isempty(out) || ~isstruct(out) || ~isfield(out, 'trajectories') || isempty(out.trajectories)
        summarySkip.emptySamples = summarySkip.emptySamples + 1;
        continue;
    end

    activationOpts.pixelSize = caseDef.pixelSize;
    metrics = trackmate_case_metrics(out, qcOpts, flowOpts, activationOpts);
    collapseData = analyze_collapse_events(out, caseDef.pixelSize, caseDef.dt, collapseOpts);
    [samplePairRows, sampleEventRows, sampleSkip] = analyze_one_sample( ...
        out, metrics, collapseData, caseDef, si, string(xmlFiles{si}), opts);
    pairRows = [pairRows; samplePairRows]; %#ok<AGROW>
    eventRows = [eventRows; sampleEventRows]; %#ok<AGROW>
    summarySkip = add_skip_counts(summarySkip, sampleSkip);
end

proxData.pairTable = pair_rows_to_table(pairRows);
proxData.eventTable = event_rows_to_table(eventRows);
proxData.binnedStats = build_binned_stats(proxData.pairTable, opts);
proxData.summaryTable = build_summary_table(proxData, opts);
proxData.skipCounts = summarySkip;

fprintf('  Proximity activation: %d primary event(s), %d neighbor pair(s), %d secondary activation(s)\n', ...
    height_or_zero(proxData.eventTable), height_or_zero(proxData.pairTable), ...
    count_true_table_var(proxData.pairTable, 'neighborSecondaryActivated'));
end


% =========================================================================
function opts = apply_defaults(opts)
opts = default_field(opts, 'maxGamma', 20);
opts = default_field(opts, 'extendedMaxGamma', 40);
opts = default_field(opts, 'gammaBins', [0 2 5 10 20 40]);
opts = default_field(opts, 'similarSizeRatioRange', [0.5 2.0]);
opts = default_field(opts, 'baselineFrames', 5);
opts = default_field(opts, 'postCollapseFrames', 5);
opts = default_field(opts, 'secondaryStartLagFrames', 1);
opts = default_field(opts, 'secondaryPostCollapseFrames', 5);
opts = default_field(opts, 'referenceFrameTolerance', 2);
opts = default_field(opts, 'minNeighborFramesInWindow', 3);
opts = default_field(opts, 'minCorrelationFrames', 5);
opts = default_field(opts, 'smoothingWindowFrames', 5);
opts = default_field(opts, 'smoothingPolyOrder', 2);
opts = default_field(opts, 'maxCrossCorrelationLagFrames', 5);
opts = default_field(opts, 'useOnlyActivatedPrimaryCollapses', true);
opts = default_field(opts, 'requireBasicValidNeighbors', true);
opts = default_field(opts, 'requireNeighborBaselineBeforePeak', true);
opts = default_field(opts, 'excludeAlreadyActivatedNeighbors', true);
opts = default_field(opts, 'randomSeed', 42);

opts.maxGamma = max(0, double(opts.maxGamma));
opts.extendedMaxGamma = max(opts.maxGamma, double(opts.extendedMaxGamma));
opts.gammaBins = unique(double(opts.gammaBins(:))).';
if numel(opts.gammaBins) < 2
    opts.gammaBins = [0 opts.maxGamma opts.extendedMaxGamma];
end
if opts.gammaBins(1) > 0
    opts.gammaBins = [0 opts.gammaBins];
end
if opts.gammaBins(end) < opts.extendedMaxGamma
    opts.gammaBins(end+1) = opts.extendedMaxGamma;
end
opts.similarSizeRatioRange = sort(double(opts.similarSizeRatioRange(1:2)));
opts.baselineFrames = max(1, round(opts.baselineFrames));
opts.postCollapseFrames = max(0, round(opts.postCollapseFrames));
opts.secondaryStartLagFrames = max(1, round(opts.secondaryStartLagFrames));
opts.secondaryPostCollapseFrames = max(0, round(opts.secondaryPostCollapseFrames));
opts.referenceFrameTolerance = max(0, round(opts.referenceFrameTolerance));
opts.minNeighborFramesInWindow = max(1, round(opts.minNeighborFramesInWindow));
opts.minCorrelationFrames = max(3, round(opts.minCorrelationFrames));
opts.smoothingWindowFrames = max(3, round(opts.smoothingWindowFrames));
if mod(opts.smoothingWindowFrames, 2) == 0
    opts.smoothingWindowFrames = opts.smoothingWindowFrames + 1;
end
opts.smoothingPolyOrder = max(1, min(round(opts.smoothingPolyOrder), opts.smoothingWindowFrames - 1));
opts.maxCrossCorrelationLagFrames = max(0, round(opts.maxCrossCorrelationLagFrames));
opts.useOnlyActivatedPrimaryCollapses = logical(opts.useOnlyActivatedPrimaryCollapses);
opts.requireBasicValidNeighbors = logical(opts.requireBasicValidNeighbors);
opts.requireNeighborBaselineBeforePeak = logical(opts.requireNeighborBaselineBeforePeak);
opts.excludeAlreadyActivatedNeighbors = logical(opts.excludeAlreadyActivatedNeighbors);
end


% =========================================================================
function s = default_field(s, name, value)
if ~isfield(s, name) || isempty(s.(name))
    s.(name) = value;
end
end


% =========================================================================
function proxData = make_result_template()
proxData = struct();
proxData.caseName = "";
proxData.Re = NaN;
proxData.kD = NaN;
proxData.pixelSize_mm = NaN;
proxData.pairTable = table();
proxData.eventTable = table();
proxData.binnedStats = table();
proxData.summaryTable = table();
proxData.skipCounts = make_skip_counts();
proxData.opts = struct();
end


% =========================================================================
function skip = make_skip_counts()
skip = struct( ...
    'emptySamples', 0, ...
    'collapseEventsTotal', 0, ...
    'primaryNoTrack', 0, ...
    'primaryNotActivated', 0, ...
    'primaryBadRadius', 0, ...
    'primaryBadTiming', 0, ...
    'neighborNoOverlap', 0, ...
    'neighborAlreadyActivated', 0, ...
    'neighborNoBaseline', 0, ...
    'neighborBadSize', 0, ...
    'neighborOutsideGamma', 0, ...
    'neighborTooFewResponseFrames', 0);
end


% =========================================================================
function a = add_skip_counts(a, b)
names = fieldnames(a);
for i = 1:numel(names)
    if isfield(b, names{i})
        a.(names{i}) = a.(names{i}) + b.(names{i});
    end
end
end


% =========================================================================
function [pairRows, eventRows, skip] = analyze_one_sample(out, metrics, collapseData, caseDef, sampleIndex, xmlFile, opts)
pairRows = {};
eventRows = {};
skip = make_skip_counts();

trackRecords = build_track_records(out, metrics);
if isempty(trackRecords)
    return;
end

[collapseFrame, collapseX, collapseY, collapseTrackId] = collapse_arrays(collapseData);
skip.collapseEventsTotal = numel(collapseFrame);
if isempty(collapseFrame)
    return;
end

for ci = 1:numel(collapseFrame)
    primaryIdx = find_track_record(trackRecords, collapseTrackId(ci));
    if ~isfinite(primaryIdx)
        skip.primaryNoTrack = skip.primaryNoTrack + 1;
        continue;
    end
    primary = trackRecords(primaryIdx);
    if opts.useOnlyActivatedPrimaryCollapses && ~primary.isActivated
        skip.primaryNotActivated = skip.primaryNotActivated + 1;
        continue;
    end

    primaryEvent = characterize_primary(primary, collapseFrame(ci), collapseX(ci), collapseY(ci), opts);
    if ~primaryEvent.hasGoodRadius
        skip.primaryBadRadius = skip.primaryBadRadius + 1;
        continue;
    end
    if ~primaryEvent.hasGoodTiming
        skip.primaryBadTiming = skip.primaryBadTiming + 1;
        continue;
    end

    primaryPairRows = {};
    primaryEventStats = make_primary_event_stats();
    eventStart = primaryEvent.activationFrame;
    if ~isfinite(eventStart)
        eventStart = primaryEvent.peakFrame;
    end
    eventStop = primaryEvent.collapseFrame;
    secondaryStart = primaryEvent.activationFrame + opts.secondaryStartLagFrames;
    secondaryStop = primaryEvent.collapseFrame + opts.secondaryPostCollapseFrames;
    responseStart = primaryEvent.peakFrame;
    responseStop = primaryEvent.collapseFrame + opts.postCollapseFrames;

    for ni = 1:numel(trackRecords)
        if ni == primaryIdx
            continue;
        end
        neighbor = trackRecords(ni);
        if opts.requireBasicValidNeighbors && ~neighbor.isBasicValid
            continue;
        end
        if opts.excludeAlreadyActivatedNeighbors && neighbor.isActivated && ...
                isfinite(neighbor.activationFrame) && neighbor.activationFrame < secondaryStart
            skip.neighborAlreadyActivated = skip.neighborAlreadyActivated + 1;
            continue;
        end

        if count_frames_in_window(neighbor.frame, eventStart, eventStop) < 1
            skip.neighborNoOverlap = skip.neighborNoOverlap + 1;
            continue;
        end

        neighborR0 = baseline_radius(neighbor, primaryEvent.peakFrame, opts.baselineFrames);
        if opts.requireNeighborBaselineBeforePeak && ~(isfinite(neighborR0) && neighborR0 > 0)
            skip.neighborNoBaseline = skip.neighborNoBaseline + 1;
            continue;
        end
        if ~(isfinite(neighborR0) && neighborR0 > 0)
            neighborR0 = radius_near_frame(neighbor, primaryEvent.peakFrame, opts.referenceFrameTolerance);
        end
        if ~(isfinite(neighborR0) && neighborR0 > 0)
            neighborR0 = first_finite_positive(neighbor.radius_mm);
        end
        sizeRatio = neighborR0 / primaryEvent.primaryR0_mm;
        if ~(isfinite(sizeRatio) && sizeRatio >= opts.similarSizeRatioRange(1) && ...
                sizeRatio <= opts.similarSizeRatioRange(2))
            skip.neighborBadSize = skip.neighborBadSize + 1;
            continue;
        end

        [minDist, nOverlap] = min_track_distance_in_window(primary, neighbor, eventStart, eventStop);
        gammaMin = minDist / primaryEvent.primaryRmax_mm;
        if ~(isfinite(gammaMin) && gammaMin <= opts.extendedMaxGamma)
            skip.neighborOutsideGamma = skip.neighborOutsideGamma + 1;
            continue;
        end

        [distAct, gammaAct] = reference_distance(primary, neighbor, primaryEvent.activationFrame, ...
            primaryEvent.primaryRmax_mm, opts.referenceFrameTolerance);
        [distPeak, gammaPeak] = reference_distance(primary, neighbor, primaryEvent.peakFrame, ...
            primaryEvent.primaryRmax_mm, opts.referenceFrameTolerance);
        [distCollapse, gammaCollapse] = reference_distance(primary, neighbor, primaryEvent.collapseFrame, ...
            primaryEvent.primaryRmax_mm, opts.referenceFrameTolerance);
        gammaForPlot = gammaPeak;
        if ~isfinite(gammaForPlot)
            gammaForPlot = gammaMin;
        end

        [maxAbsResp, maxExpansion, maxCompression, endResp, nResponse] = ...
            radius_response_metrics(neighbor, responseStart, responseStop, neighborR0);
        if nResponse < opts.minNeighborFramesInWindow
            skip.neighborTooFewResponseFrames = skip.neighborTooFewResponseFrames + 1;
            continue;
        end
        randomResp = random_window_response(neighbor, responseStart, responseStop, opts);

        [isSecondary, afterPeak, afterCollapse, lagA, lagP, lagC] = ...
            secondary_activation_flags(neighbor, primaryEvent, secondaryStart, secondaryStop);
        [rhoZero, rhoMax, lagAtMax, rdotRatio, nCorr] = ...
            radius_rate_coupling(primary, neighbor, eventStart, responseStop, opts);

        row = { ...
            string(caseDef.name), caseDef.Re, caseDef.kD, sampleIndex, xmlFile, ...
            primary.TRACK_ID, primaryEvent.activationFrame, primaryEvent.peakFrame, primaryEvent.collapseFrame, ...
            primaryEvent.primaryR0_mm, primaryEvent.primaryRmax_mm, primaryEvent.peakArea_px2, ...
            neighbor.TRACK_ID, neighborR0, sizeRatio, ...
            distAct, distPeak, distCollapse, minDist, ...
            gammaAct, gammaPeak, gammaCollapse, gammaMin, gammaForPlot, ...
            maxAbsResp, maxExpansion, maxCompression, endResp, randomResp, ...
            neighbor.isActivated, neighbor.activationFrame, isSecondary, afterPeak, afterCollapse, ...
            lagA, lagP, lagC, ...
            rhoZero, rhoMax, lagAtMax, rdotRatio, nCorr, nResponse};
        primaryPairRows(end+1, :) = row; %#ok<AGROW>

        primaryEventStats = update_primary_event_stats(primaryEventStats, gammaMin, isSecondary);
    end

    pairRows = [pairRows; primaryPairRows]; %#ok<AGROW>
    eventRows(end+1, :) = primary_event_row(caseDef, sampleIndex, xmlFile, primary, primaryEvent, primaryEventStats); %#ok<AGROW>
end
end


% =========================================================================
function records = build_track_records(out, metrics)
records = make_track_record_template();
records = records([]);
if isempty(out) || ~isstruct(out) || ~isfield(out, 'trajectories') || isempty(out.trajectories)
    return;
end
areaMap = build_spot_area_map(out);
catalog = [];
if isstruct(metrics) && isfield(metrics, 'trackCatalog')
    catalog = metrics.trackCatalog(:);
end

traj = out.trajectories(:);
records = repmat(make_track_record_template(), numel(traj), 1);
for i = 1:numel(traj)
    tr = traj(i);
    rec = make_track_record_template();
    rec.TRACK_ID = get_scalar_field(tr, 'TRACK_ID', NaN);
    if isfield(tr, 'frame'), rec.frame = double(tr.frame(:)); end
    if isfield(tr, 't'), rec.t_s = double(tr.t(:)); end
    if isfield(tr, 'x_phys'), rec.x_mm = double(tr.x_phys(:)); end
    if isfield(tr, 'y_phys'), rec.y_mm = double(tr.y_phys(:)); end
    if isfield(tr, 'spotIds'), rec.spotIds = double(tr.spotIds(:)); end
    rec.area_px2 = area_values_for_spots(rec.spotIds, areaMap);
    rec.radius_mm = sqrt(max(rec.area_px2, 0) ./ pi) .* get_pixel_size_from_metrics(metrics);
    if numel(rec.radius_mm) ~= numel(rec.frame)
        rec.radius_mm = align_vector_length(rec.radius_mm, numel(rec.frame));
    end

    catIdx = find_catalog_record(catalog, rec.TRACK_ID);
    if isfinite(catIdx)
        cat = catalog(catIdx);
        rec.isBasicValid = get_bool_field(cat, 'isBasicValid');
        rec.isActivated = get_bool_field(cat, 'isActivated');
        rec.isStrictActivated = get_bool_field(cat, 'isStrictActivated');
        rec.activationFrame = get_scalar_field(cat, 'activationFrame', NaN);
        rec.activationIndex = get_scalar_field(cat, 'activationIndex', NaN);
        if ~isfinite(rec.activationFrame)
            rec.activationFrame = get_scalar_field(cat, 'strictActivationFrame', NaN);
        end
        if ~isfinite(rec.activationIndex)
            rec.activationIndex = get_scalar_field(cat, 'strictActivationIndex', NaN);
        end
    end
    if ~isfinite(rec.activationFrame) && isfinite(rec.activationIndex)
        idx = round(rec.activationIndex);
        if idx >= 1 && idx <= numel(rec.frame)
            rec.activationFrame = rec.frame(idx);
        end
    end
    records(i) = rec;
end
end


% =========================================================================
function rec = make_track_record_template()
rec = struct( ...
    'TRACK_ID', NaN, ...
    'frame', nan(0,1), ...
    't_s', nan(0,1), ...
    'x_mm', nan(0,1), ...
    'y_mm', nan(0,1), ...
    'spotIds', nan(0,1), ...
    'area_px2', nan(0,1), ...
    'radius_mm', nan(0,1), ...
    'isBasicValid', false, ...
    'isActivated', false, ...
    'isStrictActivated', false, ...
    'activationFrame', NaN, ...
    'activationIndex', NaN);
end


% =========================================================================
function px = get_pixel_size_from_metrics(metrics)
px = NaN;
if isstruct(metrics) && isfield(metrics, 'activationOpts') && isstruct(metrics.activationOpts) && ...
        isfield(metrics.activationOpts, 'pixelSize')
    px = metrics.activationOpts.pixelSize;
end
if ~(isfinite(px) && px > 0)
    px = 1;
end
end


% =========================================================================
function areaMap = build_spot_area_map(out)
areaMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
if ~isstruct(out) || ~isfield(out, 'spots') || ~istable(out.spots) || ...
        ~ismember('ID', out.spots.Properties.VariableNames) || ...
        ~ismember('AREA', out.spots.Properties.VariableNames)
    return;
end
spots = out.spots;
for i = 1:height(spots)
    sid = double(spots.ID(i));
    if isfinite(sid)
        areaMap(sid) = double(spots.AREA(i));
    end
end
end


% =========================================================================
function areaVals = area_values_for_spots(spotIds, areaMap)
areaVals = nan(numel(spotIds), 1);
for i = 1:numel(spotIds)
    sid = double(spotIds(i));
    if isfinite(sid) && isKey(areaMap, sid)
        areaVals(i) = areaMap(sid);
    end
end
end


% =========================================================================
function v = align_vector_length(v, n)
v = v(:);
if numel(v) >= n
    v = v(1:n);
else
    v(end+1:n, 1) = NaN;
end
end


% =========================================================================
function idx = find_catalog_record(catalog, trackId)
idx = NaN;
if isempty(catalog) || ~isfinite(trackId) || ~isfield(catalog, 'TRACK_ID')
    return;
end
ids = [catalog.TRACK_ID];
j = find(ids == trackId, 1, 'first');
if ~isempty(j)
    idx = j;
end
end


% =========================================================================
function idx = find_track_record(records, trackId)
idx = NaN;
if isempty(records) || ~isfinite(trackId)
    return;
end
ids = [records.TRACK_ID];
j = find(ids == trackId, 1, 'first');
if ~isempty(j)
    idx = j;
end
end


% =========================================================================
function primaryEvent = characterize_primary(primary, collapseFrame, collapseX, collapseY, opts)
primaryEvent = struct();
primaryEvent.activationFrame = primary.activationFrame;
primaryEvent.peakFrame = NaN;
primaryEvent.collapseFrame = collapseFrame;
primaryEvent.primaryR0_mm = NaN;
primaryEvent.primaryRmax_mm = NaN;
primaryEvent.peakArea_px2 = NaN;
primaryEvent.activationX_mm = NaN;
primaryEvent.activationY_mm = NaN;
primaryEvent.peakX_mm = NaN;
primaryEvent.peakY_mm = NaN;
primaryEvent.collapseX_mm = collapseX;
primaryEvent.collapseY_mm = collapseY;
primaryEvent.hasGoodRadius = false;
primaryEvent.hasGoodTiming = false;

radius = primary.radius_mm(:);
areaVals = primary.area_px2(:);
validR = isfinite(radius) & radius > 0;
if ~any(validR)
    return;
end
peakIdx = find(radius == max(radius(validR)), 1, 'first');
if isempty(peakIdx) || peakIdx < 1 || peakIdx > numel(primary.frame)
    return;
end
primaryEvent.primaryRmax_mm = radius(peakIdx);
primaryEvent.peakFrame = primary.frame(peakIdx);
primaryEvent.peakX_mm = primary.x_mm(peakIdx);
primaryEvent.peakY_mm = primary.y_mm(peakIdx);
if peakIdx <= numel(areaVals)
    primaryEvent.peakArea_px2 = areaVals(peakIdx);
end

if isfinite(primary.activationFrame)
    primaryEvent.primaryR0_mm = baseline_radius(primary, primary.activationFrame, opts.baselineFrames);
    [primaryEvent.activationX_mm, primaryEvent.activationY_mm] = ...
        position_near_frame(primary, primary.activationFrame, opts.referenceFrameTolerance);
end
if ~(isfinite(primaryEvent.primaryR0_mm) && primaryEvent.primaryR0_mm > 0)
    primaryEvent.primaryR0_mm = first_finite_positive(radius);
end

primaryEvent.hasGoodRadius = isfinite(primaryEvent.primaryR0_mm) && primaryEvent.primaryR0_mm > 0 && ...
    isfinite(primaryEvent.primaryRmax_mm) && primaryEvent.primaryRmax_mm > 0;
primaryEvent.hasGoodTiming = isfinite(primaryEvent.peakFrame) && isfinite(primaryEvent.collapseFrame) && ...
    primaryEvent.collapseFrame >= primaryEvent.peakFrame && ...
    (~isfinite(primaryEvent.activationFrame) || primaryEvent.peakFrame >= primaryEvent.activationFrame);
end


% =========================================================================
function r0 = baseline_radius(track, refFrame, nFrames)
r0 = NaN;
if ~(isfinite(refFrame) && nFrames > 0)
    return;
end
frame = track.frame(:);
radius = track.radius_mm(:);
n = min(numel(frame), numel(radius));
frame = frame(1:n);
radius = radius(1:n);
mask = isfinite(frame) & frame >= refFrame - nFrames & frame < refFrame & ...
    isfinite(radius) & radius > 0;
vals = radius(mask);
if isempty(vals)
    return;
end
r0 = median(vals);
end


% =========================================================================
function r = radius_near_frame(track, refFrame, toleranceFrames)
r = NaN;
if ~isfinite(refFrame)
    return;
end
frame = track.frame(:);
radius = track.radius_mm(:);
n = min(numel(frame), numel(radius));
frame = frame(1:n);
radius = radius(1:n);
mask = isfinite(frame) & abs(frame - refFrame) <= toleranceFrames & isfinite(radius) & radius > 0;
if ~any(mask)
    return;
end
[~, j] = min(abs(frame(mask) - refFrame));
vals = radius(mask);
r = vals(j);
end


% =========================================================================
function v = first_finite_positive(x)
x = x(:);
idx = find(isfinite(x) & x > 0, 1, 'first');
if isempty(idx)
    v = NaN;
else
    v = x(idx);
end
end


% =========================================================================
function n = count_frames_in_window(frame, f0, f1)
if ~(isfinite(f0) && isfinite(f1))
    n = 0;
    return;
end
frame = frame(:);
n = sum(isfinite(frame) & frame >= f0 & frame <= f1);
end


% =========================================================================
function [minDist, nOverlap] = min_track_distance_in_window(primary, neighbor, f0, f1)
minDist = NaN;
nOverlap = 0;
if ~(isfinite(f0) && isfinite(f1))
    return;
end
[commonFrames, ip, in] = intersect(primary.frame(:), neighbor.frame(:), 'stable');
mask = isfinite(commonFrames) & commonFrames >= f0 & commonFrames <= f1;
ip = ip(mask);
in = in(mask);
nOverlap = numel(ip);
if nOverlap == 0
    return;
end
dx = primary.x_mm(ip) - neighbor.x_mm(in);
dy = primary.y_mm(ip) - neighbor.y_mm(in);
dist = hypot(dx, dy);
dist = dist(isfinite(dist));
if isempty(dist)
    nOverlap = 0;
    return;
end
minDist = min(dist);
end


% =========================================================================
function [dist, gamma] = reference_distance(primary, neighbor, refFrame, rmax, toleranceFrames)
dist = NaN;
gamma = NaN;
if ~(isfinite(refFrame) && isfinite(rmax) && rmax > 0)
    return;
end
[xp, yp] = position_near_frame(primary, refFrame, toleranceFrames);
[xn, yn] = position_near_frame(neighbor, refFrame, toleranceFrames);
if all(isfinite([xp yp xn yn]))
    dist = hypot(xp - xn, yp - yn);
    gamma = dist / rmax;
end
end


% =========================================================================
function [x, y] = position_near_frame(track, refFrame, toleranceFrames)
x = NaN;
y = NaN;
if ~isfinite(refFrame)
    return;
end
frame = track.frame(:);
n = min([numel(frame), numel(track.x_mm), numel(track.y_mm)]);
frame = frame(1:n);
xv = track.x_mm(1:n);
yv = track.y_mm(1:n);
mask = isfinite(frame) & abs(frame - refFrame) <= toleranceFrames & isfinite(xv) & isfinite(yv);
if ~any(mask)
    return;
end
idxAll = find(mask);
[~, j] = min(abs(frame(idxAll) - refFrame));
idx = idxAll(j);
x = xv(idx);
y = yv(idx);
end


% =========================================================================
function [maxAbsResp, maxExpansion, maxCompression, endResp, nResponse] = ...
    radius_response_metrics(track, f0, f1, r0)
maxAbsResp = NaN;
maxExpansion = NaN;
maxCompression = NaN;
endResp = NaN;
nResponse = 0;
if ~(isfinite(f0) && isfinite(f1) && isfinite(r0) && r0 > 0)
    return;
end
frame = track.frame(:);
radius = track.radius_mm(:);
n = min(numel(frame), numel(radius));
frame = frame(1:n);
radius = radius(1:n);
mask = isfinite(frame) & frame >= f0 & frame <= f1 & isfinite(radius) & radius > 0;
r = radius(mask);
nResponse = numel(r);
if isempty(r)
    return;
end
dr = (r - r0) ./ r0;
maxAbsResp = max(abs(dr));
maxExpansion = max(dr);
maxCompression = min(dr);
endResp = dr(end);
end


% =========================================================================
function randomResp = random_window_response(track, responseStart, responseStop, opts)
randomResp = NaN;
if ~(isfinite(responseStart) && isfinite(responseStop))
    return;
end
winLen = max(1, round(responseStop - responseStart + 1));
frame = track.frame(:);
frame = frame(isfinite(frame));
if numel(frame) < winLen + opts.baselineFrames
    return;
end
fMin = min(frame);
fMax = max(frame);
startCandidates = (ceil(fMin + opts.baselineFrames):floor(fMax - winLen + 1)).';
if isempty(startCandidates)
    return;
end
eventBuffer0 = responseStart - opts.baselineFrames;
eventBuffer1 = responseStop;
startCandidates = startCandidates((startCandidates + winLen - 1) < eventBuffer0 | startCandidates > eventBuffer1);
if isempty(startCandidates)
    return;
end
startCandidates = startCandidates(randperm(numel(startCandidates)));
for i = 1:numel(startCandidates)
    f0 = startCandidates(i);
    r0 = baseline_radius(track, f0, opts.baselineFrames);
    if ~(isfinite(r0) && r0 > 0)
        continue;
    end
    [resp, ~, ~, ~, nResp] = radius_response_metrics(track, f0, f0 + winLen - 1, r0);
    if nResp >= opts.minNeighborFramesInWindow && isfinite(resp)
        randomResp = resp;
        return;
    end
end
end


% =========================================================================
function [isSecondary, afterPeak, afterCollapse, lagA, lagP, lagC] = ...
    secondary_activation_flags(neighbor, primaryEvent, secondaryStart, secondaryStop)
isSecondary = false;
afterPeak = false;
afterCollapse = false;
lagA = NaN;
lagP = NaN;
lagC = NaN;
if ~(neighbor.isActivated && isfinite(neighbor.activationFrame))
    return;
end
af = neighbor.activationFrame;
lagA = af - primaryEvent.activationFrame;
lagP = af - primaryEvent.peakFrame;
lagC = af - primaryEvent.collapseFrame;
isSecondary = isfinite(secondaryStart) && isfinite(secondaryStop) && af >= secondaryStart && af <= secondaryStop;
afterPeak = isSecondary && isfinite(primaryEvent.peakFrame) && af >= primaryEvent.peakFrame;
afterCollapse = isSecondary && isfinite(primaryEvent.collapseFrame) && af >= primaryEvent.collapseFrame;
end


% =========================================================================
function [rhoZero, rhoMax, lagAtMax, rmsRatio, nCommon] = radius_rate_coupling(primary, neighbor, f0, f1, opts)
rhoZero = NaN;
rhoMax = NaN;
lagAtMax = NaN;
rmsRatio = NaN;
nCommon = 0;
if ~(isfinite(f0) && isfinite(f1))
    return;
end
[pRdot, pValid] = smooth_radius_rate(primary, opts);
[nRdot, nValid] = smooth_radius_rate(neighbor, opts);
[commonFrames, ip, in] = intersect(primary.frame(:), neighbor.frame(:), 'stable');
inBounds = ip <= numel(pValid) & in <= numel(nValid) & ip <= numel(pRdot) & in <= numel(nRdot);
mask = inBounds & commonFrames >= f0 & commonFrames <= f1 & pValid(ip) & nValid(in) & ...
    isfinite(pRdot(ip)) & isfinite(nRdot(in));
ip = ip(mask);
in = in(mask);
nCommon = numel(ip);
if nCommon < opts.minCorrelationFrames
    return;
end
p = pRdot(ip);
n = nRdot(in);
rhoZero = pearson_corr(p, n);
rmsP = sqrt(mean(p.^2));
rmsN = sqrt(mean(n.^2));
if isfinite(rmsP) && rmsP > 0
    rmsRatio = rmsN / rmsP;
end

bestAbs = -Inf;
for lag = -opts.maxCrossCorrelationLagFrames:opts.maxCrossCorrelationLagFrames
    if lag < 0
        pp = p((1-lag):end);
        nn = n(1:(end+lag));
    elseif lag > 0
        pp = p(1:(end-lag));
        nn = n((1+lag):end);
    else
        pp = p;
        nn = n;
    end
    if numel(pp) < opts.minCorrelationFrames || numel(nn) < opts.minCorrelationFrames
        continue;
    end
    rho = pearson_corr(pp, nn);
    if isfinite(rho) && abs(rho) > bestAbs
        bestAbs = abs(rho);
        rhoMax = rho;
        lagAtMax = lag;
    end
end
end


% =========================================================================
function [rdot, valid] = smooth_radius_rate(track, opts)
radius = track.radius_mm(:);
frame = track.frame(:);
t = track.t_s(:);
nFull = numel(frame);
rdot = nan(nFull,1);
valid = false(nFull,1);
n = min([numel(radius), nFull, numel(t)]);
radiusUse = radius(1:n);
tUse = t(1:n);
if n < opts.smoothingWindowFrames || any(~isfinite(radiusUse)) || any(~isfinite(tUse))
    return;
end
[~, rdotUse, ~, validUse] = local_poly_kinematics(tUse, radiusUse, opts.smoothingWindowFrames, opts.smoothingPolyOrder);
rdot(1:n) = rdotUse;
valid(1:n) = validUse;
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
if n < windowFrames || numel(t) ~= n
    return;
end
halfWindow = floor(windowFrames / 2);
for i = 1:n
    i0 = max(1, i - halfWindow);
    i1 = min(n, i + halfWindow);
    if (i1 - i0 + 1) < windowFrames
        continue;
    end
    tt = t(i0:i1) - t(i);
    xx = x(i0:i1);
    if any(~isfinite(tt)) || any(~isfinite(xx)) || numel(unique(tt)) <= polyOrder
        continue;
    end
    try
        p = polyfit(tt, xx, polyOrder);
    catch
        continue;
    end
    xSmooth(i) = polyval(p, 0);
    dp = polyder(p);
    velocity(i) = polyval(dp, 0);
    if numel(dp) >= 2
        ddp = polyder(dp);
        acceleration(i) = polyval(ddp, 0);
    else
        acceleration(i) = 0;
    end
    valid(i) = true;
end
end


% =========================================================================
function rho = pearson_corr(x, y)
x = x(:);
y = y(:);
n = min(numel(x), numel(y));
x = x(1:n);
y = y(1:n);
valid = isfinite(x) & isfinite(y);
x = x(valid);
y = y(valid);
if numel(x) < 3
    rho = NaN;
    return;
end
x = x - mean(x);
y = y - mean(y);
den = sqrt(sum(x.^2) * sum(y.^2));
if den <= 0
    rho = NaN;
else
    rho = sum(x .* y) / den;
end
end


% =========================================================================
function stats = make_primary_event_stats()
stats = struct( ...
    'nViableNeighbors', 0, ...
    'nViableWithin5Rmax', 0, ...
    'nViableWithin10Rmax', 0, ...
    'nViableWithin20Rmax', 0, ...
    'nSecondaryActivatedWithin5Rmax', 0, ...
    'nSecondaryActivatedWithin10Rmax', 0, ...
    'nSecondaryActivatedWithin20Rmax', 0, ...
    'nearestViableGamma', NaN, ...
    'nearestStableGamma', NaN, ...
    'nearestActivatedGamma', NaN);
end


% =========================================================================
function stats = update_primary_event_stats(stats, gammaMin, isSecondary)
if ~(isfinite(gammaMin) && gammaMin >= 0)
    return;
end
stats.nViableNeighbors = stats.nViableNeighbors + 1;
if gammaMin <= 5
    stats.nViableWithin5Rmax = stats.nViableWithin5Rmax + 1;
end
if gammaMin <= 10
    stats.nViableWithin10Rmax = stats.nViableWithin10Rmax + 1;
end
if gammaMin <= 20
    stats.nViableWithin20Rmax = stats.nViableWithin20Rmax + 1;
end
if isSecondary && gammaMin <= 5
    stats.nSecondaryActivatedWithin5Rmax = stats.nSecondaryActivatedWithin5Rmax + 1;
end
if isSecondary && gammaMin <= 10
    stats.nSecondaryActivatedWithin10Rmax = stats.nSecondaryActivatedWithin10Rmax + 1;
end
if isSecondary && gammaMin <= 20
    stats.nSecondaryActivatedWithin20Rmax = stats.nSecondaryActivatedWithin20Rmax + 1;
end
stats.nearestViableGamma = finite_min_pair(stats.nearestViableGamma, gammaMin);
if isSecondary
    stats.nearestActivatedGamma = finite_min_pair(stats.nearestActivatedGamma, gammaMin);
else
    stats.nearestStableGamma = finite_min_pair(stats.nearestStableGamma, gammaMin);
end
end


% =========================================================================
function v = finite_min_pair(a, b)
if ~isfinite(a)
    v = b;
elseif ~isfinite(b)
    v = a;
else
    v = min(a, b);
end
end


% =========================================================================
function row = primary_event_row(caseDef, sampleIndex, xmlFile, primary, primaryEvent, stats)
row = { ...
    string(caseDef.name), caseDef.Re, caseDef.kD, sampleIndex, xmlFile, ...
    primary.TRACK_ID, primaryEvent.activationFrame, primaryEvent.peakFrame, primaryEvent.collapseFrame, ...
    primaryEvent.primaryR0_mm, primaryEvent.primaryRmax_mm, primaryEvent.peakArea_px2, ...
    stats.nViableNeighbors, stats.nViableWithin5Rmax, stats.nViableWithin10Rmax, stats.nViableWithin20Rmax, ...
    stats.nSecondaryActivatedWithin5Rmax, stats.nSecondaryActivatedWithin10Rmax, stats.nSecondaryActivatedWithin20Rmax, ...
    stats.nearestViableGamma, stats.nearestStableGamma, stats.nearestActivatedGamma};
end


% =========================================================================
function T = pair_rows_to_table(rows)
if isempty(rows)
    T = table();
    return;
end
T = cell2table(rows, 'VariableNames', { ...
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
end


% =========================================================================
function T = event_rows_to_table(rows)
if isempty(rows)
    T = table();
    return;
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
function T = build_binned_stats(pairTable, opts)
T = table();
if isempty(pairTable) || height(pairTable) == 0
    return;
end
caseNames = string(pairTable.Case);
ReVals = double(pairTable.Re);
kDVals = double(pairTable.kD);
groups = unique([caseNames, string(ReVals), string(kDVals)], 'rows', 'stable');
rows = {};
for gi = 1:size(groups, 1)
    cName = groups(gi, 1);
    Re = str2double(groups(gi, 2));
    kD = str2double(groups(gi, 3));
    gMask = caseNames == cName & ReVals == Re & kDVals == kD;
    for bi = 1:(numel(opts.gammaBins) - 1)
        lo = opts.gammaBins(bi);
        hi = opts.gammaBins(bi + 1);
        gamma = double(pairTable.gammaMinDuringEvent);
        inBin = gMask & gamma >= lo & gamma < hi;
        if bi == (numel(opts.gammaBins) - 1)
            inBin = gMask & gamma >= lo & gamma <= hi;
        end
        nPairs = sum(inBin);
        if nPairs == 0
            continue;
        end
        secondary = logical(pairTable.neighborSecondaryActivated(inBin));
        nSecondary = sum(secondary);
        p = nSecondary / nPairs;
        [ciLow, ciHigh] = wilson_ci_local(nSecondary, nPairs);
        resp = double(pairTable.maxAbsDeltaR_over_R0(inBin));
        randomResp = double(pairTable.randomMaxAbsDeltaR_over_R0(inBin));
        rdotRatio = double(pairTable.rdotRmsRatio(inBin));
        rows(end+1,:) = { ...
            cName, Re, kD, lo, hi, (lo + hi) / 2, nPairs, nSecondary, p, ciLow, ciHigh, ...
            finite_median(resp), finite_prctile(resp, 90), ...
            finite_median(randomResp), finite_prctile(randomResp, 90), ...
            finite_median(rdotRatio), finite_prctile(rdotRatio, 90)}; %#ok<AGROW>
    end
end
if ~isempty(rows)
    T = cell2table(rows, 'VariableNames', { ...
        'Case','Re','kD','gammaBinLow','gammaBinHigh','gammaBinCenter', ...
        'nPairs','nSecondaryActivated','secondaryActivationProbability','secondaryActivationCI_low','secondaryActivationCI_high', ...
        'medianMaxAbsDeltaR_over_R0','p90MaxAbsDeltaR_over_R0', ...
        'medianRandomMaxAbsDeltaR_over_R0','p90RandomMaxAbsDeltaR_over_R0', ...
        'medianRdotRmsRatio','p90RdotRmsRatio'});
end
end


% =========================================================================
function T = build_summary_table(proxData, opts)
pairTable = proxData.pairTable;
eventTable = proxData.eventTable;
nPrimary = height_or_zero(eventTable);
nPairs = height_or_zero(pairTable);
nSecondary = count_true_table_var(pairTable, 'neighborSecondaryActivated');
nPrimaryWithNeighbors = 0;
nWithin10 = 0;
nSecondaryWithin10 = 0;
nearestStableMedian = NaN;
medianResponse = NaN;
medianRandom = NaN;
medianRdotRatio = NaN;

if nPrimary > 0
    nPrimaryWithNeighbors = sum(double(eventTable.nViableNeighbors) > 0);
end
if nPairs > 0
    gamma = double(pairTable.gammaMinDuringEvent);
    sec = logical(pairTable.neighborSecondaryActivated);
    nWithin10 = sum(isfinite(gamma) & gamma <= 10);
    nSecondaryWithin10 = sum(isfinite(gamma) & gamma <= 10 & sec);
    nearestStableMedian = finite_median(double(eventTable.nearestStableGamma));
    medianResponse = finite_median(double(pairTable.maxAbsDeltaR_over_R0));
    medianRandom = finite_median(double(pairTable.randomMaxAbsDeltaR_over_R0));
    medianRdotRatio = finite_median(double(pairTable.rdotRmsRatio));
end
pSecondary = safe_div(nSecondary, nPairs);
pSecondaryWithin10 = safe_div(nSecondaryWithin10, nWithin10);

T = table( ...
    proxData.caseName, proxData.Re, proxData.kD, ...
    nPrimary, nPrimaryWithNeighbors, nPairs, nSecondary, pSecondary, ...
    nWithin10, nSecondaryWithin10, pSecondaryWithin10, ...
    nearestStableMedian, medianResponse, medianRandom, medianRdotRatio, ...
    opts.maxGamma, opts.extendedMaxGamma, opts.secondaryPostCollapseFrames, ...
    'VariableNames', {'Case','Re','kD', ...
    'nPrimaryEvents','nPrimaryEventsWithNeighbors','nNeighborPairs','nSecondaryActivations','secondaryActivationProbability', ...
    'nPairsWithin10Rmax','nSecondaryWithin10Rmax','secondaryActivationProbabilityWithin10Rmax', ...
    'medianNearestStableGamma','medianMaxAbsDeltaR_over_R0','medianRandomMaxAbsDeltaR_over_R0','medianRdotRmsRatio', ...
    'maxGamma','extendedMaxGamma','secondaryPostCollapseFrames'});
end


% =========================================================================
function [frame, x, y, trackId] = collapse_arrays(collapseData)
frame = nan(0,1);
x = nan(0,1);
y = nan(0,1);
trackId = nan(0,1);
if isempty(collapseData) || ~isstruct(collapseData)
    return;
end
if isfield(collapseData, 'collapseFrame'), frame = collapseData.collapseFrame(:); end
if isfield(collapseData, 'collapseX_mm'), x = collapseData.collapseX_mm(:); end
if isfield(collapseData, 'collapseY_mm'), y = collapseData.collapseY_mm(:); end
if isfield(collapseData, 'collapseTrackId'), trackId = collapseData.collapseTrackId(:); end
n = min([numel(frame), numel(x), numel(y), numel(trackId)]);
frame = frame(1:n);
x = x(1:n);
y = y(1:n);
trackId = trackId(1:n);
valid = isfinite(frame) & isfinite(x) & isfinite(y) & isfinite(trackId);
frame = frame(valid);
x = x(valid);
y = y(valid);
trackId = trackId(valid);
end


% =========================================================================
function v = get_scalar_field(s, fieldName, defaultValue)
v = defaultValue;
if isstruct(s) && isfield(s, fieldName)
    raw = s.(fieldName);
    if ~isempty(raw) && isnumeric(raw) && isfinite(raw(1))
        v = double(raw(1));
    end
end
end


% =========================================================================
function tf = get_bool_field(s, fieldName)
tf = false;
if isstruct(s) && isfield(s, fieldName)
    raw = s.(fieldName);
    if islogical(raw)
        tf = any(raw(:));
    elseif isnumeric(raw)
        tf = any(raw(:) ~= 0);
    end
end
end


% =========================================================================
function [ciLow, ciHigh] = wilson_ci_local(k, n)
ciLow = NaN;
ciHigh = NaN;
if ~(isfinite(k) && isfinite(n) && n > 0 && k >= 0 && k <= n)
    return;
end
z = 1.95996398454005;
phat = k / n;
denom = 1 + z^2 / n;
center = (phat + z^2 / (2 * n)) / denom;
halfWidth = z * sqrt((phat * (1 - phat) / n) + (z^2 / (4 * n^2))) / denom;
ciLow = max(0, center - halfWidth);
ciHigh = min(1, center + halfWidth);
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


% =========================================================================
function y = safe_div(a, b)
if isfinite(a) && isfinite(b) && b > 0
    y = a / b;
else
    y = NaN;
end
end


% =========================================================================
function n = height_or_zero(T)
if istable(T)
    n = height(T);
else
    n = 0;
end
end


% =========================================================================
function n = count_true_table_var(T, varName)
n = 0;
if istable(T) && height(T) > 0 && ismember(varName, T.Properties.VariableNames)
    n = sum(logical(T.(varName)));
end
end
