function recircData = analyze_collapse_recirculation(out, metrics, collapseData, caseDef, opts)
%ANALYZE_COLLAPSE_RECIRCULATION  Attribute newly born microbubble tracks to collapse events.
%
% A collapse-generated microbubble is defined here as a small-start track
% born within a short frame window after a qualified collapse and within a
% local distance of the collapse location.

if nargin < 5 || isempty(opts)
    opts = struct();
end
opts = apply_defaults(opts, metrics);

[collapseFrame, collapseX, collapseY, collapseTrackId] = collapse_arrays(collapseData);
nCollapse = numel(collapseFrame);

[trackInfo, activationEventIds] = track_info_from_metrics(metrics);
if isempty(trackInfo) || nCollapse == 0
    recircData = pack_result(caseDef, opts, nCollapse, table(), metrics, activationEventIds);
    return;
end

matchRows = {};
for ti = 1:numel(trackInfo)
    tr = trackInfo(ti);
    if opts.requireBasicValid && ~tr.isBasicValid
        continue;
    end
    if ~isfinite(tr.startFrame) || ~isfinite(tr.startX_mm) || ~isfinite(tr.startY_mm)
        continue;
    end
    if ~isfinite(tr.startArea_px2) || tr.startArea_px2 < opts.microbubbleStartAreaRange_px2(1) || ...
            tr.startArea_px2 > opts.microbubbleStartAreaRange_px2(2)
        continue;
    end

    bestIdx = NaN;
    bestDist = Inf;
    bestLag = Inf;
    for ci = 1:nCollapse
        if isfinite(tr.TRACK_ID) && isfinite(collapseTrackId(ci)) && tr.TRACK_ID == collapseTrackId(ci)
            continue;
        end
        lag = tr.startFrame - collapseFrame(ci);
        if ~isfinite(lag) || lag < opts.minFrameLag || lag > opts.maxFrameLag
            continue;
        end
        dist = hypot(tr.startX_mm - collapseX(ci), tr.startY_mm - collapseY(ci));
        if ~isfinite(dist) || dist > opts.maxDistance_mm
            continue;
        end
        if dist < bestDist || (abs(dist - bestDist) <= eps && lag < bestLag)
            bestIdx = ci;
            bestDist = dist;
            bestLag = lag;
        end
    end

    if ~isfinite(bestIdx)
        continue;
    end

    isActivationEvent = false;
    if ~isempty(activationEventIds) && isfinite(tr.TRACK_ID)
        isActivationEvent = any(activationEventIds == tr.TRACK_ID);
    elseif isempty(activationEventIds)
        isActivationEvent = tr.isActivated;
    end

    matchRows(end+1,:) = { ... %#ok<AGROW>
        string(caseDef.name), caseDef.Re, caseDef.kD, ...
        bestIdx, collapseFrame(bestIdx), collapseX(bestIdx), collapseY(bestIdx), collapseTrackId(bestIdx), ...
        tr.TRACK_ID, tr.startFrame, bestLag, bestDist, tr.startArea_px2, ...
        tr.isStrictPrimary, isActivationEvent, tr.isStrictActivated};
end

if isempty(matchRows)
    matchTable = table();
else
    matchTable = cell2table(matchRows, 'VariableNames', { ...
        'Case','Re','kD', ...
        'collapseIndex','collapseFrame','collapseX_mm','collapseY_mm','collapseTrackId', ...
        'generatedTrackId','generatedStartFrame','frameLag','distance_mm','startArea_px2', ...
        'isLeftMovingRecirculating','isActivationEvent','isStrictActivated'});
end

recircData = pack_result(caseDef, opts, nCollapse, matchTable, metrics, activationEventIds);
end


% =========================================================================
function opts = apply_defaults(opts, metrics)
if ~isfield(opts, 'minFrameLag') || isempty(opts.minFrameLag) || ~isfinite(opts.minFrameLag)
    opts.minFrameLag = 0;
end
if ~isfield(opts, 'maxFrameLag') || isempty(opts.maxFrameLag) || ~isfinite(opts.maxFrameLag)
    opts.maxFrameLag = 3;
end
opts.minFrameLag = round(opts.minFrameLag);
opts.maxFrameLag = max(opts.minFrameLag, round(opts.maxFrameLag));

if ~isfield(opts, 'maxDistance_mm') || isempty(opts.maxDistance_mm) || ~isfinite(opts.maxDistance_mm)
    opts.maxDistance_mm = 0.15;
end
opts.maxDistance_mm = max(opts.maxDistance_mm, 0);

if ~isfield(opts, 'microbubbleStartAreaRange_px2') || numel(opts.microbubbleStartAreaRange_px2) < 2
    if isfield(metrics, 'microbubbleActivationAreaRange_px2') && numel(metrics.microbubbleActivationAreaRange_px2) >= 2
        opts.microbubbleStartAreaRange_px2 = metrics.microbubbleActivationAreaRange_px2(1:2);
    else
        opts.microbubbleStartAreaRange_px2 = [1 120];
    end
end
opts.microbubbleStartAreaRange_px2 = sort(double(opts.microbubbleStartAreaRange_px2(1:2)));

if ~isfield(opts, 'requireBasicValid') || isempty(opts.requireBasicValid)
    opts.requireBasicValid = true;
end
opts.requireBasicValid = logical(opts.requireBasicValid);
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

n = min([numel(frame), numel(x), numel(y)]);
frame = frame(1:n);
x = x(1:n);
y = y(1:n);
if numel(trackId) < n
    trackId(end+1:n,1) = NaN;
else
    trackId = trackId(1:n);
end
valid = isfinite(frame) & isfinite(x) & isfinite(y);
frame = frame(valid);
x = x(valid);
y = y(valid);
trackId = trackId(valid);
end


% =========================================================================
function [trackInfo, activationEventIds] = track_info_from_metrics(metrics)
trackInfo = make_track_info_template();
trackInfo = trackInfo([]);
activationEventIds = activation_ids_from_metrics(metrics);
if ~isstruct(metrics) || ~isfield(metrics, 'trackCatalog') || isempty(metrics.trackCatalog)
    return;
end

catalog = metrics.trackCatalog(:);
trackInfo = repmat(make_track_info_template(), numel(catalog), 1);
for k = 1:numel(catalog)
    tr = catalog(k);
    trackInfo(k).TRACK_ID = get_scalar_field(tr, 'TRACK_ID', NaN);
    [trackInfo(k).startFrame, trackInfo(k).startX_mm, trackInfo(k).startY_mm] = ...
        first_track_point(tr);
    trackInfo(k).startArea_px2 = get_scalar_field(tr, 'microbubbleSeedArea_px2', NaN);
    trackInfo(k).isBasicValid = get_bool_field(tr, 'isBasicValid');
    trackInfo(k).isStrictPrimary = get_bool_field(tr, 'isStrictPrimary');
    trackInfo(k).isStrictActivated = get_bool_field(tr, 'isStrictActivated');
    trackInfo(k).isActivated = get_bool_field(tr, 'isActivated');
end
end


% =========================================================================
function tpl = make_track_info_template()
tpl = struct( ...
    'TRACK_ID', NaN, ...
    'startFrame', NaN, ...
    'startX_mm', NaN, ...
    'startY_mm', NaN, ...
    'startArea_px2', NaN, ...
    'isBasicValid', false, ...
    'isStrictPrimary', false, ...
    'isStrictActivated', false, ...
    'isActivated', false);
end


% =========================================================================
function ids = activation_ids_from_metrics(metrics)
ids = nan(0,1);
if ~isstruct(metrics)
    return;
end
fields = {'strictActivationEvent_trackId', 'microbubbleActivationEvent_trackId_nonLeft'};
for i = 1:numel(fields)
    if isfield(metrics, fields{i}) && ~isempty(metrics.(fields{i}))
        ids = [ids; metrics.(fields{i})(:)]; %#ok<AGROW>
    end
end
ids = unique(ids(isfinite(ids)), 'stable');
end


% =========================================================================
function [frame0, x0, y0] = first_track_point(tr)
frame0 = NaN;
x0 = NaN;
y0 = NaN;
if ~isstruct(tr) || ~isfield(tr, 'frame') || ~isfield(tr, 'x') || ~isfield(tr, 'y')
    return;
end
frame = tr.frame(:);
x = tr.x(:);
y = tr.y(:);
n = min([numel(frame), numel(x), numel(y)]);
if n < 1
    return;
end
valid = isfinite(frame(1:n)) & isfinite(x(1:n)) & isfinite(y(1:n));
idx = find(valid, 1, 'first');
if isempty(idx)
    return;
end
frame0 = frame(idx);
x0 = x(idx);
y0 = y(idx);
end


% =========================================================================
function recircData = pack_result(caseDef, opts, nCollapse, matchTable, metrics, activationEventIds)
nGenerated = height_or_zero(matchTable);
if nGenerated > 0
    nCollapseWithGenerated = numel(unique(matchTable.collapseIndex(isfinite(matchTable.collapseIndex))));
    nGeneratedLeftMoving = sum(logical(matchTable.isLeftMovingRecirculating));
    nGeneratedActivation = sum(logical(matchTable.isActivationEvent));
    nGeneratedStrictActivated = sum(logical(matchTable.isStrictActivated));
else
    nCollapseWithGenerated = 0;
    nGeneratedLeftMoving = 0;
    nGeneratedActivation = 0;
    nGeneratedStrictActivated = 0;
end

nTotalLeftMoving = get_metric_scalar(metrics, 'nStrictPrimaryTracks', 'nLeftMovingTracks', 0);
nTotalActivationEvents = get_metric_scalar(metrics, 'strictActivationEventsTotal', 'activationEventsTotal', 0);
if isstruct(metrics) && isfield(metrics, 'microbubbleActivationEvent_trackId_nonLeft') && ~isempty(metrics.microbubbleActivationEvent_trackId_nonLeft)
    nTotalActivationEvents = nTotalActivationEvents + numel(metrics.microbubbleActivationEvent_trackId_nonLeft);
elseif ~isempty(activationEventIds)
    nTotalActivationEvents = max(nTotalActivationEvents, numel(activationEventIds));
end

recircData = struct();
recircData.caseName = string(caseDef.name);
recircData.Re = caseDef.Re;
recircData.kD = caseDef.kD;
recircData.nCollapseEvents = nCollapse;
recircData.nCollapseEventsWithGeneratedMicrobubbles = nCollapseWithGenerated;
recircData.nCollapseGeneratedMicrobubbles = nGenerated;
recircData.nCollapseGeneratedLeftMoving = nGeneratedLeftMoving;
recircData.nCollapseGeneratedActivationEvents = nGeneratedActivation;
recircData.nCollapseGeneratedStrictActivated = nGeneratedStrictActivated;
recircData.nTotalLeftMovingTracks = nTotalLeftMoving;
recircData.nTotalActivationEvents = nTotalActivationEvents;
recircData.collapseGeneratedPerCollapse = safe_div(nGenerated, nCollapse);
recircData.pctCollapseEventsWithGeneratedMicrobubbles = safe_pct(nCollapseWithGenerated, nCollapse);
recircData.pctCollapseGeneratedLeftMoving_ofGenerated = safe_pct(nGeneratedLeftMoving, nGenerated);
recircData.pctLeftMovingFromCollapseEvents = safe_pct(nGeneratedLeftMoving, nTotalLeftMoving);
recircData.pctCollapseGeneratedActivation_ofGenerated = safe_pct(nGeneratedActivation, nGenerated);
recircData.pctActivationFromCollapseEvents = safe_pct(nGeneratedActivation, nTotalActivationEvents);
recircData.matchMinFrameLag = opts.minFrameLag;
recircData.matchMaxFrameLag = opts.maxFrameLag;
recircData.matchMaxDistance_mm = opts.maxDistance_mm;
recircData.microbubbleStartAreaMin_px2 = opts.microbubbleStartAreaRange_px2(1);
recircData.microbubbleStartAreaMax_px2 = opts.microbubbleStartAreaRange_px2(2);
recircData.matchTable = matchTable;
end


% =========================================================================
function n = height_or_zero(tbl)
if istable(tbl)
    n = height(tbl);
else
    n = 0;
end
end


% =========================================================================
function val = get_metric_scalar(metrics, fieldPrimary, fieldFallback, defaultVal)
val = defaultVal;
if ~isstruct(metrics)
    return;
end
if isfield(metrics, fieldPrimary) && ~isempty(metrics.(fieldPrimary))
    raw = metrics.(fieldPrimary);
elseif isfield(metrics, fieldFallback) && ~isempty(metrics.(fieldFallback))
    raw = metrics.(fieldFallback);
else
    return;
end
raw = raw(:);
raw = raw(isfinite(raw));
if ~isempty(raw)
    val = raw(1);
end
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
function out = safe_div(num, den)
if ~isfinite(den) || den <= 0
    out = NaN;
else
    out = num ./ den;
end
end


% =========================================================================
function out = safe_pct(num, den)
out = 100 * safe_div(num, den);
end
