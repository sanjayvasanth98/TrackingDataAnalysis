function rateData = analyze_growth_collapse_rates(out, metrics, caseDef, opts)
%ANALYZE_GROWTH_COLLAPSE_RATES  Track-based bubble growth/collapse rates.
%
%   rateData = analyze_growth_collapse_rates(out, metrics, caseDef, opts)
%
%   Computes normalized axial-length rates from TrackMate tracks:
%     growth   = dL/dt / U, positive segments only
%     collapse = dL/dt / U, negative post-peak segments only
%
%   Growth is computed from strict activated upstream-moving tracks by
%   default. Collapse is computed from collapse-qualified tracks using the
%   same peak-area and end-shrinkage gates as the collapse analysis.

if nargin < 4 || isempty(opts), opts = struct(); end
opts = apply_defaults(opts);

rateData = make_empty_rate_data(caseDef, opts);

if ~isstruct(out) || ~isfield(out, 'trajectories') || isempty(out.trajectories) || ...
        ~isfield(out, 'spots') || ~istable(out.spots) || height(out.spots) == 0
    return;
end

pixelSize = caseDef.pixelSize; % mm/px
U_mm_s = opts.U_m_s * 1000;
if ~(isfinite(U_mm_s) && U_mm_s > 0)
    warning('analyze_growth_collapse_rates: invalid U_m_s. Skipping.');
    return;
end

[areaById, arById, majorById, minorById] = build_spot_maps(out.spots);
growthTrackIds = select_growth_track_ids(metrics, opts);

growthSeg = make_empty_segments();
collapseSeg = make_empty_segments();

traj = out.trajectories;
for k = 1:numel(traj)
    tk = traj(k);
    if ~isfield(tk, 'spotIds') || ~isfield(tk, 't') || ~isfield(tk, 'frame')
        continue;
    end

    [trackID, frame, time_s, area_px2, arVals, length_px] = ...
        extract_track_series(tk, areaById, arById, majorById, minorById);
    if numel(time_s) < 2 || numel(area_px2) ~= numel(time_s) || numel(length_px) ~= numel(time_s)
        continue;
    end

    valid = isfinite(time_s) & isfinite(area_px2) & area_px2 > 0 & ...
        isfinite(length_px) & length_px > 0;
    if sum(valid) < opts.minTrackSpots
        continue;
    end

    length_mm = length_px * pixelSize;
    dt_s = diff(time_s);
    dL_mm = diff(length_mm);
    rate_mm_s = dL_mm ./ dt_s;
    arMid = pairwise_nanmean(arVals(1:end-1), arVals(2:end));
    frameMid = pairwise_nanmean(frame(1:end-1), frame(2:end));
    lengthMid_mm = pairwise_nanmean(length_mm(1:end-1), length_mm(2:end));

    validRate = isfinite(rate_mm_s) & isfinite(dt_s) & dt_s > 0 & ...
        isfinite(arMid) & isfinite(frameMid) & isfinite(lengthMid_mm);
    if ~any(validRate)
        continue;
    end

    if ismember(trackID, growthTrackIds)
        gMask = validRate & rate_mm_s > 0;
        growthSeg = append_segments(growthSeg, trackID, frameMid(gMask), arMid(gMask), ...
            lengthMid_mm(gMask), rate_mm_s(gMask), rate_mm_s(gMask) ./ U_mm_s);
    end

    if is_collapse_qualified_track(tk, area_px2, opts)
        validLengthIdx = find(isfinite(length_mm) & length_mm > 0);
        if isempty(validLengthIdx)
            continue;
        end
        [~, relPeakIdx] = max(length_mm(validLengthIdx));
        peakIdx = validLengthIdx(relPeakIdx);
        segmentStartIdx = (1:numel(rate_mm_s)).';
        cMask = validRate & segmentStartIdx >= peakIdx & rate_mm_s < 0;
        collapseSeg = append_segments(collapseSeg, trackID, frameMid(cMask), arMid(cMask), ...
            lengthMid_mm(cMask), rate_mm_s(cMask), rate_mm_s(cMask) ./ U_mm_s);
    end
end

rateData.growthSegments = growthSeg;
rateData.collapseSegments = collapseSeg;
rateData.growth = summarize_segments_by_ar(growthSeg, opts.aspectRatioEdges, opts.aspectRatioLabels);
rateData.collapse = summarize_segments_by_ar(collapseSeg, opts.aspectRatioEdges, opts.aspectRatioLabels);
end


% =========================================================================
function opts = apply_defaults(opts)
if ~isfield(opts, 'aspectRatioEdges') || isempty(opts.aspectRatioEdges)
    opts.aspectRatioEdges = [0, 2, 5, Inf];
end
opts.aspectRatioEdges = double(opts.aspectRatioEdges(:)).';
if numel(opts.aspectRatioEdges) < 2
    opts.aspectRatioEdges = [0, Inf];
end
opts.aspectRatioEdges = sort(opts.aspectRatioEdges);

if ~isfield(opts, 'aspectRatioLabels') || isempty(opts.aspectRatioLabels)
    opts.aspectRatioLabels = make_ar_labels(opts.aspectRatioEdges);
else
    opts.aspectRatioLabels = string(opts.aspectRatioLabels(:));
end

if ~isfield(opts, 'U_m_s') || isempty(opts.U_m_s)
    opts.U_m_s = 13.32;
end
if ~isfield(opts, 'growthTrackMode') || isempty(opts.growthTrackMode)
    opts.growthTrackMode = "strictActivated";
end
if ~isfield(opts, 'minTrackSpots') || isempty(opts.minTrackSpots)
    opts.minTrackSpots = 3;
end
if ~isfield(opts, 'collapseMinPeakArea_px2') || isempty(opts.collapseMinPeakArea_px2)
    opts.collapseMinPeakArea_px2 = 50;
end
if ~isfield(opts, 'collapseTruncationFactor') || isempty(opts.collapseTruncationFactor)
    opts.collapseTruncationFactor = 0.6;
end
if ~isfield(opts, 'roiUnwantedMask')
    opts.roiUnwantedMask = [];
end
if ~isfield(opts, 'roiPixelSize') || isempty(opts.roiPixelSize)
    opts.roiPixelSize = 0;
end
end


% =========================================================================
function rateData = make_empty_rate_data(caseDef, opts)
rateData = struct();
rateData.caseName = string(caseDef.name);
rateData.Re = caseDef.Re;
rateData.kD = caseDef.kD;
rateData.pixelSize = caseDef.pixelSize;
rateData.dt = caseDef.dt;
rateData.U_m_s = opts.U_m_s;
rateData.U_mm_s = opts.U_m_s * 1000;
rateData.aspectRatioEdges = opts.aspectRatioEdges;
rateData.aspectRatioLabels = opts.aspectRatioLabels;
rateData.growthTrackMode = string(opts.growthTrackMode);
rateData.growthSegments = make_empty_segments();
rateData.collapseSegments = make_empty_segments();
rateData.growth = summarize_segments_by_ar(rateData.growthSegments, opts.aspectRatioEdges, opts.aspectRatioLabels);
rateData.collapse = summarize_segments_by_ar(rateData.collapseSegments, opts.aspectRatioEdges, opts.aspectRatioLabels);
end


% =========================================================================
function seg = make_empty_segments()
seg = struct();
seg.trackID = nan(0,1);
seg.frame = nan(0,1);
seg.aspectRatio = nan(0,1);
seg.length_mm = nan(0,1);
seg.rate_mm_s = nan(0,1);
seg.rate_over_U = nan(0,1);
end


% =========================================================================
function [areaById, arById, majorById, minorById] = build_spot_maps(spots)
areaById = containers.Map('KeyType', 'double', 'ValueType', 'double');
arById = containers.Map('KeyType', 'double', 'ValueType', 'double');
majorById = containers.Map('KeyType', 'double', 'ValueType', 'double');
minorById = containers.Map('KeyType', 'double', 'ValueType', 'double');

hasArea = ismember('AREA', spots.Properties.VariableNames);
hasAR = ismember('ELLIPSE_ASPECTRATIO', spots.Properties.VariableNames);
hasMajor = ismember('ELLIPSE_MAJOR', spots.Properties.VariableNames);
hasMinor = ismember('ELLIPSE_MINOR', spots.Properties.VariableNames);

for si = 1:height(spots)
    sid = double(spots.ID(si));
    if ~isfinite(sid), continue; end
    if hasArea
        areaById(sid) = double(spots.AREA(si));
    end
    if hasAR
        arById(sid) = double(spots.ELLIPSE_ASPECTRATIO(si));
    end
    if hasMajor
        majorById(sid) = double(spots.ELLIPSE_MAJOR(si));
    end
    if hasMinor
        minorById(sid) = double(spots.ELLIPSE_MINOR(si));
    end
end
end


% =========================================================================
function trackIds = select_growth_track_ids(metrics, opts)
trackIds = nan(0,1);
mode = lower(string(opts.growthTrackMode));

if strcmp(mode, "strictprimary") && isfield(metrics, 'strictPrimaryTrackIds')
    trackIds = metrics.strictPrimaryTrackIds(:);
elseif isfield(metrics, 'strictActivatedTrackIds')
    trackIds = metrics.strictActivatedTrackIds(:);
elseif isfield(metrics, 'activationEvent_trackId')
    trackIds = metrics.activationEvent_trackId(:);
elseif isfield(metrics, 'strictPrimaryTrackIds')
    trackIds = metrics.strictPrimaryTrackIds(:);
end

trackIds = unique(trackIds(isfinite(trackIds)), 'stable');
end


% =========================================================================
function [trackID, frame, time_s, area_px2, arVals, length_px] = ...
    extract_track_series(tk, areaById, arById, majorById, minorById)
trackID = NaN;
if isfield(tk, 'TRACK_ID') && ~isempty(tk.TRACK_ID)
    trackID = double(tk.TRACK_ID);
end

spotIds = tk.spotIds(:);
n = numel(spotIds);
frame = nan(n,1);
time_s = nan(n,1);
area_px2 = nan(n,1);
arVals = nan(n,1);
length_px = nan(n,1);

if isfield(tk, 'frame') && numel(tk.frame) == n
    frame = double(tk.frame(:));
end
if isfield(tk, 't') && numel(tk.t) == n
    time_s = double(tk.t(:));
end

for ii = 1:n
    sid = double(spotIds(ii));
    if isfinite(sid) && isKey(areaById, sid)
        area_px2(ii) = areaById(sid);
    end
    if isfinite(sid) && isKey(arById, sid)
        arVals(ii) = arById(sid);
    end
    if isfinite(sid) && isKey(majorById, sid)
        length_px(ii) = majorById(sid);
    end
    if ~isfinite(arVals(ii)) && isfinite(sid) && isKey(majorById, sid) && isKey(minorById, sid)
        majorVal = majorById(sid);
        minorVal = minorById(sid);
        if isfinite(majorVal) && isfinite(minorVal) && minorVal > 0
            arVals(ii) = majorVal / minorVal;
        end
    end
    if ~isfinite(length_px(ii)) && isfinite(area_px2(ii)) && area_px2(ii) > 0
        % Fallback when TrackMate ellipse features are unavailable.
        length_px(ii) = sqrt(4 * area_px2(ii) / pi);
    end
end
end


% =========================================================================
function out = pairwise_nanmean(a, b)
out = nan(size(a));
validA = isfinite(a);
validB = isfinite(b);
both = validA & validB;
out(both) = 0.5 .* (a(both) + b(both));
onlyA = validA & ~validB;
out(onlyA) = a(onlyA);
onlyB = ~validA & validB;
out(onlyB) = b(onlyB);
end


% =========================================================================
function seg = append_segments(seg, trackID, frame, arVals, length_mm, rate_mm_s, rateNorm)
n = numel(rateNorm);
if n == 0
    return;
end
seg.trackID = [seg.trackID; repmat(trackID, n, 1)];
seg.frame = [seg.frame; frame(:)];
seg.aspectRatio = [seg.aspectRatio; arVals(:)];
seg.length_mm = [seg.length_mm; length_mm(:)];
seg.rate_mm_s = [seg.rate_mm_s; rate_mm_s(:)];
seg.rate_over_U = [seg.rate_over_U; rateNorm(:)];
end


% =========================================================================
function tf = is_collapse_qualified_track(tk, area_px2, opts)
tf = false;

validA = area_px2(isfinite(area_px2) & area_px2 > 0);
if numel(validA) < opts.minTrackSpots
    return;
end

peakA = max(validA);
if peakA < opts.collapseMinPeakArea_px2
    return;
end

lastA = validA(end);
if lastA > opts.collapseTruncationFactor * peakA
    return;
end

if ~isempty(opts.roiUnwantedMask) && isfinite(opts.roiPixelSize) && opts.roiPixelSize > 0 && ...
        isfield(tk, 'x_phys') && isfield(tk, 'y_phys') && ~isempty(tk.x_phys) && ~isempty(tk.y_phys)
    if pt_in_roi(double(tk.x_phys(end)), double(tk.y_phys(end)), opts.roiUnwantedMask, opts.roiPixelSize)
        return;
    end
end

tf = true;
end


% =========================================================================
function result = pt_in_roi(x_mm, y_mm_image, mask, pixelSize)
[nRows, nCols] = size(mask);
c = round(x_mm / pixelSize);
r = round(y_mm_image / pixelSize);
result = isfinite(c) && isfinite(r) && c >= 1 && c <= nCols && r >= 1 && r <= nRows && mask(r, c);
end


% =========================================================================
function summary = summarize_segments_by_ar(seg, arEdges, arLabels)
nAR = max(numel(arEdges) - 1, 0);
summary = struct();
summary.aspectRatioEdges = arEdges(:);
summary.aspectRatioLabels = arLabels(:);
summary.nSegments = zeros(nAR, 1);
summary.nTracks = zeros(nAR, 1);
summary.rateMedian_over_U = nan(nAR, 1);
summary.rateMean_over_U = nan(nAR, 1);
summary.rateStd_over_U = nan(nAR, 1);
summary.rateSem_over_U = nan(nAR, 1);
summary.rateMedian_mm_s = nan(nAR, 1);

for ai = 1:nAR
    lo = arEdges(ai);
    hi = arEdges(ai + 1);
    if isinf(hi)
        arMask = seg.aspectRatio >= lo;
    else
        arMask = seg.aspectRatio >= lo & seg.aspectRatio < hi;
    end
    mask = isfinite(seg.rate_over_U) & isfinite(seg.aspectRatio) & arMask;
    vals = seg.rate_over_U(mask);
    rawVals = seg.rate_mm_s(mask);
    ids = seg.trackID(mask);
    vals = vals(isfinite(vals));
    rawVals = rawVals(isfinite(rawVals));

    summary.nSegments(ai) = numel(vals);
    summary.nTracks(ai) = numel(unique(ids(isfinite(ids))));
    if isempty(vals)
        continue;
    end
    summary.rateMedian_over_U(ai) = median(vals);
    summary.rateMean_over_U(ai) = mean(vals);
    if numel(vals) > 1
        summary.rateStd_over_U(ai) = std(vals, 0);
        summary.rateSem_over_U(ai) = summary.rateStd_over_U(ai) / sqrt(numel(vals));
    else
        summary.rateStd_over_U(ai) = 0;
        summary.rateSem_over_U(ai) = NaN;
    end
    if ~isempty(rawVals)
        summary.rateMedian_mm_s(ai) = median(rawVals);
    end
end
end


% =========================================================================
function labels = make_ar_labels(edges)
n = max(numel(edges) - 1, 0);
labels = strings(n, 1);
for i = 1:n
    lo = edges(i);
    hi = edges(i + 1);
    if isinf(hi)
        labels(i) = sprintf('AR >= %.1f', lo);
    elseif lo <= 0
        labels(i) = sprintf('AR < %.1f', hi);
    else
        labels(i) = sprintf('%.1f <= AR < %.1f', lo, hi);
    end
end
end
