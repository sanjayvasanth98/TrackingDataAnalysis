function save_video_overlay_gif_from_avi(caseDef, metrics, outDir, plotOpts)

if ~isstruct(metrics) || ~isfield(metrics, 'trackCatalog') || isempty(metrics.trackCatalog)
    warning('Track catalog missing for case %s. Skipping AVI-overlay GIF export.', char(caseDef.name));
    return;
end

if ~isfield(plotOpts, 'saveVideoOverlayGifs') || ~logical(plotOpts.saveVideoOverlayGifs)
    return;
end

if ~isfolder(outDir)
    mkdir(outDir);
end

videoFile = resolve_case_video_file(caseDef);
if strlength(videoFile) < 1 || ~isfile(videoFile)
    warning('AVI file not found for case %s. Set cases(i).videoFile to a valid .avi path.', char(caseDef.name));
    return;
end

trackCatalog = metrics.trackCatalog;
trackRequests = [];
if isfield(caseDef, 'diagnosticTrackIds') && ~isempty(caseDef.diagnosticTrackIds)
    trackRequests = caseDef.diagnosticTrackIds(:).';
end
[selectedIdx, selectedTrackIds] = resolve_diagnostic_track_indices(trackCatalog, trackRequests, caseDef, 'AVI overlay GIF');
if isempty(selectedIdx)
    warning('AVI overlay GIF: no tracks selected for case %s.', char(caseDef.name));
    return;
end

[leftIdx, leftTrackIds] = resolve_leftmoving_subset(trackCatalog, selectedIdx, selectedTrackIds, metrics);
[microIdx, microTrackIds] = resolve_microbubble_subset(trackCatalog, selectedIdx, selectedTrackIds, metrics);
if isempty(leftIdx) && isempty(microIdx)
    warning('AVI overlay GIF: no left-moving or microbubble-rescue tracks selected for case %s.', char(caseDef.name));
    return;
end

[actXY, actFrames, actTrackIds] = activation_events_left_and_micro(metrics);
if ~isempty(actTrackIds)
    finiteIds = isfinite(actTrackIds);
    if any(finiteIds)
        keepIds = unique([leftTrackIds(:); microTrackIds(:)]);
        keepAct = ismember(actTrackIds, keepIds);
        actXY = actXY(keepAct, :);
        actFrames = actFrames(keepAct);
        actTrackIds = actTrackIds(keepAct);
    end
end

drawIdx = unique([leftIdx(:); microIdx(:)]);
drawTrackFrames = collect_track_frames(trackCatalog, drawIdx);
allTrackFrames = [drawTrackFrames(:); actFrames(:)];
allTrackFrames = allTrackFrames(isfinite(allTrackFrames));
if isempty(allTrackFrames)
    warning('AVI overlay GIF: no finite frame values for case %s.', char(caseDef.name));
    return;
end

drawTrackIds = nan(numel(drawIdx), 1);
for t = 1:numel(drawIdx)
    drawTrackIds(t) = get_track_id(trackCatalog(drawIdx(t)));
end
drawStrictMask = ismember(drawTrackIds, leftTrackIds);
drawMicroMask = ismember(drawTrackIds, microTrackIds);
microBlueState = false(numel(drawIdx), 1);

isZeroBased = min(allTrackFrames) <= 0;

v = VideoReader(char(videoFile));
nFramesEst = max(1, floor(v.Duration * v.FrameRate));

trackFrameIdx = map_track_frames_to_video_idx(drawTrackFrames, isZeroBased, nFramesEst);
actFrameIdx = map_track_frames_to_video_idx(actFrames, isZeroBased, nFramesEst);
baseFrameIdx = unique([trackFrameIdx(:); actFrameIdx(:)]);
baseFrameIdx = baseFrameIdx(isfinite(baseFrameIdx));
if isempty(baseFrameIdx)
    warning('AVI overlay GIF: no mappable video frames for case %s.', char(caseDef.name));
    return;
end

useContiguous = false;
if isfield(plotOpts, 'videoOverlayUseContiguousRange')
    useContiguous = logical(plotOpts.videoOverlayUseContiguousRange);
end
if useContiguous
    frameIdxToRender = (min(baseFrameIdx):max(baseFrameIdx)).';
else
    frameIdxToRender = unique(baseFrameIdx(:));
end

maxFrames = Inf;
if isfield(plotOpts, 'videoOverlayMaxFrames') && ~isempty(plotOpts.videoOverlayMaxFrames)
    maxFrames = double(plotOpts.videoOverlayMaxFrames);
end
if isfinite(maxFrames) && maxFrames > 0 && numel(frameIdxToRender) > maxFrames
    step = ceil(numel(frameIdxToRender) / maxFrames);
    frameIdxToRender = frameIdxToRender(1:step:end);
    if frameIdxToRender(end) ~= baseFrameIdx(end)
        frameIdxToRender(end+1,1) = baseFrameIdx(end); %#ok<AGROW>
    end
end
frameIdxToRender = unique(frameIdxToRender(:));

trailLength = get_numeric_opt(plotOpts, 'videoOverlayTrailLength', 10);
trailLength = max(1, round(trailLength));
delayTime = get_numeric_opt(plotOpts, 'videoOverlayDelayTime', 0.08);
delayTime = max(0.01, delayTime);
fadeHalfLifeFrames = get_numeric_opt(plotOpts, 'videoOverlayFadeHalfLifeFrames', 30);
fadeHalfLifeFrames = max(1, fadeHalfLifeFrames);
markerSize = get_numeric_opt(plotOpts, 'videoOverlayMarkerSize', 26);
markerSize = max(4, markerSize);

xLim = [0 4.8];
if isfield(plotOpts, 'videoOverlayXLim_mm') && numel(plotOpts.videoOverlayXLim_mm) >= 2
    xLim = double(plotOpts.videoOverlayXLim_mm(1:2));
elseif isfield(plotOpts, 'inceptionXLim_mm') && numel(plotOpts.inceptionXLim_mm) >= 2
    xLim = double(plotOpts.inceptionXLim_mm(1:2));
end
if ~all(isfinite(xLim)) || xLim(2) <= xLim(1)
    xLim = [0 4.8];
end
yLim = [0 1.2];
if isfield(plotOpts, 'inceptionYLim_mm') && numel(plotOpts.inceptionYLim_mm) >= 2
    yLim = double(plotOpts.inceptionYLim_mm(1:2));
end
xTicks = fixed_ticks(xLim(1), xLim(2), 0.5);
yTicks = fixed_ticks(yLim(1), yLim(2), 0.2);

figPos = [120 120 1280 320];
if isfield(plotOpts, 'inceptionImageSize_px') && numel(plotOpts.inceptionImageSize_px) >= 2
    sz = double(plotOpts.inceptionImageSize_px(1:2));
    figPos = [120 120 round(sz(1)) round(sz(2))];
end

gifPath = fullfile(outDir, sprintf('%s_Re_%g_kDh_%g_video_overlay.gif', char(caseDef.name), caseDef.Re, caseDef.kDh));

f = figure('Color', 'w', 'Position', figPos, 'Visible', 'off');
ax = axes(f);
savedFrames = 0;
for ii = 1:numel(frameIdxToRender)
    frameIdx = frameIdxToRender(ii);
    frameRgb = safe_read_video_frame(v, frameIdx);
    if isempty(frameRgb)
        continue;
    end

    [hPx, wPx, ~] = size(frameRgb);
    yExtent_mm = hPx * caseDef.pixelSize;
    currentTrackFrame = video_idx_to_track_frame(frameIdx, isZeroBased);

    cla(ax);
    hold(ax, 'on');

    xVec = (0:wPx-1) * caseDef.pixelSize;
    yVec = (0:hPx-1) * caseDef.pixelSize;
    image(ax, ...
        'XData', [xVec(1), xVec(end)], ...
        'YData', [yVec(1), yVec(end)], ...
        'CData', flipud(frameRgb));
    set(ax, 'YDir', 'normal');

    for t = 1:numel(drawIdx)
        tr = trackCatalog(drawIdx(t));
        [xTail, yTailPlot] = tail_xy_at_frame(tr, currentTrackFrame, trailLength, yExtent_mm);
        if isempty(xTail)
            continue;
        end

        if drawStrictMask(t)
            cLine = [0 0 1];
        elseif drawMicroMask(t)
            motionSignal = recent_motion_signal(xTail, 3);
            if motionSignal < 0
                microBlueState(t) = true;
            elseif motionSignal > 0
                microBlueState(t) = false;
            end

            if microBlueState(t)
                cLine = [0 0 1];
            else
                cLine = [0 0.60 0.10];
            end
        else
            cLine = [0.4 0.4 0.4];
        end
        plot(ax, xTail, yTailPlot, '-', 'Color', cLine, 'LineWidth', 1.0);
    end

    for a = 1:size(actXY, 1)
        fAct = actFrames(a);
        if ~isfinite(fAct)
            continue;
        end
        age = currentTrackFrame - fAct;
        if age < 0
            continue;
        end
        alpha = exp(log(0.5) * (age / fadeHalfLifeFrames));
        if alpha < 0.04
            continue;
        end
        cAct = alpha * [0.9 0.15 0.1] + (1 - alpha) * [1 1 1];
        scatter(ax, actXY(a,1), yExtent_mm - actXY(a,2), markerSize, ...
            'Marker', 'o', ...
            'MarkerFaceColor', cAct, ...
            'MarkerEdgeColor', 'none');
    end

    xlim(ax, xLim);
    ylim(ax, yLim);
    set(ax, ...
        'XTick', xTicks, ...
        'YTick', yTicks, ...
        'XLimMode', 'manual', ...
        'YLimMode', 'manual', ...
        'XTickMode', 'manual', ...
        'YTickMode', 'manual', ...
        'DataAspectRatioMode', 'auto', ...
        'PlotBoxAspectRatioMode', 'auto');
    xlabel(ax, '$x\;(\mathrm{mm})$', 'Interpreter', 'latex');
    ylabel(ax, '$y\;(\mathrm{mm})$', 'Interpreter', 'latex');
    grid(ax, 'on');
    box(ax, 'on');

    drawnow;
    fr = getframe(f);
    rgb = frame2im(fr);
    if size(rgb,3) < 3
        rgb = repmat(rgb, [1 1 3]);
    end
    [imInd, map] = rgb2ind(rgb, 256);
    if savedFrames == 0
        imwrite(imInd, map, gifPath, 'gif', 'LoopCount', Inf, 'DelayTime', delayTime);
    else
        imwrite(imInd, map, gifPath, 'gif', 'WriteMode', 'append', 'DelayTime', delayTime);
    end
    savedFrames = savedFrames + 1;
end

close(f);
if savedFrames == 0
    warning('AVI overlay GIF: no frames were saved for case %s.', char(caseDef.name));
    return;
end
fprintf('Saved AVI-overlay GIF: %s (%d frames)\n', gifPath, savedFrames);
end

function [leftIdx, leftTrackIds] = resolve_leftmoving_subset(trackCatalog, selectedIdx, selectedTrackIds, metrics)
leftIdx = zeros(0,1);
leftTrackIds = nan(0,1);

metricIds = nan(0,1);
if isfield(metrics, 'strictPrimaryTrackIds') && ~isempty(metrics.strictPrimaryTrackIds)
    metricIds = metrics.strictPrimaryTrackIds(:);
elseif isfield(metrics, 'leftMovingTrackIds_netLeftLegacy') && ~isempty(metrics.leftMovingTrackIds_netLeftLegacy)
    metricIds = metrics.leftMovingTrackIds_netLeftLegacy(:);
end
metricIds = unique(metricIds(isfinite(metricIds)), 'stable');

if ~isempty(metricIds)
    mask = ismember(selectedTrackIds, metricIds);
else
    mask = false(numel(selectedIdx),1);
    hasStrictFlag = false;
    for i = 1:numel(selectedIdx)
        hasStrictFlag = hasStrictFlag || is_true_field(trackCatalog(selectedIdx(i)), 'isStrictPrimary');
    end
    for i = 1:numel(selectedIdx)
        tr = trackCatalog(selectedIdx(i));
        if hasStrictFlag
            mask(i) = is_true_field(tr, 'isStrictPrimary');
        else
            mask(i) = is_true_field(tr, 'isLeftMoving');
        end
    end
end

leftIdx = selectedIdx(mask);
leftTrackIds = selectedTrackIds(mask);
leftTrackIds = unique(leftTrackIds(isfinite(leftTrackIds)), 'stable');
end

function [microIdx, microTrackIds] = resolve_microbubble_subset(trackCatalog, selectedIdx, selectedTrackIds, metrics)
microIdx = zeros(0,1);
microTrackIds = nan(0,1);

metricIds = nan(0,1);
if isfield(metrics, 'microbubbleRescueTrackIds_nonLeft') && ~isempty(metrics.microbubbleRescueTrackIds_nonLeft)
    metricIds = metrics.microbubbleRescueTrackIds_nonLeft(:);
end
metricIds = unique(metricIds(isfinite(metricIds)), 'stable');
if isempty(metricIds)
    return;
end

mask = ismember(selectedTrackIds, metricIds);
microIdx = selectedIdx(mask);
microTrackIds = selectedTrackIds(mask);
microTrackIds = unique(microTrackIds(isfinite(microTrackIds)), 'stable');
end

function [actXY, actFrames, actTrackIds] = activation_events_left_and_micro(metrics)
actXY = zeros(0,2);
actFrames = nan(0,1);
actTrackIds = nan(0,1);

[xyA, frA, idA] = get_activation_triplet(metrics, ...
    'strictActivationEvent_xy', ...
    'strictActivationEvent_frame', ...
    'strictActivationEvent_trackId');
if isempty(xyA)
    [xyA, frA, idA] = get_activation_triplet(metrics, ...
        'activationEvent_xy', ...
        'activationEvent_frame', ...
        'activationEvent_trackId');
end
if isempty(xyA)
    [xyA, frA, idA] = get_activation_triplet(metrics, ...
        'activationEvent_xy_netLeftLegacy', ...
        'activationEvent_frame_netLeftLegacy', ...
        'activationEvent_trackId_netLeftLegacy');
end

[xyB, frB, idB] = get_activation_triplet(metrics, ...
    'microbubbleActivationEvent_xy_nonLeft', ...
    'microbubbleActivationEvent_frame_nonLeft', ...
    'microbubbleActivationEvent_trackId_nonLeft');

[actXY, actFrames, actTrackIds] = concat_unique_events(xyA, frA, idA, xyB, frB, idB);
end

function frameVals = collect_track_frames(trackCatalog, trackIdx)
frameVals = nan(0,1);
for i = 1:numel(trackIdx)
    f = trackCatalog(trackIdx(i)).frame(:);
    f = f(isfinite(f));
    if ~isempty(f)
        frameVals = [frameVals; f]; %#ok<AGROW>
    end
end
end

function idx = map_track_frames_to_video_idx(frameVals, isZeroBased, nFrames)
idx = round(frameVals(:));
if isZeroBased
    idx = idx + 1;
end
idx(~isfinite(idx)) = NaN;
idx = max(1, min(nFrames, idx));
idx = unique(idx(isfinite(idx)));
end

function frameVal = video_idx_to_track_frame(frameIdx, isZeroBased)
if isZeroBased
    frameVal = frameIdx - 1;
else
    frameVal = frameIdx;
end
end

function [xTail, yTailPlot] = tail_xy_at_frame(tr, frameVal, trailLength, yExtent_mm)
xTail = nan(0,1);
yTailPlot = nan(0,1);

f = tr.frame(:);
x = tr.x(:);
y = tr.y(:);
if isempty(f) || isempty(x) || isempty(y) || numel(f) ~= numel(x) || numel(y) ~= numel(x)
    return;
end
valid = isfinite(f) & isfinite(x) & isfinite(y);
if ~any(valid)
    return;
end
f = f(valid);
x = x(valid);
y = y(valid);

[f, ord] = sort(f);
x = x(ord);
y = y(ord);

if frameVal < f(1) || frameVal > f(end)
    return;
end

idxNow = find(f <= frameVal, 1, 'last');
if isempty(idxNow)
    return;
end
tailStart = max(1, idxNow - trailLength + 1);
xTail = x(tailStart:idxNow);
yTailPlot = yExtent_mm - y(tailStart:idxNow);
end

function rgb = safe_read_video_frame(v, frameIdx)
rgb = [];
try
    frm = read(v, frameIdx);
catch
    try
        t = max(0, (frameIdx - 1) / max(v.FrameRate, eps));
        v.CurrentTime = min(max(0, t), max(0, v.Duration - eps));
        if hasFrame(v)
            frm = readFrame(v);
        else
            frm = [];
        end
    catch
        frm = [];
    end
end
if isempty(frm)
    return;
end
rgb = ensure_uint8_rgb(frm);
end

function rgb = ensure_uint8_rgb(frm)
if ndims(frm) == 2
    frm = repmat(frm, [1 1 3]);
elseif ndims(frm) == 3 && size(frm,3) == 1
    frm = repmat(frm, [1 1 3]);
end

if isa(frm, 'uint8')
    rgb = frm;
    return;
end
if isa(frm, 'uint16')
    rgb = uint8(double(frm) / 257);
    return;
end

frm = double(frm);
if isempty(frm)
    rgb = uint8([]);
    return;
end
if max(frm(:)) <= 1.0
    frm = frm * 255;
end
frm = max(0, min(255, frm));
rgb = uint8(round(frm));
end

function s = resolve_case_video_file(caseDef)
s = "";
if isfield(caseDef, 'videoFile') && strlength(string(caseDef.videoFile)) > 0
    s = string(caseDef.videoFile);
    return;
end

if ~isfield(caseDef, 'xmlFile') || strlength(string(caseDef.xmlFile)) < 1
    return;
end
[p, n] = fileparts(char(string(caseDef.xmlFile)));
candidates = string({fullfile(p, [n '.avi']), fullfile(p, [n '.AVI'])});
for i = 1:numel(candidates)
    if isfile(candidates(i))
        s = candidates(i);
        return;
    end
end
end

function out = get_numeric_opt(opts, fieldName, defaultVal)
out = defaultVal;
if isfield(opts, fieldName)
    v = opts.(fieldName);
    if isscalar(v) && ~isempty(v) && isfinite(v)
        out = double(v);
    end
end
end

function tid = get_track_id(tr)
tid = NaN;
if isstruct(tr) && isfield(tr, 'TRACK_ID') && isfinite(tr.TRACK_ID)
    tid = tr.TRACK_ID;
end
end

function [xy, fr, tid] = get_activation_triplet(metrics, fieldXY, fieldFr, fieldTid)
xy = zeros(0,2);
fr = nan(0,1);
tid = nan(0,1);

if isfield(metrics, fieldXY) && ~isempty(metrics.(fieldXY))
    xy = metrics.(fieldXY);
end
if isfield(metrics, fieldFr) && ~isempty(metrics.(fieldFr))
    fr = metrics.(fieldFr)(:);
end
if isfield(metrics, fieldTid) && ~isempty(metrics.(fieldTid))
    tid = metrics.(fieldTid)(:);
end

n = min(size(xy,1), numel(fr));
if isempty(tid)
    tid = nan(n,1);
else
    n = min(n, numel(tid));
end
if n < 1
    xy = zeros(0,2);
    fr = nan(0,1);
    tid = nan(0,1);
    return;
end
xy = xy(1:n, :);
fr = fr(1:n);
tid = tid(1:n);
end

function [xyOut, frOut, tidOut] = concat_unique_events(xyA, frA, tidA, xyB, frB, tidB)
xy = [xyA; xyB];
fr = [frA; frB];
tid = [tidA; tidB];

n = min([size(xy,1), numel(fr), numel(tid)]);
if n < 1
    xyOut = zeros(0,2);
    frOut = nan(0,1);
    tidOut = nan(0,1);
    return;
end

xy = xy(1:n, :);
fr = fr(1:n);
tid = tid(1:n);

keep = true(n,1);
seen = containers.Map('KeyType', 'char', 'ValueType', 'logical');
for i = 1:n
    if isfinite(tid(i)) && isfinite(fr(i))
        key = sprintf('tid:%0.0f|fr:%0.0f', tid(i), fr(i));
    else
        key = sprintf('xy:%0.6f|%0.6f|fr:%0.6f', xy(i,1), xy(i,2), fr(i));
    end
    if isKey(seen, key)
        keep(i) = false;
    else
        seen(key) = true;
    end
end

xyOut = xy(keep, :);
frOut = fr(keep);
tidOut = tid(keep);
end

function ticks = fixed_ticks(lo, hi, step)
if ~(isfinite(lo) && isfinite(hi) && isfinite(step) && step > 0)
    ticks = [];
    return;
end
if hi < lo
    t = lo;
    lo = hi;
    hi = t;
end
ticks = lo:step:hi;
if isempty(ticks) || ticks(end) < hi
    ticks = [ticks, hi]; %#ok<AGROW>
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

function motionSignal = recent_motion_signal(xVals, nSteps)
% motionSignal: -1 => last nSteps are leftward, +1 => last nSteps are rightward, 0 => no switch.
motionSignal = 0;
if nargin < 2 || ~isfinite(nSteps) || nSteps < 1
    nSteps = 3;
end
nSteps = round(nSteps);

xVals = xVals(:);
xVals = xVals(isfinite(xVals));
if numel(xVals) < (nSteps + 1)
    return;
end

dx = diff(xVals);
if numel(dx) < nSteps
    return;
end
recent = dx(end - nSteps + 1:end);
if all(recent < 0)
    motionSignal = -1;
elseif all(recent > 0)
    motionSignal = 1;
end
end
