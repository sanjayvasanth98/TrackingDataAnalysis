function save_diagnostic_track_gifs(caseDef, metrics, outDir, plotOpts, varargin) %#ok<INUSD>

if ~isstruct(metrics) || ~isfield(metrics, 'trackCatalog') || isempty(metrics.trackCatalog)
    warning('Track catalog missing for case %s. Skipping diagnostic GIF export.', char(caseDef.name));
    return;
end

if ~isfolder(outDir), mkdir(outDir); end

trackCatalog = metrics.trackCatalog;
trackRequests = [];
if isfield(caseDef, 'diagnosticTrackIds') && ~isempty(caseDef.diagnosticTrackIds)
    trackRequests = caseDef.diagnosticTrackIds(:).';
end
[selectedIdx, ~] = resolve_diagnostic_track_indices(trackCatalog, trackRequests, caseDef, 'Diagnostic GIF');
if isempty(selectedIdx)
    return;
end

yExtent_mm = plotOpts.inceptionImageSize_px(2) * caseDef.pixelSize;
delayTime = plotOpts.diagnosticGifDelayTime;
trailLength = max(1, round(plotOpts.diagnosticGifTrailLength));

allFrames = collect_selected_frames(trackCatalog, selectedIdx);
if isempty(allFrames)
    warning('No valid frame values found for case %s. Skipping diagnostic GIF export.', char(caseDef.name));
    return;
end

gifPath = fullfile(outDir, sprintf('%s_Re_%g_kDh_%g_all_tracks.gif', char(caseDef.name), caseDef.Re, caseDef.kDh));

[activationXY, activationFrames, activationTrackIds] = activation_arrays(metrics);
[leftFrameAxis, ~, leftFrameCumExposure] = ...
    build_leftmoving_frame_counts(trackCatalog, selectedIdx);
selectedLeftTrackIds = collect_leftmoving_track_ids(trackCatalog, selectedIdx);
if ~isempty(activationTrackIds)
    finiteIds = isfinite(activationTrackIds);
    if any(finiteIds)
        keepAct = ismember(activationTrackIds, selectedLeftTrackIds);
        activationXY = activationXY(keepAct, :);
        activationFrames = activationFrames(keepAct);
        activationTrackIds = activationTrackIds(keepAct);
    end
end

xLim = [0 5];
yLim = [0 1.3];
xTicks = fixed_ticks(xLim(1), xLim(2), 0.5);
yTicks = fixed_ticks(yLim(1), yLim(2), 0.2);

f = figure('Color', 'w', 'Position', [120 120 plotOpts.inceptionImageSize_px], 'Visible', 'off');
ax = axes(f);

for fi = 1:numel(allFrames)
    frameVal = allFrames(fi);
    cla(ax);
    hold(ax, 'on');

    nVisibleAll = 0;
    nVisibleLeftMoving = 0;
    for i = 1:numel(selectedIdx)
        tr = trackCatalog(selectedIdx(i));
        [idxNow, xTail, yTail] = visible_tail_until_frame(tr, frameVal, trailLength);
        if isempty(idxNow)
            continue;
        end

        nVisibleAll = nVisibleAll + 1;
        if is_true_field(tr, 'isLeftMoving')
            nVisibleLeftMoving = nVisibleLeftMoving + 1;
            lineColor = [0 0 1];
            lineWidth = 1.6;
            spotColor = [0 0 0.85];
            markerSize = 16;
        else
            lineColor = [0.55 0.55 0.55];
            lineWidth = 1.0;
            spotColor = [0 0 0];
            markerSize = 12;
        end

        yTailPlot = yExtent_mm - yTail;
        yNowPlot = yExtent_mm - tr.y(idxNow);
        plot(ax, xTail, yTailPlot, '-', 'Color', lineColor, 'LineWidth', lineWidth);
        scatter(ax, tr.x(idxNow), yNowPlot, markerSize, ...
            'Marker', 'o', ...
            'MarkerFaceColor', spotColor, ...
            'MarkerEdgeColor', 'none');
    end

    idxActVisible = isfinite(activationFrames) & (activationFrames <= frameVal);
    if any(idxActVisible)
        xAct = activationXY(idxActVisible, 1);
        yAct = yExtent_mm - activationXY(idxActVisible, 2);
        scatter(ax, xAct, yAct, 24, ...
            'Marker', 'o', ...
            'MarkerFaceColor', [0.9 0.15 0.1], ...
            'MarkerEdgeColor', 'none');
    end

    nActThisFrame = sum(isfinite(activationFrames) & (activationFrames == frameVal));
    cumAct = sum(idxActVisible);
    cumExposure = frame_lookup_arrays(leftFrameAxis, leftFrameCumExposure, frameVal, true);

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
    grid(ax, 'on');
    box(ax, 'on');
    xlabel(ax, '$x\;(\mathrm{mm})$', 'Interpreter', 'latex');
    ylabel(ax, '$y\;(\mathrm{mm})$', 'Interpreter', 'latex');
    title(ax, sprintf(['Diagnostic GIF: %s | frame %g (%d/%d) | visible=%d | ', ...
        'left-moving=%d | act@frame=%d | cum exp=%d | cum act=%d | A/I=%.4g'], ...
        char(caseDef.name), frameVal, fi, numel(allFrames), nVisibleAll, ...
        nVisibleLeftMoving, nActThisFrame, cumExposure, cumAct, metrics.A_over_I));
    apply_plot_theme(ax, 'normal');

    drawnow;
    frame = getframe(f);
    [imInd, map] = rgb2ind(frame2im(frame), 256);
    if fi == 1
        imwrite(imInd, map, gifPath, 'gif', 'LoopCount', Inf, 'DelayTime', delayTime);
    else
        imwrite(imInd, map, gifPath, 'gif', 'WriteMode', 'append', 'DelayTime', delayTime);
    end
end

close(f);
end

function allFrames = collect_selected_frames(trackCatalog, selectedIdx)
allFrames = nan(0,1);
for i = 1:numel(selectedIdx)
    frameVals = trackCatalog(selectedIdx(i)).frame(:);
    frameVals = frameVals(isfinite(frameVals));
    if ~isempty(frameVals)
        allFrames = [allFrames; frameVals(:)]; %#ok<AGROW>
    end
end
allFrames = unique(allFrames);
end

function [activationXY, activationFrames, activationTrackIds] = activation_arrays(metrics)
activationXY = zeros(0,2);
activationFrames = nan(0,1);
activationTrackIds = nan(0,1);

if isfield(metrics, 'activationEvent_xy_netLeftLegacy') && ~isempty(metrics.activationEvent_xy_netLeftLegacy)
    activationXY = metrics.activationEvent_xy_netLeftLegacy;
elseif isfield(metrics, 'activationEvent_xy') && ~isempty(metrics.activationEvent_xy)
    activationXY = metrics.activationEvent_xy;
elseif isfield(metrics, 'strictActivationEvent_xy') && ~isempty(metrics.strictActivationEvent_xy)
    activationXY = metrics.strictActivationEvent_xy;
end

if isfield(metrics, 'activationEvent_frame_netLeftLegacy') && ~isempty(metrics.activationEvent_frame_netLeftLegacy)
    activationFrames = metrics.activationEvent_frame_netLeftLegacy(:);
elseif isfield(metrics, 'activationEvent_frame') && ~isempty(metrics.activationEvent_frame)
    activationFrames = metrics.activationEvent_frame(:);
elseif isfield(metrics, 'strictActivationEvent_frame') && ~isempty(metrics.strictActivationEvent_frame)
    activationFrames = metrics.strictActivationEvent_frame(:);
end

if isfield(metrics, 'activationEvent_trackId_netLeftLegacy') && ~isempty(metrics.activationEvent_trackId_netLeftLegacy)
    activationTrackIds = metrics.activationEvent_trackId_netLeftLegacy(:);
elseif isfield(metrics, 'activationEvent_trackId') && ~isempty(metrics.activationEvent_trackId)
    activationTrackIds = metrics.activationEvent_trackId(:);
elseif isfield(metrics, 'strictActivationEvent_trackId') && ~isempty(metrics.strictActivationEvent_trackId)
    activationTrackIds = metrics.strictActivationEvent_trackId(:);
end

nAct = min(size(activationXY, 1), numel(activationFrames));
if isempty(activationTrackIds)
    activationTrackIds = nan(nAct,1);
else
    nAct = min(nAct, numel(activationTrackIds));
end
if nAct < 1
    activationXY = zeros(0,2);
    activationFrames = nan(0,1);
    activationTrackIds = nan(0,1);
    return;
end
activationXY = activationXY(1:nAct, :);
activationFrames = activationFrames(1:nAct);
activationTrackIds = activationTrackIds(1:nAct);
end

function trackIds = collect_leftmoving_track_ids(trackCatalog, selectedIdx)
trackIds = nan(0,1);
for i = 1:numel(selectedIdx)
    tr = trackCatalog(selectedIdx(i));
    if ~is_true_field(tr, 'isLeftMoving')
        continue;
    end
    if isfield(tr, 'TRACK_ID')
        tid = tr.TRACK_ID;
        if isfinite(tid)
            trackIds(end+1,1) = tid; %#ok<AGROW>
        end
    end
end
trackIds = unique(trackIds);
end

function [frameAxis, nVisible, cumExposure] = build_leftmoving_frame_counts(trackCatalog, selectedIdx)
allLeftFrames = nan(0,1);
for i = 1:numel(selectedIdx)
    tr = trackCatalog(selectedIdx(i));
    if ~is_true_field(tr, 'isLeftMoving')
        continue;
    end
    frameVals = tr.frame(:);
    frameVals = frameVals(isfinite(frameVals));
    if isempty(frameVals)
        continue;
    end
    allLeftFrames = [allLeftFrames; unique(frameVals)]; %#ok<AGROW>
end

frameAxis = unique(allLeftFrames);
if isempty(frameAxis)
    nVisible = nan(0,1);
    cumExposure = nan(0,1);
    return;
end

nFrames = numel(frameAxis);
nVisible = zeros(nFrames,1);
[~, loc] = ismember(allLeftFrames, frameAxis);
loc = loc(loc > 0);
if ~isempty(loc)
    nVisible = accumarray(loc, 1, [nFrames, 1]);
end
cumExposure = cumsum(nVisible);
end

function ticks = fixed_ticks(lo, hi, step)
if ~(isfinite(lo) && isfinite(hi) && isfinite(step) && step > 0)
    ticks = [];
    return;
end
if hi < lo
    tmp = lo;
    lo = hi;
    hi = tmp;
end
ticks = lo:step:hi;
if isempty(ticks) || ticks(end) < hi
    ticks = [ticks, hi]; %#ok<AGROW>
end
end

function [idxNow, xTail, yTail] = visible_tail_until_frame(tr, frameVal, trailLength)
idxNow = [];
xTail = nan(0,1);
yTail = nan(0,1);

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

idxNow = find(f <= frameVal, 1, 'last');
if isempty(idxNow) || idxNow < 1
    idxNow = [];
    return;
end

tailStart = max(1, idxNow - trailLength + 1);
xTail = x(tailStart:idxNow);
yTail = y(tailStart:idxNow);
end

function outVal = frame_lookup_arrays(frameAxis, vals, frameVal, cumulativeMode)
outVal = 0;
if isempty(frameAxis) || isempty(vals)
    return;
end
frameAxis = frameAxis(:);
vals = vals(:);
if numel(frameAxis) ~= numel(vals)
    return;
end

if cumulativeMode
    idx = find(frameAxis <= frameVal, 1, 'last');
else
    idx = find(frameAxis == frameVal, 1, 'first');
end
if isempty(idx) || ~isfinite(vals(idx))
    return;
end
outVal = vals(idx);
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
