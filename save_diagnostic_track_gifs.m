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

maxFrames = Inf;
if isfield(plotOpts, 'diagnosticGifMaxFrames') && ~isempty(plotOpts.diagnosticGifMaxFrames)
    maxFrames = double(plotOpts.diagnosticGifMaxFrames);
end
if isfinite(maxFrames) && maxFrames > 0 && numel(allFrames) > maxFrames
    allFrames = allFrames(1:maxFrames);
end

gifPath = fullfile(outDir, sprintf('%s_Re_%g_kD_%g_all_tracks.gif', char(caseDef.name), caseDef.Re, caseDef.kD));

[activationXY, activationFrames, activationTrackIds] = activation_arrays(metrics);
[leftMaskSelected, selectedLeftTrackIds] = resolve_leftmoving_selection(metrics, trackCatalog, selectedIdx);
[microMaskSelected, selectedMicroTrackIds] = resolve_microbubble_selection(metrics, trackCatalog, selectedIdx);
if ~any(leftMaskSelected) && ~any(microMaskSelected)
    warning('No left-moving or microbubble-rescue tracks selected for case %s. Skipping diagnostic GIF export.', char(caseDef.name));
    return;
end
microBlueState = false(numel(selectedIdx), 1);
if ~isempty(activationTrackIds)
    finiteIds = isfinite(activationTrackIds);
    if any(finiteIds)
        keepIds = unique([selectedLeftTrackIds(:); selectedMicroTrackIds(:)]);
        keepAct = ismember(activationTrackIds, keepIds);
        activationXY = activationXY(keepAct, :);
        activationFrames = activationFrames(keepAct);
        activationTrackIds = activationTrackIds(keepAct);
    end
end

xLim = [0 4.8];
yLim = [0 1.2];
if isfield(plotOpts, 'inceptionYLim_mm') && numel(plotOpts.inceptionYLim_mm) >= 2
    yLim = double(plotOpts.inceptionYLim_mm(1:2));
end
xTicks = fixed_ticks(xLim(1), xLim(2), 0.5);
yTicks = fixed_ticks(yLim(1), yLim(2), 0.2);

f = figure('Color', 'w', 'Position', [120 120 plotOpts.inceptionImageSize_px], 'Visible', 'off');
ax = axes(f);

for fi = 1:numel(allFrames)
    frameVal = allFrames(fi);
    cla(ax);
    hold(ax, 'on');

    minSpotsForDisplay = 5;
    nVisibleLeftMoving = 0;
    nVisibleMicroRescue = 0;
    for i = 1:numel(selectedIdx)
        tr = trackCatalog(selectedIdx(i));
        if numel(tr.frame) < minSpotsForDisplay, continue; end
        [idxNow, xTail, yTail] = visible_tail_until_frame(tr, frameVal, trailLength);
        if isempty(idxNow)
            continue;
        end

        if leftMaskSelected(i)
            nVisibleLeftMoving = nVisibleLeftMoving + 1;
            lineColor = [0 0 1];
            lineWidth = 1.6;
            spotColor = [0 0 0.85];
            markerSize = 16;
        elseif microMaskSelected(i)
            nVisibleMicroRescue = nVisibleMicroRescue + 1;
            motionSignal = recent_motion_signal(xTail, 3);
            if motionSignal < 0
                microBlueState(i) = true;
            elseif motionSignal > 0
                microBlueState(i) = false;
            end

            if microBlueState(i)
                lineColor = [0 0 1];
                lineWidth = 1.6;
                spotColor = [0 0 0.85];
                markerSize = 16;
            else
                lineColor = [0 0.60 0.10];
                lineWidth = 1.6;
                spotColor = [0 0.55 0.10];
                markerSize = 16;
            end
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

    % Legend proxy handles
    hBlue  = plot(ax, NaN, NaN, '-', 'Color', [0 0 1], 'LineWidth', 1.6);
    hGreen = plot(ax, NaN, NaN, '-', 'Color', [0 0.60 0.10], 'LineWidth', 1.6);
    hRed   = plot(ax, NaN, NaN, 'o', 'MarkerFaceColor', [0.9 0.15 0.1], 'MarkerEdgeColor', 'none', 'MarkerSize', 6);
    hGray  = plot(ax, NaN, NaN, '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.0);
    legend(ax, [hBlue, hGreen, hRed, hGray], ...
        {'Upstream tracks', 'Microbubbles (1-120\mum)', 'Activation event', 'Other tracks'}, ...
        'Location', 'northwest', 'Box', 'off', ...
        'FontName', 'Times New Roman', 'FontSize', 8);

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
    grid(ax, 'off');
    box(ax, 'on');
    xlabel(ax, '$x\;(\mathrm{mm})$', 'Interpreter', 'latex');
    ylabel(ax, '$y\;(\mathrm{mm})$', 'Interpreter', 'latex');
    title(ax, sprintf(['Diagnostic GIF: k/d=%g | frame %g (%d/%d) | ', ...
        'left-moving=%d | micro-rescue=%d | act@frame=%d | AE %d'], ...
        caseDef.kD, frameVal, fi, numel(allFrames), ...
        nVisibleLeftMoving, nVisibleMicroRescue, nActThisFrame, numel(activationFrames)));
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

[activationXY, activationFrames, activationTrackIds] = concat_unique_events( ...
    xyA, frA, idA, xyB, frB, idB);
end

function [leftMaskSelected, leftTrackIds] = resolve_leftmoving_selection(metrics, trackCatalog, selectedIdx)
leftMaskSelected = false(numel(selectedIdx), 1);
leftTrackIds = nan(0,1);

metricIds = nan(0,1);
if isfield(metrics, 'strictPrimaryTrackIds') && ~isempty(metrics.strictPrimaryTrackIds)
    metricIds = metrics.strictPrimaryTrackIds(:);
elseif isfield(metrics, 'leftMovingTrackIds_netLeftLegacy') && ~isempty(metrics.leftMovingTrackIds_netLeftLegacy)
    metricIds = metrics.leftMovingTrackIds_netLeftLegacy(:);
end
metricIds = unique(metricIds(isfinite(metricIds)), 'stable');

if ~isempty(metricIds)
    for i = 1:numel(selectedIdx)
        tid = get_track_id(trackCatalog(selectedIdx(i)));
        if isfinite(tid) && ismember(tid, metricIds)
            leftMaskSelected(i) = true;
        end
    end
else
    hasStrictFlag = false;
    for i = 1:numel(selectedIdx)
        hasStrictFlag = hasStrictFlag || is_true_field(trackCatalog(selectedIdx(i)), 'isStrictPrimary');
    end
    for i = 1:numel(selectedIdx)
        if hasStrictFlag
            leftMaskSelected(i) = is_true_field(trackCatalog(selectedIdx(i)), 'isStrictPrimary');
        else
            leftMaskSelected(i) = is_true_field(trackCatalog(selectedIdx(i)), 'isLeftMoving');
        end
    end
end

for i = 1:numel(selectedIdx)
    if ~leftMaskSelected(i)
        continue;
    end
    tid = get_track_id(trackCatalog(selectedIdx(i)));
    if isfinite(tid)
        leftTrackIds(end+1,1) = tid; %#ok<AGROW>
    end
end
leftTrackIds = unique(leftTrackIds, 'stable');
end

function [microMaskSelected, microTrackIds] = resolve_microbubble_selection(metrics, trackCatalog, selectedIdx)
microMaskSelected = false(numel(selectedIdx), 1);
microTrackIds = nan(0,1);

metricIds = nan(0,1);
if isfield(metrics, 'microbubbleRescueTrackIds_nonLeft') && ~isempty(metrics.microbubbleRescueTrackIds_nonLeft)
    metricIds = metrics.microbubbleRescueTrackIds_nonLeft(:);
end
metricIds = unique(metricIds(isfinite(metricIds)), 'stable');

if isempty(metricIds)
    return;
end

for i = 1:numel(selectedIdx)
    tid = get_track_id(trackCatalog(selectedIdx(i)));
    if isfinite(tid) && ismember(tid, metricIds)
        microMaskSelected(i) = true;
        microTrackIds(end+1,1) = tid; %#ok<AGROW>
    end
end
microTrackIds = unique(microTrackIds(isfinite(microTrackIds)), 'stable');
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
