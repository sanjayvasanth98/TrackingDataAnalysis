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
[selectedIdx, selectedTrackIds] = resolve_diagnostic_track_indices(trackCatalog, trackRequests, caseDef, 'Diagnostic GIF');
if isempty(selectedIdx)
    warning('Diagnostic GIF: no tracks selected for case %s.', char(caseDef.name));
    return;
end

strictMask = false(numel(selectedIdx), 1);
for i = 1:numel(selectedIdx)
    strictMask(i) = is_true_field(trackCatalog(selectedIdx(i)), 'isStrictPrimary');
end
strictIdx = selectedIdx(strictMask);
strictTrackIds = selectedTrackIds(strictMask);

if isempty(strictIdx)
    warning('Diagnostic GIF: selected set has no strict recirculation tracks for case %s.', char(caseDef.name));
    return;
end

[activationXY, activationFrames, activationTrackIds] = strict_activation_events(metrics);
if ~isempty(activationTrackIds)
    keepAct = ismember(activationTrackIds, strictTrackIds);
    activationXY = activationXY(keepAct, :);
    activationFrames = activationFrames(keepAct);
    activationTrackIds = activationTrackIds(keepAct);
end

[frameAxisCounts, nVisibleCounts, nActCounts, cumExpCounts, cumActCounts] = ...
    strict_subset_framewise_counts(trackCatalog, strictIdx, activationFrames);

if isempty(trackRequests)
    verify_strict_parity_warnings(metrics, strictTrackIds, activationTrackIds, caseDef);
end

allFrames = frameAxisCounts(:);
if isempty(allFrames)
    warning('Diagnostic GIF: no valid frame values found for case %s.', char(caseDef.name));
    return;
end

yExtent_mm = plotOpts.inceptionImageSize_px(2) * caseDef.pixelSize;
delayTime = plotOpts.diagnosticGifDelayTime;
trailLength = max(1, round(plotOpts.diagnosticGifTrailLength));

gifPath = fullfile(outDir, sprintf('%s_Re_%g_kDh_%g_strict_tracks.gif', char(caseDef.name), caseDef.Re, caseDef.kDh));

f = figure('Color', 'w', 'Position', [120 120 plotOpts.inceptionImageSize_px], 'Visible', 'off');
ax = axes(f);

for fi = 1:numel(allFrames)
    frameVal = allFrames(fi);
    cla(ax);
    hold(ax, 'on');

    nStrictVisibleNow = 0;
    for i = 1:numel(strictIdx)
        tr = trackCatalog(strictIdx(i));
        [idxNow, xTail, yTail] = visible_tail_until_frame(tr, frameVal, trailLength);
        if isempty(idxNow)
            continue;
        end

        nStrictVisibleNow = nStrictVisibleNow + 1;
        yTailPlot = yExtent_mm - yTail;
        yNowPlot = yExtent_mm - tr.y(idxNow);

        plot(ax, xTail, yTailPlot, '-', 'Color', [0 0 1], 'LineWidth', 1.6);
        scatter(ax, tr.x(idxNow), yNowPlot, 16, ...
            'Marker', 'o', ...
            'MarkerFaceColor', [0 0 0.85], ...
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

    nActThisFrame = frame_lookup_arrays(frameAxisCounts, nActCounts, frameVal, false);
    cumExposure = frame_lookup_arrays(frameAxisCounts, cumExpCounts, frameVal, true);
    cumAct = frame_lookup_arrays(frameAxisCounts, cumActCounts, frameVal, true);

    xlim(ax, plotOpts.inceptionXLim_mm);
    ylim(ax, [0 yExtent_mm]);
    axis(ax, 'equal');
    grid(ax, 'on');
    box(ax, 'on');
    xlabel(ax, '$x\;(\mathrm{mm})$', 'Interpreter', 'latex');
    ylabel(ax, '$y\;(\mathrm{mm})$', 'Interpreter', 'latex');
    title(ax, sprintf(['Diagnostic GIF (strict): %s | frame %g (%d/%d) | visible=%d | ', ...
        'strict-visible=%d | act@frame=%d | cum exp=%d | cum act=%d | A/I=%.4g'], ...
        char(caseDef.name), frameVal, fi, numel(allFrames), nStrictVisibleNow, ...
        frame_lookup_arrays(frameAxisCounts, nVisibleCounts, frameVal, false), ...
        nActThisFrame, cumExposure, cumAct, metrics.A_over_I));
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

function [activationXY, activationFrames, activationTrackIds] = strict_activation_events(metrics)
activationXY = zeros(0,2);
activationFrames = nan(0,1);
activationTrackIds = nan(0,1);

if isfield(metrics, 'strictActivationEvent_xy') && ~isempty(metrics.strictActivationEvent_xy)
    activationXY = metrics.strictActivationEvent_xy;
elseif isfield(metrics, 'activationEvent_xy') && ~isempty(metrics.activationEvent_xy)
    activationXY = metrics.activationEvent_xy;
end

if isfield(metrics, 'strictActivationEvent_frame') && ~isempty(metrics.strictActivationEvent_frame)
    activationFrames = metrics.strictActivationEvent_frame(:);
elseif isfield(metrics, 'activationEvent_frame') && ~isempty(metrics.activationEvent_frame)
    activationFrames = metrics.activationEvent_frame(:);
end

if isfield(metrics, 'strictActivationEvent_trackId') && ~isempty(metrics.strictActivationEvent_trackId)
    activationTrackIds = metrics.strictActivationEvent_trackId(:);
elseif isfield(metrics, 'activationEvent_trackId') && ~isempty(metrics.activationEvent_trackId)
    activationTrackIds = metrics.activationEvent_trackId(:);
end

nAct = min([size(activationXY, 1), numel(activationFrames), numel(activationTrackIds)]);
if isempty(nAct) || ~isfinite(nAct) || nAct < 1
    activationXY = zeros(0,2);
    activationFrames = nan(0,1);
    activationTrackIds = nan(0,1);
    return;
end
activationXY = activationXY(1:nAct, :);
activationFrames = activationFrames(1:nAct);
activationTrackIds = activationTrackIds(1:nAct);
end

function [frameAxis, nVisible, nActEvents, cumExposure, cumActEvents] = ...
    strict_subset_framewise_counts(trackCatalog, strictIdx, activationFrames)

if isempty(strictIdx)
    frameAxis = nan(0,1);
    nVisible = nan(0,1);
    nActEvents = nan(0,1);
    cumExposure = nan(0,1);
    cumActEvents = nan(0,1);
    return;
end

allFrames = nan(0,1);
for i = 1:numel(strictIdx)
    tr = trackCatalog(strictIdx(i));
    frameVals = tr.frame(:);
    frameVals = frameVals(isfinite(frameVals));
    if ~isempty(frameVals)
        allFrames = [allFrames; unique(frameVals)]; %#ok<AGROW>
    end
end

activationFrames = activationFrames(:);
activationFrames = activationFrames(isfinite(activationFrames));
frameAxis = unique([allFrames; activationFrames]);
if isempty(frameAxis)
    nVisible = nan(0,1);
    nActEvents = nan(0,1);
    cumExposure = nan(0,1);
    cumActEvents = nan(0,1);
    return;
end

nFrames = numel(frameAxis);
nVisible = zeros(nFrames,1);
nActEvents = zeros(nFrames,1);

if ~isempty(allFrames)
    [~, locVisible] = ismember(allFrames, frameAxis);
    locVisible = locVisible(locVisible > 0);
    if ~isempty(locVisible)
        nVisible = accumarray(locVisible, 1, [nFrames, 1]);
    end
end

if ~isempty(activationFrames)
    [~, locAct] = ismember(activationFrames, frameAxis);
    locAct = locAct(locAct > 0);
    if ~isempty(locAct)
        nActEvents = accumarray(locAct, 1, [nFrames, 1]);
    end
end

cumExposure = cumsum(nVisible);
cumActEvents = cumsum(nActEvents);
end

function verify_strict_parity_warnings(metrics, strictTrackIds, activationTrackIds, caseDef)
if isfield(metrics, 'nStrictPrimaryTracks') && isfinite(metrics.nStrictPrimaryTracks)
    if numel(strictTrackIds) ~= metrics.nStrictPrimaryTracks
        warning(['Diagnostic GIF parity check failed for %s: strict blue track count (%d) ', ...
            'does not match metrics.nStrictPrimaryTracks (%d).'], ...
            char(caseDef.name), numel(strictTrackIds), metrics.nStrictPrimaryTracks);
    end
end

if isfield(metrics, 'strictActivationEventsTotal') && isfinite(metrics.strictActivationEventsTotal)
    if numel(activationTrackIds) ~= metrics.strictActivationEventsTotal
        warning(['Diagnostic GIF parity check failed for %s: strict activation event count (%d) ', ...
            'does not match metrics.strictActivationEventsTotal (%d).'], ...
            char(caseDef.name), numel(activationTrackIds), metrics.strictActivationEventsTotal);
    end
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
