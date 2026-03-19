function save_diagnostic_track_gifs(caseDef, metrics, outDir, plotOpts, varargin) %#ok<INUSD>

if ~isstruct(metrics) || ~isfield(metrics, 'trackCatalog') || isempty(metrics.trackCatalog)
    warning('Track catalog missing for case %s. Skipping diagnostic GIF export.', char(caseDef.name));
    return;
end

if ~isfolder(outDir), mkdir(outDir); end

trackCatalog = metrics.trackCatalog;
trackRequests = caseDef.diagnosticTrackIds(:).';
selectedIdx = resolve_requested_tracks(trackCatalog, trackRequests, caseDef);
if isempty(selectedIdx)
    return;
end

yExtent_mm = plotOpts.inceptionImageSize_px(2) * caseDef.pixelSize;
delayTime = plotOpts.diagnosticGifDelayTime;
trailLength = max(1, round(plotOpts.diagnosticGifTrailLength));

allFrames = nan(0,1);
for i = 1:numel(selectedIdx)
    frameVals = trackCatalog(selectedIdx(i)).frame(:);
    frameVals = frameVals(isfinite(frameVals));
    if ~isempty(frameVals)
        allFrames = [allFrames; frameVals(:)]; %#ok<AGROW>
    end
end
allFrames = unique(allFrames);
if isempty(allFrames)
    warning('No valid frame values found for case %s. Skipping diagnostic GIF export.', char(caseDef.name));
    return;
end

gifPath = fullfile(outDir, sprintf('%s_Re_%g_kDh_%g_all_tracks.gif', char(caseDef.name), caseDef.Re, caseDef.kDh));

activationXY = zeros(0,2);
activationFrames = nan(0,1);
if isfield(metrics, 'activationEvent_xy') && ~isempty(metrics.activationEvent_xy)
    activationXY = metrics.activationEvent_xy;
end
if isfield(metrics, 'activationEvent_frame') && ~isempty(metrics.activationEvent_frame)
    activationFrames = metrics.activationEvent_frame(:);
end
if ~isempty(activationXY) || ~isempty(activationFrames)
    nAct = min(size(activationXY, 1), numel(activationFrames));
    activationXY = activationXY(1:nAct, :);
    activationFrames = activationFrames(1:nAct);
end

f = figure('Color', 'w', 'Position', [120 120 plotOpts.inceptionImageSize_px], 'Visible', 'off');
ax = axes(f);

for fi = 1:numel(allFrames)
    frameVal = allFrames(fi);
    cla(ax);
    hold(ax, 'on');

    nVisibleAll = 0;
    for i = 1:numel(selectedIdx)
        tr = trackCatalog(selectedIdx(i));
        [idxNow, xTail, yTail] = visible_tail_until_frame(tr, frameVal, trailLength);
        if isempty(idxNow)
            continue;
        end

        nVisibleAll = nVisibleAll + 1;

        if tr.isLeftMoving
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

    nLeftVisible = frame_lookup(metrics, 'frame_nLeftMovingVisible', frameVal, false);
    nActThisFrame = frame_lookup(metrics, 'frame_nActivationEvents', frameVal, false);
    cumExposure = frame_lookup(metrics, 'frame_cumExposure', frameVal, true);
    cumAct = frame_lookup(metrics, 'frame_cumActivationEvents', frameVal, true);

    xlim(ax, plotOpts.inceptionXLim_mm);
    ylim(ax, [0 yExtent_mm]);
    axis(ax, 'equal');
    grid(ax, 'on');
    box(ax, 'on');
    xlabel(ax, '$x\;(\mathrm{mm})$', 'Interpreter', 'latex');
    ylabel(ax, '$y\;(\mathrm{mm})$', 'Interpreter', 'latex');
    title(ax, sprintf(['Diagnostic GIF: %s | frame %g (%d/%d) | visible=%d | ', ...
        'left-moving=%d | act@frame=%d | cum exp=%d | cum act=%d | A/I=%.4g'], ...
        char(caseDef.name), frameVal, fi, numel(allFrames), nVisibleAll, ...
        nLeftVisible, nActThisFrame, cumExposure, cumAct, metrics.A_over_I));
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

function selectedIdx = resolve_requested_tracks(trackCatalog, trackRequests, caseDef)
if isempty(trackRequests)
    selectedIdx = (1:numel(trackCatalog)).';
    return;
end

trackIds = [trackCatalog.TRACK_ID];
selectedIdx = nan(0,1);
for req = trackRequests
    trackIdx = find(trackIds == req, 1, 'first');
    if isempty(trackIdx) && isfinite(req) && req >= 1 && req <= numel(trackCatalog) && abs(req - round(req)) < 1e-9
        trackIdx = round(req);
    end

    if isempty(trackIdx)
        warning('Diagnostic GIF track request %g not found for case %s.', req, char(caseDef.name));
        continue;
    end

    selectedIdx(end+1,1) = trackIdx; %#ok<AGROW>
end

selectedIdx = unique(selectedIdx, 'stable');
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

function outVal = frame_lookup(metrics, fieldName, frameVal, cumulativeMode)
outVal = 0;
if ~isfield(metrics, 'frame_axis') || isempty(metrics.frame_axis)
    return;
end
if ~isfield(metrics, fieldName) || isempty(metrics.(fieldName))
    return;
end

frameAxis = metrics.frame_axis(:);
vals = metrics.(fieldName)(:);
if numel(vals) ~= numel(frameAxis)
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
