function save_diagnostic_track_gifs(caseDef, out, outDir, plotOpts, minNetDx, minTrackSpots, areaJumpFactor)

traj = out.trajectories;
if isempty(traj)
    warning('No trajectories found for case %s. Skipping diagnostic GIF export.', char(caseDef.name));
    return;
end

if ~isfolder(outDir), mkdir(outDir); end

trackRequests = caseDef.diagnosticTrackIds(:).';
selectedIdx = resolve_requested_tracks(traj, trackRequests, caseDef);
if isempty(selectedIdx)
    return;
end

yExtent_mm = plotOpts.inceptionImageSize_px(2) * caseDef.pixelSize;
delayTime = plotOpts.diagnosticGifDelayTime;
trailLength = max(1, round(plotOpts.diagnosticGifTrailLength));

spots = out.spots;
[spotIdSorted, sortIdx] = sort(spots.ID);
spotRowById = containers.Map('KeyType', 'double', 'ValueType', 'double');
for k = 1:numel(spotIdSorted)
    spotRowById(spotIdSorted(k)) = sortIdx(k);
end

trackData = repmat(struct( ...
    'TRACK_ID', NaN, ...
    'frame', [], ...
    'x', [], ...
    'y', [], ...
    'isRecirculating', false, ...
    'isActivated', false), numel(selectedIdx), 1);

allFrames = nan(0,1);
nRecirculating = 0;
nActivated = 0;

for i = 1:numel(selectedIdx)
    thisTraj = traj(selectedIdx(i));
    [isRecirculating, isActivated] = classify_track(thisTraj, spotRowById, spots, minNetDx, minTrackSpots, areaJumpFactor);

    frameVals = thisTraj.frame(:);
    if isempty(frameVals) || ~any(isfinite(frameVals))
        frameVals = (1:numel(thisTraj.x_phys)).';
    end

    trackData(i).TRACK_ID = thisTraj.TRACK_ID;
    trackData(i).frame = frameVals(:);
    trackData(i).x = thisTraj.x_phys(:);
    trackData(i).y = yExtent_mm - thisTraj.y_phys(:);
    trackData(i).isRecirculating = isRecirculating;
    trackData(i).isActivated = isActivated;

    allFrames = [allFrames; frameVals(:)]; %#ok<AGROW>
    nRecirculating = nRecirculating + double(isRecirculating);
    nActivated = nActivated + double(isActivated);
end

allFrames = unique(allFrames(isfinite(allFrames)));
if isempty(allFrames)
    warning('No valid frame values found for case %s. Skipping diagnostic GIF export.', char(caseDef.name));
    return;
end

gifPath = fullfile(outDir, sprintf('%s_Re_%g_kDh_%g_all_tracks.gif', char(caseDef.name), caseDef.Re, caseDef.kDh));

    f = figure('Color', 'w', 'Position', [120 120 plotOpts.inceptionImageSize_px], ...
        'Visible', 'off');
ax = axes(f);

for fi = 1:numel(allFrames)
    frameVal = allFrames(fi);
    cla(ax);
    hold(ax, 'on');

    nVisible = 0;
    for i = 1:numel(trackData)
        idxNow = find(trackData(i).frame <= frameVal, 1, 'last');
        if isempty(idxNow) || idxNow < 1
            continue;
        end

        nVisible = nVisible + 1;
        tailStart = max(1, idxNow - trailLength + 1);
        xTail = trackData(i).x(tailStart:idxNow);
        yTail = trackData(i).y(tailStart:idxNow);

        if trackData(i).isRecirculating
            lineColor = [0 0 1];
            spotColor = [1 0 0];
            lineWidth = 1.6;
            markerSize = 18;
        else
            lineColor = [0.55 0.55 0.55];
            spotColor = [0 0 0];
            lineWidth = 1.0;
            markerSize = 12;
        end

        plot(ax, xTail, yTail, '-', 'Color', lineColor, 'LineWidth', lineWidth);
        scatter(ax, trackData(i).x(idxNow), trackData(i).y(idxNow), markerSize, ...
            'Marker', 'o', ...
            'MarkerFaceColor', spotColor, ...
            'MarkerEdgeColor', 'none');
    end

    xlim(ax, plotOpts.inceptionXLim_mm);
    ylim(ax, [0 yExtent_mm]);
    axis(ax, 'equal');
    grid(ax, 'on');
    box(ax, 'on');
    xlabel(ax, '$x\;(\mathrm{mm})$', 'Interpreter', 'latex');
    ylabel(ax, '$y\;(\mathrm{mm})$', 'Interpreter', 'latex');
    title(ax, sprintf(['Diagnostic GIF: %s | frame %g (%d/%d) | visible=%d | ', ...
        'recirc=%d/%d | activated=%d/%d'], ...
        char(caseDef.name), frameVal, fi, numel(allFrames), nVisible, ...
        nRecirculating, numel(trackData), nActivated, numel(trackData)));
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

function selectedIdx = resolve_requested_tracks(traj, trackRequests, caseDef)
if isempty(trackRequests)
    selectedIdx = 1:numel(traj);
    return;
end

selectedIdx = nan(0,1);
for req = trackRequests
    trackIdx = find([traj.TRACK_ID] == req, 1, 'first');
    if isempty(trackIdx) && isfinite(req) && req >= 1 && req <= numel(traj) && abs(req - round(req)) < 1e-9
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

function [isRecirculating, isActivated] = classify_track(thisTraj, spotRowById, spots, minNetDx, minTrackSpots, areaJumpFactor)
isRecirculating = false;
isActivated = false;

if numel(thisTraj.spotIds) < minTrackSpots || numel(thisTraj.x_phys) < 2
    return;
end

netDx = thisTraj.x_phys(end) - thisTraj.x_phys(1);
isRecirculating = (netDx <= -minNetDx);

areaVals = nan(numel(thisTraj.spotIds),1);
for ii = 1:numel(thisTraj.spotIds)
    sid = thisTraj.spotIds(ii);
    if isKey(spotRowById, sid)
        areaVals(ii) = spots.AREA(spotRowById(sid));
    end
end

a1 = areaVals(1:end-1);
a2 = areaVals(2:end);
ratio = a2 ./ a1;
ratio(~isfinite(ratio)) = NaN;
isActivated = ~isempty(find(ratio >= areaJumpFactor, 1, 'first'));
end
