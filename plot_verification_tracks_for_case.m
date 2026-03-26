function plot_verification_tracks_for_case(caseDef, metrics, outDir, plotOpts)

if nargin < 4 || ~isfield(plotOpts, 'themes') || isempty(plotOpts.themes)
    plotOpts.themes = "normal";
end
if ~isstruct(metrics) || ~isfield(metrics, 'trackCatalog') || isempty(metrics.trackCatalog)
    warning('Track catalog missing for case %s. Skipping diagnostic track plot.', char(caseDef.name));
    return;
end

trackCatalog = metrics.trackCatalog;
trackRequests = [];
if isfield(caseDef, 'diagnosticTrackIds') && ~isempty(caseDef.diagnosticTrackIds)
    trackRequests = caseDef.diagnosticTrackIds(:).';
end
[selectedIdx, selectedTrackIds] = resolve_diagnostic_track_indices(trackCatalog, trackRequests, caseDef, 'Track diagnostics');
if isempty(selectedIdx)
    warning('Track diagnostics: no tracks selected for case %s.', char(caseDef.name));
    return;
end

[leftIdx, leftTrackIds] = resolve_leftmoving_subset(trackCatalog, selectedIdx, selectedTrackIds, metrics);
[microIdx, microTrackIds] = resolve_microbubble_subset(trackCatalog, selectedIdx, selectedTrackIds, metrics);
if isempty(leftIdx) && isempty(microIdx)
    warning('No left-moving or microbubble-rescue tracks found for case %s. Skipping diagnostic track plot.', char(caseDef.name));
    return;
end

[actXY, actTrackIds] = activation_points_left_and_micro(metrics);
if ~isempty(actTrackIds)
    finiteIds = isfinite(actTrackIds);
    if any(finiteIds)
        keepIds = unique([leftTrackIds(:); microTrackIds(:)]);
        keepAct = ismember(actTrackIds, keepIds);
        actXY = actXY(keepAct, :);
        actTrackIds = actTrackIds(keepAct);
    end
end

yExtent_mm = plotOpts.inceptionImageSize_px(2) * caseDef.pixelSize;
trailLength = max(1, round(plotOpts.diagnosticGifTrailLength));
xLim = [0 5];
yLim = [0 1.2];
if isfield(plotOpts, 'inceptionYLim_mm') && numel(plotOpts.inceptionYLim_mm) >= 2
    yLim = double(plotOpts.inceptionYLim_mm(1:2));
end

for theme = reshape(plotOpts.themes, 1, [])
    f = figure('Color', 'w', 'Position', [120 120 plotOpts.inceptionImageSize_px]);
    ax = axes(f);
    hold(ax, 'on');

    trackColor = [0 0 1];
    activationColor = [0.9 0.15 0.1];

    hTrack = gobjects(0,1);
    for i = 1:numel(leftIdx)
        tr = trackCatalog(leftIdx(i));
        [xComp, yComp] = trail_composite_xy(tr, yExtent_mm, trailLength);
        if isempty(xComp)
            continue;
        end
        hTrack(end+1,1) = plot(ax, xComp, yComp, '-', ...
            'Color', trackColor, ...
            'LineWidth', 1.6); %#ok<AGROW>
    end

    microTrackColor = [0 0.60 0.10];
    hMicro = gobjects(0,1);
    for i = 1:numel(microIdx)
        tr = trackCatalog(microIdx(i));
        [xComp, yComp] = trail_composite_xy(tr, yExtent_mm, trailLength);
        if isempty(xComp)
            continue;
        end
        hMicro(end+1,1) = plot(ax, xComp, yComp, '-', ...
            'Color', microTrackColor, ...
            'LineWidth', 1.6); %#ok<AGROW>
    end

    hAct = gobjects(0,1);
    if ~isempty(actXY)
        yAct = yExtent_mm - actXY(:,2);
        hAct = scatter(ax, actXY(:,1), yAct, 36, ...
            'Marker', 'o', ...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceColor', activationColor);
    end

    xlim(ax, xLim);
    ylim(ax, yLim);
    set(ax, ...
        'XLimMode', 'manual', ...
        'YLimMode', 'manual', ...
        'DataAspectRatioMode', 'auto', ...
        'PlotBoxAspectRatioMode', 'auto');
    xlabel(ax, '$x\;(\mathrm{mm})$', 'Interpreter', 'latex');
    ylabel(ax, '$y\;(\mathrm{mm})$', 'Interpreter', 'latex');

    actTotal = size(actXY, 1);
    nTrackTotal = numel(unique([leftTrackIds(:); microTrackIds(:)]));

    title(ax, sprintf('Track diagnostics | k/Dh %.4g | Re %g | AE %d | tracks %d', ...
        caseDef.kD, caseDef.Re, actTotal, nTrackTotal));
    grid(ax, 'on');
    box(ax, 'on');

    lgdHandles = gobjects(0,1);
    lgdLabels = strings(0,1);
    if ~isempty(hTrack)
        lgdHandles(end+1,1) = hTrack(1); %#ok<AGROW>
        lgdLabels(end+1,1) = "Left moving"; %#ok<AGROW>
    end
    if ~isempty(hMicro)
        lgdHandles(end+1,1) = hMicro(1); %#ok<AGROW>
        lgdLabels(end+1,1) = "Microbubble start 1-120"; %#ok<AGROW>
    end
    if ~isempty(hAct)
        lgdHandles(end+1,1) = hAct; %#ok<AGROW>
        lgdLabels(end+1,1) = "Activation events"; %#ok<AGROW>
    end

    if ~isempty(lgdHandles)
        leg = legend(ax, lgdHandles, cellstr(lgdLabels), 'Location', 'northwest', 'Box', 'off');
    else
        leg = [];
    end

    apply_plot_theme(ax, char(theme));
    style_legend_for_theme(leg, char(theme));

    outBase = fullfile(outDir, sprintf('TrackDiagnostics_%s_Re_%g_kD_%g_%s', char(caseDef.name), caseDef.Re, caseDef.kD, char(theme)));
    save_fig_dual_safe(f, outBase, plotOpts);
    close(f);
end

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
    mask = false(numel(selectedIdx), 1);
    hasStrictFlag = false;
    for i = 1:numel(selectedIdx)
        hasStrictFlag = hasStrictFlag || is_true_field(trackCatalog(selectedIdx(i)), 'isStrictPrimary');
    end
    for i = 1:numel(selectedIdx)
        if hasStrictFlag
            mask(i) = is_true_field(trackCatalog(selectedIdx(i)), 'isStrictPrimary');
        else
            mask(i) = is_true_field(trackCatalog(selectedIdx(i)), 'isLeftMoving');
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

function [actXY, actTrackIds] = activation_points_left_and_micro(metrics)
actXY = zeros(0,2);
actTrackIds = nan(0,1);

[xyA, idA] = get_event_pair(metrics, 'strictActivationEvent_xy', 'strictActivationEvent_trackId');
if isempty(xyA)
    [xyA, idA] = get_event_pair(metrics, 'activationEvent_xy', 'activationEvent_trackId');
end
if isempty(xyA)
    [xyA, idA] = get_event_pair(metrics, 'activationEvent_xy_netLeftLegacy', 'activationEvent_trackId_netLeftLegacy');
end

[xyB, idB] = get_event_pair(metrics, 'microbubbleActivationEvent_xy_nonLeft', 'microbubbleActivationEvent_trackId_nonLeft');

[actXY, actTrackIds] = concat_unique_event_pairs(xyA, idA, xyB, idB);
end

function [xComp, yComp] = trail_composite_xy(tr, yExtent_mm, trailLength)
xComp = nan(0,1);
yComp = nan(0,1);

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
yPlot = yExtent_mm - y;

n = numel(x);
if n == 0
    return;
end

sampleIdx = unique([trailLength:trailLength:n, n]);
for i = 1:numel(sampleIdx)
    idxNow = sampleIdx(i);
    tailStart = max(1, idxNow - trailLength + 1);
    xComp = [xComp; x(tailStart:idxNow); NaN]; %#ok<AGROW>
    yComp = [yComp; yPlot(tailStart:idxNow); NaN]; %#ok<AGROW>
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

function [xy, tid] = get_event_pair(metrics, fieldXY, fieldTid)
xy = zeros(0,2);
tid = nan(0,1);
if isfield(metrics, fieldXY) && ~isempty(metrics.(fieldXY))
    xy = metrics.(fieldXY);
end
if isfield(metrics, fieldTid) && ~isempty(metrics.(fieldTid))
    tid = metrics.(fieldTid)(:);
end
n = min(size(xy,1), numel(tid));
if isempty(tid)
    n = size(xy,1);
    tid = nan(n,1);
end
if n < 1
    xy = zeros(0,2);
    tid = nan(0,1);
    return;
end
xy = xy(1:n, :);
tid = tid(1:n);
end

function [xyOut, tidOut] = concat_unique_event_pairs(xyA, tidA, xyB, tidB)
xy = [xyA; xyB];
tid = [tidA; tidB];

n = min(size(xy,1), numel(tid));
if n < 1
    xyOut = zeros(0,2);
    tidOut = nan(0,1);
    return;
end
xy = xy(1:n, :);
tid = tid(1:n);

keep = true(n,1);
seen = containers.Map('KeyType', 'char', 'ValueType', 'logical');
for i = 1:n
    if isfinite(tid(i))
        key = sprintf('tid:%0.0f', tid(i));
    else
        key = sprintf('xy:%0.6f|%0.6f', xy(i,1), xy(i,2));
    end
    if isKey(seen, key)
        keep(i) = false;
    else
        seen(key) = true;
    end
end

xyOut = xy(keep, :);
tidOut = tid(keep);
end

function style_legend_for_theme(leg, theme)
if isempty(leg) || ~isgraphics(leg)
    return;
end

if strcmp(theme, 'poster')
    leg.TextColor = [1 1 1];
    leg.Color = [0 0 0];
else
    leg.TextColor = [0 0 0];
    leg.Color = 'none';
end
end
