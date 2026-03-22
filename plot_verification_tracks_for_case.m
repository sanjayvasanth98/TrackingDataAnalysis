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

leftMask = false(numel(selectedIdx), 1);
for i = 1:numel(selectedIdx)
    leftMask(i) = is_true_field(trackCatalog(selectedIdx(i)), 'isLeftMoving');
end
leftIdx = selectedIdx(leftMask);
leftTrackIds = selectedTrackIds(leftMask);
if isempty(leftIdx)
    warning('No left-moving tracks found for case %s. Skipping diagnostic track plot.', char(caseDef.name));
    return;
end

[actXY, actTrackIds] = leftmoving_activation_points(metrics);
if ~isempty(actTrackIds)
    finiteIds = isfinite(actTrackIds);
    if any(finiteIds)
        keepAct = ismember(actTrackIds, leftTrackIds);
        actXY = actXY(keepAct, :);
        actTrackIds = actTrackIds(keepAct);
    end
end

yExtent_mm = plotOpts.inceptionImageSize_px(2) * caseDef.pixelSize;
trailLength = max(1, round(plotOpts.diagnosticGifTrailLength));
xLim = [0 5];
yLim = [0 1.3];

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

    actTotal = 0;
    if isfield(metrics, 'activationEventsTotal_netLeftLegacy') && isfinite(metrics.activationEventsTotal_netLeftLegacy)
        actTotal = metrics.activationEventsTotal_netLeftLegacy;
    elseif isfield(metrics, 'activationEventsTotal') && isfinite(metrics.activationEventsTotal)
        actTotal = metrics.activationEventsTotal;
    end

    title(ax, sprintf('Track diagnostics | k/Dh %.4g | Re %g | AE %d | tracks %d', ...
        caseDef.kDh, caseDef.Re, actTotal, numel(leftTrackIds)));
    grid(ax, 'on');
    box(ax, 'on');

    lgdHandles = gobjects(0,1);
    lgdLabels = strings(0,1);
    if ~isempty(hTrack)
        lgdHandles(end+1,1) = hTrack(1); %#ok<AGROW>
        lgdLabels(end+1,1) = "Left moving"; %#ok<AGROW>
    end
    if ~isempty(hAct)
        lgdHandles(end+1,1) = hAct; %#ok<AGROW>
        lgdLabels(end+1,1) = "Activation events"; %#ok<AGROW>
    end

    if ~isempty(lgdHandles)
        leg = legend(ax, lgdHandles, cellstr(lgdLabels), 'Location', 'best', 'Box', 'off');
    else
        leg = [];
    end

    apply_plot_theme(ax, char(theme));
    style_legend_for_theme(leg, char(theme));

    outBase = fullfile(outDir, sprintf('TrackDiagnostics_%s_Re_%g_kDh_%g_%s', char(caseDef.name), caseDef.Re, caseDef.kDh, char(theme)));
    save_fig_dual_safe(f, outBase, plotOpts);
    close(f);
end

end

function [actXY, actTrackIds] = leftmoving_activation_points(metrics)
actXY = zeros(0,2);
actTrackIds = nan(0,1);

if isfield(metrics, 'activationEvent_xy_netLeftLegacy') && ~isempty(metrics.activationEvent_xy_netLeftLegacy)
    actXY = metrics.activationEvent_xy_netLeftLegacy;
elseif isfield(metrics, 'activationEvent_xy') && ~isempty(metrics.activationEvent_xy)
    actXY = metrics.activationEvent_xy;
end

if isfield(metrics, 'activationEvent_trackId_netLeftLegacy') && ~isempty(metrics.activationEvent_trackId_netLeftLegacy)
    actTrackIds = metrics.activationEvent_trackId_netLeftLegacy(:);
elseif isfield(metrics, 'activationEvent_trackId') && ~isempty(metrics.activationEvent_trackId)
    actTrackIds = metrics.activationEvent_trackId(:);
end

nAct = min(size(actXY,1), numel(actTrackIds));
if isempty(actTrackIds)
    nAct = size(actXY,1);
    actTrackIds = nan(nAct,1);
end
if nAct < 1
    actXY = zeros(0,2);
    actTrackIds = nan(0,1);
    return;
end
actXY = actXY(1:nAct, :);
actTrackIds = actTrackIds(1:nAct);
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
