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

strictMask = false(numel(selectedIdx), 1);
for i = 1:numel(selectedIdx)
    strictMask(i) = is_true_field(trackCatalog(selectedIdx(i)), 'isStrictPrimary');
end
strictIdx = selectedIdx(strictMask);
strictTrackIds = selectedTrackIds(strictMask);
if isempty(strictIdx)
    warning('No strict recirculation tracks found for case %s. Skipping diagnostic track plot.', char(caseDef.name));
    return;
end

[actXY, actTrackIds] = strict_activation_points(metrics);
if ~isempty(actTrackIds)
    keepAct = ismember(actTrackIds, strictTrackIds);
    actXY = actXY(keepAct, :);
    actTrackIds = actTrackIds(keepAct);
end

if isempty(trackRequests)
    verify_static_parity_warnings(metrics, strictTrackIds, actTrackIds, caseDef);
end

yExtent_mm = plotOpts.inceptionImageSize_px(2) * caseDef.pixelSize;
trailLength = max(1, round(plotOpts.diagnosticGifTrailLength));

for theme = reshape(plotOpts.themes, 1, [])
    f = figure('Color', 'w', 'Position', [120 120 plotOpts.inceptionImageSize_px]);
    ax = axes(f);
    hold(ax, 'on');

    if strcmp(theme, 'poster')
        trackColor = [0.35 0.75 1.00];
        activationColor = [1.00 0.45 0.20];
    else
        trackColor = [0.10 0.35 0.80];
        activationColor = [0.85 0.20 0.15];
    end

    hTrack = gobjects(0,1);
    for i = 1:numel(strictIdx)
        tr = trackCatalog(strictIdx(i));
        [xComp, yComp] = trail_composite_xy(tr, yExtent_mm, trailLength);
        if isempty(xComp)
            continue;
        end
        hTrack(end+1,1) = plot(ax, xComp, yComp, '-', ...
            'Color', trackColor, ...
            'LineWidth', 1.0); %#ok<AGROW>
    end

    hAct = gobjects(0,1);
    if ~isempty(actXY)
        yAct = yExtent_mm - actXY(:,2);
        hAct = scatter(ax, actXY(:,1), yAct, 36, ...
            'Marker', 'o', ...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceColor', activationColor);
    end

    xlim(ax, plotOpts.inceptionXLim_mm);
    ylim(ax, [0 yExtent_mm]);
    axis(ax, 'equal');
    xlabel(ax, '$x\;(\mathrm{mm})$', 'Interpreter', 'latex');
    ylabel(ax, '$y\;(\mathrm{mm})$', 'Interpreter', 'latex');

    exposureVal = 0;
    if isfield(metrics, 'strictTrackFrameExposure') && isfinite(metrics.strictTrackFrameExposure)
        exposureVal = metrics.strictTrackFrameExposure;
    elseif isfield(metrics, 'leftMovingTrackFrameExposure') && isfinite(metrics.leftMovingTrackFrameExposure)
        exposureVal = metrics.leftMovingTrackFrameExposure;
    end

    actTotal = 0;
    if isfield(metrics, 'strictActivationEventsTotal') && isfinite(metrics.strictActivationEventsTotal)
        actTotal = metrics.strictActivationEventsTotal;
    elseif isfield(metrics, 'activationEventsTotal') && isfinite(metrics.activationEventsTotal)
        actTotal = metrics.activationEventsTotal;
    end

    title(ax, sprintf(['Track diagnostics (strict recirculation): %s, Re=%g, k/D_h=%.4g | ', ...
        'A/I=%.4g, events=%d, exposure=%d'], ...
        char(caseDef.name), caseDef.Re, caseDef.kDh, metrics.A_over_I, actTotal, exposureVal));
    grid(ax, 'on');
    box(ax, 'on');

    lgdHandles = gobjects(0,1);
    lgdLabels = strings(0,1);
    if ~isempty(hTrack)
        lgdHandles(end+1,1) = hTrack(1); %#ok<AGROW>
        lgdLabels(end+1,1) = sprintf('Strict recirculation tracks (A/I denominator, n=%d)', numel(strictTrackIds)); %#ok<AGROW>
    end
    if ~isempty(hAct)
        lgdHandles(end+1,1) = hAct; %#ok<AGROW>
        lgdLabels(end+1,1) = sprintf('Strict activation events (A/I numerator, n=%d)', size(actXY, 1)); %#ok<AGROW>
    end

    if ~isempty(lgdHandles)
        leg = legend(ax, lgdHandles, cellstr(lgdLabels), 'Location', 'eastoutside', 'Box', 'off');
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

function [actXY, actTrackIds] = strict_activation_points(metrics)
actXY = zeros(0,2);
actTrackIds = nan(0,1);

if isfield(metrics, 'strictActivationEvent_xy') && ~isempty(metrics.strictActivationEvent_xy)
    actXY = metrics.strictActivationEvent_xy;
elseif isfield(metrics, 'activationEvent_xy') && ~isempty(metrics.activationEvent_xy)
    actXY = metrics.activationEvent_xy;
end

if isfield(metrics, 'strictActivationEvent_trackId') && ~isempty(metrics.strictActivationEvent_trackId)
    actTrackIds = metrics.strictActivationEvent_trackId(:);
elseif isfield(metrics, 'activationEvent_trackId') && ~isempty(metrics.activationEvent_trackId)
    actTrackIds = metrics.activationEvent_trackId(:);
end

nAct = min(size(actXY,1), numel(actTrackIds));
if nAct < 1
    actXY = zeros(0,2);
    actTrackIds = nan(0,1);
    return;
end
actXY = actXY(1:nAct, :);
actTrackIds = actTrackIds(1:nAct);
end

function verify_static_parity_warnings(metrics, strictTrackIds, actTrackIds, caseDef)
if isfield(metrics, 'nStrictPrimaryTracks') && isfinite(metrics.nStrictPrimaryTracks)
    if numel(strictTrackIds) ~= metrics.nStrictPrimaryTracks
        warning(['Track diagnostics parity check failed for %s: strict blue track count (%d) ', ...
            'does not match metrics.nStrictPrimaryTracks (%d).'], ...
            char(caseDef.name), numel(strictTrackIds), metrics.nStrictPrimaryTracks);
    end
end

if isfield(metrics, 'strictActivationEventsTotal') && isfinite(metrics.strictActivationEventsTotal)
    if numel(actTrackIds) ~= metrics.strictActivationEventsTotal
        warning(['Track diagnostics parity check failed for %s: strict activation event count (%d) ', ...
            'does not match metrics.strictActivationEventsTotal (%d).'], ...
            char(caseDef.name), numel(actTrackIds), metrics.strictActivationEventsTotal);
    end
end
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
