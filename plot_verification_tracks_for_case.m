function plot_verification_tracks_for_case(caseDef, metrics, outDir, plotOpts)

if nargin < 4 || ~isfield(plotOpts, 'themes') || isempty(plotOpts.themes)
    plotOpts.themes = "normal";
end

trackCells = metrics.upstreamTrack_xy;
if isempty(trackCells)
    warning('No upstream tracks found for case %s. Skipping diagnostic track plot.', char(caseDef.name));
    return;
end

yExtent_mm = plotOpts.inceptionImageSize_px(2) * caseDef.pixelSize;

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
    for i = 1:numel(trackCells)
        xy = trackCells{i};
        if isempty(xy)
            continue;
        end

        yTrack = yExtent_mm - xy(:,2);

        hTrack(end+1,1) = plot(ax, xy(:,1), yTrack, '-', ...
            'Color', trackColor, ...
            'LineWidth', 1.0); %#ok<AGROW>

        dx = diff(xy(:,1));
        dy = diff(yTrack);
        nSeg = numel(dx);
        if nSeg >= 1
            arrowIdx = unique(round(linspace(1, nSeg, min(4, nSeg))));
            quiver(ax, xy(arrowIdx,1), yTrack(arrowIdx), dx(arrowIdx), dy(arrowIdx), 0, ...
                'Color', [0 0 0], ...
                'LineWidth', 1.2, ...
                'MaxHeadSize', 2.0, ...
                'AutoScale', 'off');
        end
    end

    hAct = gobjects(0,1);
    if ~isempty(metrics.actLocation_xy)
        yAct = yExtent_mm - metrics.actLocation_xy(:,2);
        hAct = scatter(ax, metrics.actLocation_xy(:,1), yAct, 36, ...
            'Marker', 'o', ...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceColor', activationColor);
    end

    xlim(ax, plotOpts.inceptionXLim_mm);
    ylim(ax, [0 yExtent_mm]);
    axis(ax, 'equal');
    xlabel(ax, '$x\;(\mathrm{mm})$', 'Interpreter', 'latex');
    ylabel(ax, '$y\;(\mathrm{mm})$', 'Interpreter', 'latex');
    title(ax, sprintf('Track diagnostics: %s, Re=%g, k/D_h=%.4g', char(caseDef.name), caseDef.Re, caseDef.kDh));
    grid(ax, 'on');
    box(ax, 'on');

    lgdHandles = gobjects(0,1);
    lgdLabels = strings(0,1);
    if ~isempty(hTrack)
        lgdHandles(end+1,1) = hTrack(1); %#ok<AGROW>
        lgdLabels(end+1,1) = sprintf('Upstream tracks (n=%d)', numel(trackCells)); %#ok<AGROW>
    end
    if ~isempty(hAct)
        lgdHandles(end+1,1) = hAct; %#ok<AGROW>
        lgdLabels(end+1,1) = sprintf('Activated locations (n=%d)', size(metrics.actLocation_xy, 1)); %#ok<AGROW>
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
