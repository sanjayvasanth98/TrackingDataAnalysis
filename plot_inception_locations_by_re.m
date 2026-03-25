function plot_inception_locations_by_re(allLoc, outDir, plotOpts)

if nargin < 3 || ~isfield(plotOpts, 'themes') || isempty(plotOpts.themes)
    plotOpts.themes = "normal";
end

nCases = numel(allLoc.caseName);
if nCases == 0
    warning('No cases in allLoc. Skipping inception location plot.');
    return;
end

ReVals = unique(allLoc.Re(:));
if isempty(ReVals)
    warning('No Reynolds numbers found. Skipping inception location plot.');
    return;
end

% Shared calibration assumption: all cases use the same pixel size.
pixelSizeVals = allLoc.pixelSize(isfinite(allLoc.pixelSize) & allLoc.pixelSize > 0);
if isempty(pixelSizeVals)
    warning('No valid pixel sizes found. Skipping inception location plot.');
    return;
end

yExtent_mm = plotOpts.inceptionImageSize_px(2) * median(pixelSizeVals);
xLim = [0 5];
if isfield(plotOpts, 'inceptionXLim_mm') && numel(plotOpts.inceptionXLim_mm) >= 2
    xLim = double(plotOpts.inceptionXLim_mm(1:2));
end
yLim = [0 1.2];
if isfield(plotOpts, 'inceptionYLim_mm') && numel(plotOpts.inceptionYLim_mm) >= 2
    yLim = double(plotOpts.inceptionYLim_mm(1:2));
end

hasPoints = false;
for i = 1:nCases
    if ~isempty(allLoc.inception2x_xy{i})
        hasPoints = true;
        break;
    end
end
if ~hasPoints
    warning('No activation points found on left-moving tracks.');
    return;
end

for theme = reshape(plotOpts.themes, 1, [])
    for r = 1:numel(ReVals)
        Rei = ReVals(r);

        f = figure('Color', 'w', 'Position', [120 120 plotOpts.inceptionImageSize_px]);
        ax = axes(f);
        hold(ax, 'on');

        idxRe = find(allLoc.Re == Rei);
        nReCases = numel(idxRe);
        cmap = lines(max(nReCases, 1));

        lgd = gobjects(0,1);
        lgdTxt = strings(0,1);

        for j = 1:nReCases
            ci = idxRe(j);
            xy = allLoc.inception2x_xy{ci};
            if isempty(xy)
                continue;
            end

            yPlot = yExtent_mm - xy(:,2);

            h = scatter(ax, xy(:,1), yPlot, 34, ...
                'Marker', 'o', ...
                'MarkerEdgeColor', 'none', ...
                'MarkerFaceColor', cmap(j,:));

            lgd(end+1,1) = h; %#ok<AGROW>
            lgdTxt(end+1,1) = sprintf('%s, k/D_h=%.4g', allLoc.caseName(ci), allLoc.kDh(ci)); %#ok<AGROW>
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
        title(ax, sprintf('Activation locations on left-moving tracks, Re=%g', Rei));
        grid(ax, 'on');
        box(ax, 'on');

        if ~isempty(lgd)
            leg = legend(ax, lgd, cellstr(lgdTxt), 'Location', 'eastoutside', 'Box', 'off');
        else
            leg = [];
        end

        apply_plot_theme(ax, char(theme));
        style_legend_for_theme(leg, char(theme));

        outBase = fullfile(outDir, sprintf('Inception2x_locations_Re_%g_%s', Rei, theme));
        save_fig_dual_safe(f, outBase, plotOpts);
        close(f);
    end
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
