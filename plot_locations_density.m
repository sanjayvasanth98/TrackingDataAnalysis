function plot_locations_density(allLoc, outDir, binSize, useAxisLimits, xLim, yLim)

nCases = numel(allLoc.caseName);

% Helper to build bins per case consistently
% Determine global limits unless user overrides
allXY = [];
for i = 1:nCases
    if ~isempty(allLoc.inj_xy{i}), allXY = [allXY; allLoc.inj_xy{i}]; end %#ok<AGROW>
    if ~isempty(allLoc.act_xy{i}), allXY = [allXY; allLoc.act_xy{i}]; end %#ok<AGROW>
end

if isempty(allXY)
    warning("No location points found. Skipping location plots.");
    return;
end

if useAxisLimits
    xMin = xLim(1); xMax = xLim(2);
    yMin = yLim(1); yMax = yLim(2);
else
    xMin = min(allXY(:,1)); xMax = max(allXY(:,1));
    yMin = min(allXY(:,2)); yMax = max(allXY(:,2));
end

% Bin edges
xEdges = xMin:binSize:xMax;
yEdges = yMin:binSize:yMax;

% Two plot sets: injected inception, activated locations
plotOneSet("Injected inception locations (density)", allLoc.inj_xy, "InjInception", xEdges, yEdges);
plotOneSet("Activated locations (density)",         allLoc.act_xy, "ActLocations", xEdges, yEdges);

    function plotOneSet(figTitle, cellXY, fileStem, xE, yE)

        for theme = ["normal","poster"]
            f = figure('Color','w','Position',[80 80 1400 800]);

            t = tiledlayout(f, 2, 3, 'TileSpacing','compact', 'Padding','compact');
            title(t, figTitle);

            for ci = 1:nCases
                nexttile(t);
                ax = gca; hold(ax,'on');

                xy = cellXY{ci};
                if isempty(xy)
                    text(ax, 0.5, 0.5, "No points", 'HorizontalAlignment','center');
                    axis(ax,'off');
                    continue;
                end

                % 2D counts
                N = histcounts2(xy(:,1), xy(:,2), xE, yE);

                % Calculate bin centers for imagesc
                xC = xE(1:end-1) + diff(xE)/2;
                yC = yE(1:end-1) + diff(yE)/2;

                % Display as image (density)
                imagesc(ax, xC, yC, N.'); % transpose so axes align
                set(ax,'YDir','normal');
                axis(ax,'tight');
                axis(ax,'equal');

                xlabel(ax, 'x (phys)');
                ylabel(ax, 'y (phys)');
                title(ax, sprintf('%s | Sa=%.4g', allLoc.caseName(ci), allLoc.Sa(ci)));

                grid(ax,'on'); box(ax,'on');
                colormap(ax, parula); % OK for density; if you want specific, tell me
                c = colorbar(ax);
                c.Label.String = "counts/bin";

                apply_plot_theme(ax, theme);

                % For poster theme, colorbar text/edge also needs to be white-ish
                if theme == "poster"
                    c.Color = [1 1 1];
                end
            end

            outBase = fullfile(outDir, fileStem + "_" + theme);
            save_fig_dual(f, outBase);
            close(f);
        end

    end

end