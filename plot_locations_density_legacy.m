function plot_locations_density_legacy(allLoc, outDir, binSize, useAxisLimits, xLim, yLim)

nCases = numel(allLoc.caseName);
if nCases == 0
    warning('No cases in allLoc. Skipping location plots.');
    return;
end

% Gather global XY for auto limits
allXY = [];
for i = 1:nCases
    if ~isempty(allLoc.inj_xy{i}), allXY = [allXY; allLoc.inj_xy{i}]; end %#ok<AGROW>
    if ~isempty(allLoc.act_xy{i}), allXY = [allXY; allLoc.act_xy{i}]; end %#ok<AGROW>
end
if isempty(allXY)
    warning('No location points found. Skipping location plots.');
    return;
end

if useAxisLimits
    xMin = xLim(1); xMax = xLim(2);
    yMin = yLim(1); yMax = yLim(2);
else
    xMin = min(allXY(:,1)); xMax = max(allXY(:,1));
    yMin = min(allXY(:,2)); yMax = max(allXY(:,2));
end

xEdges = xMin:binSize:xMax;
yEdges = yMin:binSize:yMax;

% Two sets
doSet('Injected inception locations (density)', allLoc.inj_xy, 'InjInception', xEdges, yEdges, allLoc, outDir);
doSet('Activated locations (density)',         allLoc.act_xy, 'ActLocations', xEdges, yEdges, allLoc, outDir);

end

function doSet(figTitle, cellXY, fileStem, xE, yE, allLoc, outDir)

nCases = numel(allLoc.caseName);
themes = {'normal','poster'};

for ti = 1:numel(themes)
    theme = themes{ti};

    f = figure('Color','w','Position',[80 80 1500 850]);

    for ci = 1:nCases
        subplot(2,3,ci);

        ax = gca; hold(ax,'on');

        xy = cellXY{ci};
        if isempty(xy)
            text(0.5,0.5,'No points','Units','normalized','HorizontalAlignment','center');
            axis(ax,'off');
            continue;
        end

        % histcounts2 fallback if missing
        N = safe_histcounts2(xy(:,1), xy(:,2), xE, yE);

        imagesc(ax, xE, yE, N.'); % transpose for axis alignment
        set(ax,'YDir','normal');
        axis(ax,'tight');
        axis(ax,'equal');

        xlabel(ax,'x (phys)');
        ylabel(ax,'y (phys)');
        title(ax, sprintf('%s | Sa=%.4g', allLoc.caseName{ci}, allLoc.Sa(ci)));

        grid(ax,'on'); box(ax,'on');
        colormap(ax, parula);
        c = colorbar(ax);
        c.Label.String = 'counts/bin';

        apply_plot_theme(ax, theme);
        if strcmp(theme,'poster')
            c.Color = [1 1 1];
        end
    end

    % Suptitle compatible alternative
    annotation(f,'textbox',[0 0.94 1 0.06], 'String',figTitle, ...
        'EdgeColor','none','HorizontalAlignment','center','FontWeight','bold','FontSize',16);

    outBase = fullfile(outDir, [fileStem '_' theme]);
    save_fig_dual_safe(f, outBase);

    close(f);
end

end

function N = safe_histcounts2(x, y, xEdges, yEdges)
% Use histcounts2 if available, otherwise manual binning.

if exist('histcounts2','file') == 2
    try
        N = histcounts2(x, y, xEdges, yEdges);
        return;
    catch
        % fall through
    end
end

% Manual binning
nx = numel(xEdges)-1;
ny = numel(yEdges)-1;
N = zeros(nx, ny);

% Bin indices
ix = discretize_legacy(x, xEdges);
iy = discretize_legacy(y, yEdges);

good = isfinite(ix) & isfinite(iy) & ix>=1 & ix<=nx & iy>=1 & iy<=ny;
ix = ix(good);
iy = iy(good);

for k = 1:numel(ix)
    N(ix(k), iy(k)) = N(ix(k), iy(k)) + 1;
end
end

function bin = discretize_legacy(v, edges)
% Minimal discretize replacement
bin = nan(size(v));
for i = 1:numel(v)
    vi = v(i);
    if ~isfinite(vi), continue; end
    j = find(vi >= edges(1:end-1) & vi < edges(2:end), 1, 'first');
    if isempty(j) && vi == edges(end)
        j = numel(edges)-1;
    end
    if ~isempty(j), bin(i) = j; end
end
end