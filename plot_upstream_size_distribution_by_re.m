function plot_upstream_size_distribution_by_re(allSize, outDir, ~, plotOpts)

if nargin < 4 || ~isfield(plotOpts, 'themes') || isempty(plotOpts.themes)
    plotOpts.themes = "normal";
end

nCases = numel(allSize.caseName);
if nCases == 0
    warning('No size data found. Skipping size distribution plots.');
    return;
end

ReVals = unique(allSize.Re(:));
if isempty(ReVals)
    warning('No Reynolds values found in size data.');
    return;
end

for theme = reshape(plotOpts.themes, 1, [])
    for r = 1:numel(ReVals)
        Rei = ReVals(r);

        idxRe = find(allSize.Re == Rei);
        if isempty(idxRe)
            continue;
        end

        if ~isempty(plotOpts.upstreamSizeXLim_um)
            xLimUse = plotOpts.upstreamSizeXLim_um;
        else
            pooledX = [];
            for j = 1:numel(idxRe)
                xCase = allSize.size_eqd{idxRe(j)};
                xCase = xCase(isfinite(xCase) & xCase > 0) * 1000; %#ok<AGROW>
                if ~isempty(xCase)
                    pooledX = [pooledX; xCase(:)]; %#ok<AGROW>
                end
            end
            xLimUse = robust_size_xlim(pooledX);
        end

        f = figure('Color', 'w', 'Position', [120 120 1000 700]);
        ax = axes(f);
        hold(ax, 'on');

        cmap = lines(max(numel(idxRe),1));
        lgd = gobjects(0,1);
        lgdTxt = strings(0,1);
        yMaxPlot = 0;

        for j = 1:numel(idxRe)
            ci = idxRe(j);
            x = allSize.size_eqd{ci};
            x = x(isfinite(x) & x > 0) * 1000; % mm -> um

            if numel(x) < 5
                continue;
            end

            xi = linspace(xLimUse(1), xLimUse(2), 300);
            fhat = ksdensity(x, xi);

            hFill = fill(ax, [xi, fliplr(xi)], [fhat, zeros(size(fhat))], cmap(j,:), ...
                'FaceAlpha', 0.5, ...
                'EdgeColor', 'none');
            plot(ax, xi, fhat, '-', 'LineWidth', 2.0, 'Color', cmap(j,:));
            yMaxPlot = max(yMaxPlot, max(fhat));

            lgd(end+1,1) = hFill; %#ok<AGROW>
            lgdTxt(end+1,1) = sprintf('%s, k/D_h=%.4g', allSize.caseName(ci), allSize.kDh(ci)); %#ok<AGROW>
        end

        xlim(ax, xLimUse);
        if yMaxPlot > 0
            ylim(ax, [0, yMaxPlot * 1.10]);
        end
        xlabel(ax, 'Equivalent diameter, $d_{eq}\,(\mu\mathrm{m})$', 'Interpreter', 'latex');
        ylabel(ax, 'PDF', 'Interpreter', 'latex');
        title(ax, sprintf('Upstream-moving microbubble size distribution, Re=%g', Rei));
        grid(ax, 'on');
        box(ax, 'on');

        if ~isempty(lgd)
            leg = legend(ax, lgd, cellstr(lgdTxt), 'Location', 'eastoutside', 'Box', 'off');
        else
            leg = [];
        end

        apply_plot_theme(ax, char(theme));
        style_legend_for_theme(leg, char(theme));

        outBase = fullfile(outDir, sprintf('UpstreamSizeDist_Re_%g_%s', Rei, theme));
        save_fig_dual_safe(f, outBase, plotOpts);
        close(f);
    end
end

end

function xLim = robust_size_xlim(x)
if isempty(x)
    xLim = [0 30];
    return;
end

x = x(isfinite(x));
if isempty(x)
    xLim = [0 30];
    return;
end

xMin = min(x);
xHi = prctile(x, 99);

if ~(isfinite(xHi) && xHi > xMin)
    xHi = max(x);
end
if ~(isfinite(xHi) && xHi > xMin)
    xHi = xMin + max(abs(xMin) * 0.1, 1);
end

pad = max(0.05 * (xHi - xMin), 5);
xLo = max(0, xMin - pad);
xLim = [xLo, xHi + pad];
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
