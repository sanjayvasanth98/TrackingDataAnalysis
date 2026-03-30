function plot_breakup_gamma_vs_dratio(breakupData, figDir, plotOpts)
% PLOT_BREAKUP_GAMMA_VS_DRATIO
%   Scatter of gamma = (x_child - x_parent)/d_roughness vs d_child/d_parent.
%   One scatter point per child-parent pair, 50% transparent markers.
%   Per-case coloured mean line (binned mean of gamma vs dRatio).
%   Zero-line reference at gamma=0.
%
%   breakupData: struct array, one element per case, with fields:
%     .caseName  string
%     .kD        double
%     .events    struct array from analyze_breakup_events

if nargin < 3 || ~isfield(plotOpts, 'themes') || isempty(plotOpts.themes)
    plotOpts.themes = "normal";
end

nCases = numel(breakupData);
if nCases == 0
    warning('No breakup data to plot.');
    return;
end

% Pool all dRatio values to determine a sensible x-axis range.
allDRatio = [];
for ci = 1:nCases
    ev = breakupData(ci).events;
    if isempty(ev), continue; end
    allDRatio = [allDRatio; [ev.dRatio].']; %#ok<AGROW>
end
if isempty(allDRatio)
    warning('No breakup events found across all cases.');
    return;
end

xLimPlot = [0, max(allDRatio) * 1.05];

% Per-case colour (using lines colourmap).
cmap = lines(max(nCases, 1));

% Number of bins for the mean line.
nBins = 8;

for theme = reshape(plotOpts.themes, 1, [])
    f = figure('Color', 'w', 'Position', [100 100 900 650]);
    ax = axes(f);
    hold(ax, 'on');

    lgd    = gobjects(0,1);
    lgdTxt = strings(0,1);

    for ci = 1:nCases
        ev = breakupData(ci).events;
        if isempty(ev), continue; end

        dRatio = [ev.dRatio].';
        gamma  = [ev.gamma].';
        col    = cmap(ci, :);

        % --- Scatter (50% transparent) ---
        scatter(ax, dRatio, gamma, 40, col, 'filled', ...
            'MarkerFaceAlpha', 0.50, ...
            'MarkerEdgeColor', 'none', ...
            'HandleVisibility', 'off');

        % --- Binned mean line ---
        xEdges = linspace(xLimPlot(1), xLimPlot(2), nBins + 1);
        xCentre = 0.5 * (xEdges(1:end-1) + xEdges(2:end));
        meanGamma = nan(nBins, 1);
        for b = 1:nBins
            inBin = dRatio >= xEdges(b) & dRatio < xEdges(b+1);
            if sum(inBin) >= 1
                meanGamma(b) = mean(gamma(inBin));
            end
        end
        % Only connect bins that have data.
        hasData = isfinite(meanGamma);
        if sum(hasData) >= 2
            hLine = plot(ax, xCentre(hasData), meanGamma(hasData), '-o', ...
                'Color', col, ...
                'LineWidth', 2.0, ...
                'MarkerSize', 6, ...
                'MarkerFaceColor', col, ...
                'MarkerEdgeColor', 'none');
            lgd(end+1,1)    = hLine; %#ok<AGROW>
            lgdTxt(end+1,1) = sprintf('k/d = %.4g  (n=%d)', ...
                breakupData(ci).kD, numel(dRatio)); %#ok<AGROW>
        end
    end

    % Zero reference line.
    yL = ylim(ax);
    plot(ax, xLimPlot, [0 0], '--', 'Color', [0.5 0.5 0.5], ...
        'LineWidth', 1.2, 'HandleVisibility', 'off');
    ylim(ax, yL);
    xlim(ax, xLimPlot);

    xlabel(ax, '$d_\mathrm{child}/d_\mathrm{parent}$', 'Interpreter', 'latex');
    ylabel(ax, '$\gamma = (x_\mathrm{child} - x_\mathrm{parent})\,/\,d$', ...
        'Interpreter', 'latex');
    title(ax, '');
    grid(ax, 'off');
    box(ax, 'on');

    if ~isempty(lgd)
        leg = legend(ax, lgd, cellstr(lgdTxt), ...
            'Location', 'southoutside', 'NumColumns', 2, 'Box', 'off');
    else
        leg = [];
    end

    apply_plot_theme(ax, char(theme));
    style_legend_for_theme(leg, char(theme));

    outBase = fullfile(figDir, "Breakup_gamma_vs_dRatio_" + theme);
    save_fig_dual_safe(f, outBase, plotOpts);
    if ~isfield(plotOpts, 'keepFiguresOpen') || ~plotOpts.keepFiguresOpen
        close(f);
    end
end

fprintf('Saved breakup plot to: %s\n', figDir);
end


% =========================================================================
function style_legend_for_theme(leg, theme)
if isempty(leg) || ~isgraphics(leg)
    return;
end
if strcmp(theme, 'poster')
    leg.TextColor = [1 1 1];
    leg.Color     = [0 0 0];
else
    leg.TextColor = [0 0 0];
    leg.Color     = 'none';
end
end
