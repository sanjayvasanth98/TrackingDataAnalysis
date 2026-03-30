function plot_collapse_rate_vs_frame(allCollapse, figDir, plotOpts)
%PLOT_COLLAPSE_RATE_VS_FRAME  Per-frame collapse-count time series, all cases.
%
%   Thin raw count line + thick smoothed line per case, all on one axes.
%   X-axis: time in milliseconds.  Y-axis: collapses per frame.
%   Smoothing: moving average over ~2.5% of total frames (min 5 frames).

if nargin < 3 || ~isfield(plotOpts,'themes') || isempty(plotOpts.themes)
    plotOpts.themes = "normal";
end

nCases = numel(allCollapse.caseName);
if nCases == 0
    warning('No collapse data to plot.');
    return;
end

cmap = lines(max(nCases,1));

for theme = reshape(plotOpts.themes, 1, [])
    fontName = resolve_plot_font_name();
    f  = figure('Color','w','Position',[100 100 1100 580]);
    ax = axes(f);
    hold(ax,'on');

    lgd    = gobjects(0,1);
    lgdTxt = strings(0,1);

    for ci = 1:nCases
        cd = allCollapse.data{ci};
        if isempty(cd) || cd.nQualified == 0, continue; end

        dt      = allCollapse.dt(ci);
        tAxis   = (double(cd.frameAxis) - double(cd.frameAxis(1))) * dt * 1000;  % ms
        counts  = double(cd.collapseCount);
        col     = cmap(ci,:);

        % Smoothed with moving average
        winLen   = max(5, round(0.025 * numel(counts)));
        kernel   = ones(winLen,1) / winLen;
        smoothed = conv(counts, kernel, 'same');

        % Raw: thin, lightened color
        rawColor = col + (1 - col) * 0.55;
        plot(ax, tAxis, counts, '-', ...
            'Color', rawColor, 'LineWidth', 0.7, ...
            'HandleVisibility','off');

        % Smoothed: thick solid
        h = plot(ax, tAxis, smoothed, '-', ...
            'Color', col, 'LineWidth', 2.2);

        lgd(end+1,1)    = h; %#ok<AGROW>
        lgdTxt(end+1,1) = sprintf('k/d=%.4g  (n=%d, %.3g/s)', ...
            allCollapse.kD(ci), cd.nQualified, cd.ratePerSec); %#ok<AGROW>
    end

    xlabel(ax, 'Time (ms)', 'Interpreter','latex');
    ylabel(ax, 'Collapse events per frame', 'Interpreter','latex');
    title(ax,'');
    grid(ax,'off');
    box(ax,'on');
    set(ax,'FontName',fontName);

    if ~isempty(lgd)
        leg = legend(ax, lgd, cellstr(lgdTxt), ...
            'Location','northoutside','NumColumns',3,'Box','off');
    else
        leg = [];
    end

    apply_plot_theme(ax, char(theme));
    style_legend_for_theme(leg, char(theme));
    % Restore y >= 0 after theme may reset limits
    yl = ylim(ax);
    if yl(1) < 0, ylim(ax, [0 yl(2)]); end

    outBase = fullfile(figDir, "CollapseCount_vs_time_" + theme);
    save_fig_dual_safe(f, outBase, plotOpts);
    if ~isfield(plotOpts,'keepFiguresOpen') || ~plotOpts.keepFiguresOpen
        close(f);
    end
end
fprintf('Saved collapse rate plot to: %s\n', figDir);
end


% =========================================================================
function style_legend_for_theme(leg, theme)
if isempty(leg) || ~isgraphics(leg), return; end
if strcmp(theme,'poster')
    leg.TextColor = [1 1 1];
    leg.Color     = [0 0 0];
else
    leg.TextColor = [0 0 0];
    leg.Color     = 'none';
end
end
