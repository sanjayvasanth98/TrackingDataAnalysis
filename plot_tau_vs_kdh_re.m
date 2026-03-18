function plot_tau_vs_kdh_re(summaryTable, figDir, plotOpts)

if nargin < 3 || ~isfield(plotOpts, 'themes') || isempty(plotOpts.themes)
    plotOpts.themes = "normal";
end

if isempty(summaryTable)
    warning('Summary table is empty. Skipping Tau plot.');
    return;
end

ReVals = unique(summaryTable.Re(:));
if isempty(ReVals)
    warning('No Reynolds groups found. Skipping Tau plot.');
    return;
end

for theme = reshape(plotOpts.themes, 1, [])
    f = figure('Color', 'w', 'Position', [120 120 1000 700]);
    ax = axes(f);
    hold(ax, 'on');

    markerFaceColor = [0 0 1];
    markerEdgeColor = [0 0 0];
    errorBarColor = [1 0 0];
    lgd = gobjects(0,1);
    lgdTxt = strings(0,1);

    for r = 1:numel(ReVals)
        Rei = ReVals(r);
        sub = summaryTable(summaryTable.Re == Rei, :);
        sub = sortrows(sub, 'kDh');

        x = sub.kDh(:);
        y = sub.tau_mean(:);
        e = sub.tau_std(:);

        valid = isfinite(x) & isfinite(y);
        x = x(valid);
        y = y(valid);
        e = e(valid);

        if isempty(x)
            continue;
        end

        e(~isfinite(e)) = 0;

        h = errorbar(ax, x, y, e, 'o', ...
            'LineStyle', 'none', ...
            'LineWidth', 1.8, ...
            'MarkerSize', 7, ...
            'CapSize', 8, ...
            'Color', errorBarColor, ...
            'MarkerFaceColor', markerFaceColor, ...
            'MarkerEdgeColor', markerEdgeColor);
        lgd(end+1,1) = h; %#ok<AGROW>
        lgdTxt(end+1,1) = sprintf('Re=%g', Rei); %#ok<AGROW>
    end

    xlabel(ax, '$k/D_h$', 'Interpreter', 'latex');
    ylabel(ax, '$\tau\;(\mathrm{s})$', 'Interpreter', 'latex');
    title(ax, 'Average upstream residence time vs k/D_h');
    grid(ax, 'on');
    box(ax, 'on');

    if ~isempty(lgd)
        leg = legend(ax, lgd, cellstr(lgdTxt), 'Location', 'southoutside', 'NumColumns', 2, 'Box', 'off');
    else
        leg = [];
    end

    apply_plot_theme(ax, char(theme));
    style_legend_for_theme(leg, char(theme));

    outBase = fullfile(figDir, "Tau_vs_kDh_by_Re_" + theme);
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
