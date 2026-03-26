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

% Normalization parameters
U_throat = 13.32;       % throat velocity (m/s)
H_throat_m = 10e-3;     % throat height (m)
t_conv = H_throat_m / U_throat; % convective time scale (s)

for theme = reshape(plotOpts.themes, 1, [])

    % --- Plot 1: dimensional tau (s) ---
    plot_tau_panel(summaryTable, ReVals, figDir, theme, plotOpts, ...
        1, '$\tau\;(\mathrm{s})$', 'Average upstream residence time vs k/d', ...
        "Tau_vs_kD_by_Re_" + theme);

    % --- Plot 2: normalized tau* = tau / t_conv ---
    plot_tau_panel(summaryTable, ReVals, figDir, theme, plotOpts, ...
        1/t_conv, '$\tau^{*} = \tau \, U / H$', 'Normalized upstream residence time vs k/d', ...
        "Tau_normalized_vs_kD_by_Re_" + theme);

end

end

%% ---- shared plotting helper ----
function plot_tau_panel(summaryTable, ReVals, figDir, theme, plotOpts, ...
    yScale, yLabelStr, titleStr, outName)

f = figure('Color', 'w', 'Position', [120 120 1000 700]);
ax = axes(f);
hold(ax, 'on');

markerFaceColor = [0 0 1];
markerEdgeColor = [0 0 0];
errorBarColor = [0 0 0];
lgd = gobjects(0,1);
lgdTxt = strings(0,1);

for r = 1:numel(ReVals)
    Rei = ReVals(r);
    sub = summaryTable(summaryTable.Re == Rei, :);
    sub = sortrows(sub, 'kD');

    x = sub.kD(:);
    y = sub.tau_mean(:) * yScale;
    e = sub.tau_std(:) * yScale;

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
        'LineWidth', 0.75, ...
        'MarkerSize', 12, ...
        'CapSize', 8, ...
        'Color', errorBarColor, ...
        'MarkerFaceColor', markerFaceColor, ...
        'MarkerEdgeColor', markerEdgeColor);
    lgd(end+1,1) = h; %#ok<AGROW>
    lgdTxt(end+1,1) = sprintf('Re=%g', Rei); %#ok<AGROW>
end

xlabel(ax, '$k/d$', 'Interpreter', 'latex');
ylabel(ax, yLabelStr, 'Interpreter', 'latex');
title(ax, titleStr, 'FontName', 'Times New Roman', 'FontSize', 12);
grid(ax, 'off');
box(ax, 'on');

if ~isempty(lgd)
    leg = legend(ax, lgd, cellstr(lgdTxt), 'Location', 'southoutside', 'NumColumns', 2, 'Box', 'off');
else
    leg = [];
end

apply_plot_theme(ax, char(theme));
style_legend_for_theme(leg, char(theme));

outBase = fullfile(figDir, outName);
save_fig_dual_safe(f, outBase, plotOpts);
close(f);
end

%% ---- theme helpers ----
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
