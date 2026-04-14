function plot_tau_velocity_colormap(summaryTable, figDir, plotOpts, matDir)
%PLOT_TAU_VELOCITY_COLORMAP  Residence time coloured by activated velocity.
%
%   Creates tau-vs-k/d plots where each case marker is coloured by the
%   average upstream axial velocity of strict activated upstream-moving
%   bubbles. Velocity is in m/s and the colour scale is shared across all
%   cases in the input table.

if nargin < 3 || ~isfield(plotOpts, 'themes') || isempty(plotOpts.themes)
    plotOpts.themes = "normal";
end
if nargin < 4
    matDir = "";
end
if ~isfolder(figDir), mkdir(figDir); end

required = {'Re','kD','tau_mean','tau_sem','activatedVelocity_mean_m_s'};
missing = required(~ismember(required, summaryTable.Properties.VariableNames));
if ~isempty(missing)
    warning('plot_tau_velocity_colormap: missing columns: %s. Skipping.', strjoin(missing, ', '));
    return;
end

plotData = build_plot_data(summaryTable);
if isempty(plotData.kD)
    warning('plot_tau_velocity_colormap: no finite tau/velocity data. Skipping.');
    return;
end

if strlength(string(matDir)) > 0
    if ~isfolder(matDir), mkdir(matDir); end
    save(fullfile(matDir, "tau_velocity_colormap_plot_data.mat"), 'plotData');
end

% Normalization parameters match plot_tau_vs_kdh_re.
U_throat = 13.32;       % m/s
H_throat_m = 10e-3;     % m
t_conv = H_throat_m / U_throat;

for theme = reshape(plotOpts.themes, 1, [])
    plot_one_tau_velocity(plotData, figDir, plotOpts, char(theme), ...
        1, '$\tau\;(\mathrm{s})$', "Tau_vs_kD_velocity_colormap_" + theme);

    plot_one_tau_velocity(plotData, figDir, plotOpts, char(theme), ...
        1/t_conv, '$\tau^{*} = \tau \, U / H$', "Tau_normalized_vs_kD_velocity_colormap_" + theme);
end

fprintf('Saved velocity-coloured residence time plots to: %s\n', figDir);
end


function plotData = build_plot_data(summaryTable)
v = summaryTable.activatedVelocity_mean_m_s(:);
valid = isfinite(summaryTable.Re(:)) & isfinite(summaryTable.kD(:)) & ...
    isfinite(summaryTable.tau_mean(:)) & isfinite(v);

plotData = struct();
plotData.caseName = strings(0,1);
if ismember('Case', summaryTable.Properties.VariableNames)
    plotData.caseName = string(summaryTable.Case(valid));
end
plotData.Re = summaryTable.Re(valid);
plotData.kD = summaryTable.kD(valid);
plotData.tau_mean = summaryTable.tau_mean(valid);
plotData.tau_sem = summaryTable.tau_sem(valid);
plotData.activatedVelocity_mean_m_s = summaryTable.activatedVelocity_mean_m_s(valid);

if ismember('activatedVelocity_std_m_s', summaryTable.Properties.VariableNames)
    plotData.activatedVelocity_std_m_s = summaryTable.activatedVelocity_std_m_s(valid);
else
    plotData.activatedVelocity_std_m_s = nan(size(plotData.kD));
end
if ismember('activatedVelocity_sem_m_s', summaryTable.Properties.VariableNames)
    plotData.activatedVelocity_sem_m_s = summaryTable.activatedVelocity_sem_m_s(valid);
else
    plotData.activatedVelocity_sem_m_s = nan(size(plotData.kD));
end
if ismember('activatedVelocity_n', summaryTable.Properties.VariableNames)
    plotData.activatedVelocity_n = summaryTable.activatedVelocity_n(valid);
else
    plotData.activatedVelocity_n = nan(size(plotData.kD));
end
end


function plot_one_tau_velocity(plotData, figDir, plotOpts, theme, yScale, yLabelStr, outName)
fontName = resolve_plot_font_name();
f = figure('Color', 'w', 'Position', [120 120 1050 720]);
ax = axes(f);
hold(ax, 'on');

ReVals = unique(plotData.Re(:));
velocityVals = plotData.activatedVelocity_mean_m_s(:);
cLim = velocity_color_limits(velocityVals);

markerStyles = {'o', 's', '^', 'd', 'v', '>', '<', 'p', 'h'};
lgd = gobjects(0,1);
lgdTxt = strings(0,1);

for r = 1:numel(ReVals)
    Rei = ReVals(r);
    mask = plotData.Re == Rei;
    if ~any(mask)
        continue;
    end

    x = plotData.kD(mask);
    y = plotData.tau_mean(mask) * yScale;
    e = plotData.tau_sem(mask) * yScale;
    c = plotData.activatedVelocity_mean_m_s(mask);

    valid = isfinite(x) & isfinite(y) & isfinite(c);
    x = x(valid);
    y = y(valid);
    e = e(valid);
    c = c(valid);
    if isempty(x)
        continue;
    end

    [x, si] = sort(x);
    y = y(si);
    e = e(si);
    c = c(si);

    e(~isfinite(e)) = 0;
    eLow = min(e, y);
    markerStyle = markerStyles{mod(r - 1, numel(markerStyles)) + 1};

    plot(ax, x, y, '-', ...
        'Color', [0.55 0.55 0.55], ...
        'LineWidth', 0.8, ...
        'HandleVisibility', 'off');
    errorbar(ax, x, y, eLow, e, ...
        'LineStyle', 'none', ...
        'Color', [0.25 0.25 0.25], ...
        'LineWidth', 0.75, ...
        'CapSize', 8, ...
        'HandleVisibility', 'off');
    hSc = scatter(ax, x, y, 95, c, markerStyle, ...
        'filled', ...
        'MarkerEdgeColor', [0 0 0], ...
        'LineWidth', 1.1);
    lgd(end+1,1) = hSc; %#ok<AGROW>
    lgdTxt(end+1,1) = sprintf('Re=%g', Rei); %#ok<AGROW>
end

colormap(ax, turbo_colormap_compat(256));
caxis(ax, cLim);
cb = colorbar(ax);
cb.Label.String = 'Avg activated upstream velocity (m/s)';
cb.Label.Interpreter = 'latex';

xlabel(ax, '$k/d$', 'Interpreter', 'latex');
ylabel(ax, yLabelStr, 'Interpreter', 'latex');
title(ax, '');
grid(ax, 'off');
box(ax, 'on');
set(ax, 'FontName', fontName);

if ~isempty(lgd)
    leg = legend(ax, lgd, cellstr(lgdTxt), ...
        'Location', 'southoutside', 'NumColumns', 2, 'Box', 'off');
else
    leg = [];
end

apply_plot_theme(ax, theme);
style_legend_and_colorbar_for_theme(leg, cb, theme);

outBase = fullfile(figDir, outName);
save_fig_dual_safe(f, outBase, plotOpts);
if ~isfield(plotOpts, 'keepFiguresOpen') || ~plotOpts.keepFiguresOpen
    close(f);
end
end


function cLim = velocity_color_limits(v)
v = v(isfinite(v));
if isempty(v)
    cLim = [0 1];
    return;
end
vMin = min(v);
vMax = max(v);
if ~(isfinite(vMin) && isfinite(vMax)) || vMax <= vMin
    pad = max(abs(vMin) * 0.05, 0.01);
    cLim = [vMin - pad, vMin + pad];
else
    pad = 0.04 * (vMax - vMin);
    cLim = [vMin - pad, vMax + pad];
end
end


function cmap = turbo_colormap_compat(n)
if nargin < 1 || isempty(n)
    n = 256;
end
try
    cmap = turbo(n);
catch
    % Older MATLAB fallback. The plot still works, but newer MATLAB uses turbo.
    cmap = jet(n);
end
end


function style_legend_and_colorbar_for_theme(leg, cb, theme)
if ~isempty(leg) && isgraphics(leg)
    if strcmp(theme, 'poster')
        leg.TextColor = [1 1 1];
        leg.Color = [0 0 0];
    else
        leg.TextColor = [0 0 0];
        leg.Color = 'none';
    end
end

if ~isempty(cb) && isgraphics(cb)
    if strcmp(theme, 'poster')
        cb.Color = [1 1 1];
    else
        cb.Color = [0 0 0];
    end
end
end
