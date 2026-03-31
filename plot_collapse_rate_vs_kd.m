function plot_collapse_rate_vs_kd(allCollapse, figDir, plotOpts)
%PLOT_COLLAPSE_RATE_VS_KD  Collapse rate (events/s) vs k/d for all cases.
%
%   One marker per case.  Y-axis: collapse rate in events per second.
%   An exponential fit (log-linear in k/d) is overlaid when >= 2 valid cases.
%   Mirrors the style of plot_void_fraction_vs_kd.

if nargin < 3 || ~isfield(plotOpts,'themes') || isempty(plotOpts.themes)
    plotOpts.themes = "normal";
end

nCases = numel(allCollapse.caseName);
if nCases == 0
    warning('plot_collapse_rate_vs_kd: no cases. Skipping.');
    return;
end

% Collect valid collapse rates
kD_all   = allCollapse.kD(:);
Re_all   = allCollapse.Re(:);
rate_all = nan(nCases, 1);
for ci = 1:nCases
    d = allCollapse.data{ci};
    if ~isempty(d) && isfield(d,'ratePerSec') && isfinite(d.ratePerSec)
        rate_all(ci) = d.ratePerSec;
    end
end

valid = isfinite(kD_all) & isfinite(rate_all) & rate_all > 0;
if ~any(valid)
    warning('plot_collapse_rate_vs_kd: no cases with finite positive collapse rate. Skipping.');
    return;
end

for theme = reshape(plotOpts.themes, 1, [])
    fontName = resolve_plot_font_name();
    f  = figure('Color','w','Position',[100 100 1000 700]);
    ax = axes(f);
    hold(ax,'on');

    markerFaceColor = [0 0 1];
    markerEdgeColor = [0 0 0];
    fitColor        = [0.35 0.35 0.35];
    lgd    = gobjects(0,1);
    lgdTxt = strings(0,1);

    ReVals = unique(Re_all(valid));

    for r = 1:numel(ReVals)
        Rei   = ReVals(r);
        mask  = valid & Re_all == Rei;
        kD_v  = kD_all(mask);
        rt_v  = rate_all(mask);

        [kD_sorted, si] = sort(kD_v);
        rt_sorted = rt_v(si);

        hPts = plot(ax, kD_sorted, rt_sorted, 'o-', ...
            'Color',           markerEdgeColor, ...
            'LineWidth',       0.75, ...
            'MarkerSize',      12, ...
            'MarkerFaceColor', markerFaceColor, ...
            'MarkerEdgeColor', markerEdgeColor);
        lgd(end+1,1)    = hPts; %#ok<AGROW>
        lgdTxt(end+1,1) = sprintf('Data, Re=%g', Rei); %#ok<AGROW>

        if numel(kD_sorted) >= 2
            p = polyfit(kD_sorted, log10(rt_sorted), 1);
            xFit = linspace(min(kD_sorted), max(kD_sorted), 200).';
            yFit = 10.^(p(2) + p(1)*xFit);
            hFit = plot(ax, xFit, yFit, '--', 'LineWidth', 1.8, 'Color', fitColor);
            lgd(end+1,1)    = hFit; %#ok<AGROW>
            lgdTxt(end+1,1) = sprintf('Fit, Re=%g', Rei); %#ok<AGROW>
        end
    end

    xlabel(ax, '$k/d$', 'Interpreter','latex');
    ylabel(ax, 'Collapse rate $\dot{N}_c\;(\mathrm{s}^{-1})$', 'Interpreter','latex');
    title(ax,  '');
    set(ax, 'YScale','linear', 'FontName', fontName);
    grid(ax, 'off');
    box(ax,  'on');

    if ~isempty(lgd)
        leg = legend(ax, lgd, cellstr(lgdTxt), ...
            'Location','southoutside','NumColumns',2,'Box','off');
    else
        leg = [];
    end

    apply_plot_theme(ax, char(theme));
    style_legend_for_theme(leg, char(theme));

    outBase = fullfile(figDir, "CollapseRate_vs_kD_" + theme);
    save_fig_dual_safe(f, outBase, plotOpts);
    if ~isfield(plotOpts,'keepFiguresOpen') || ~plotOpts.keepFiguresOpen
        close(f);
    end
end
fprintf('Saved collapse rate vs k/d plot to: %s\n', figDir);
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
