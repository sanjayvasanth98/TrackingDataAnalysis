function plot_void_fraction_vs_kd(allVoidFrac, figDir, plotOpts)
%PLOT_VOID_FRACTION_VS_KD  Mean 2-D void fraction (%) vs k/d for all cases.
%
%   One marker per case.  Y-axis: mean void fraction in percent.
%   An exponential fit (log-linear in k/d) is overlaid when >= 2 valid cases.
%   Mirrors the style of plot_ai_vs_kd_capped.

if nargin < 3 || ~isfield(plotOpts,'themes') || isempty(plotOpts.themes)
    plotOpts.themes = "normal";
end

nCases = numel(allVoidFrac.caseName);
if nCases == 0
    warning('plot_void_fraction_vs_kd: no cases. Skipping.');
    return;
end

% Collect valid mean void fractions
kD_all  = allVoidFrac.kD(:);
Re_all  = allVoidFrac.Re(:);
vf_all  = nan(nCases, 1);
for ci = 1:nCases
    d = allVoidFrac.data{ci};
    if ~isempty(d) && isfield(d,'vfMean') && isfinite(d.vfMean)
        vf_all(ci) = d.vfMean * 100;  % convert to percent
    end
end

valid = isfinite(kD_all) & isfinite(vf_all) & vf_all > 0;
if ~any(valid)
    warning('plot_void_fraction_vs_kd: no cases with finite positive void fraction. Skipping.');
    return;
end

for theme = reshape(plotOpts.themes, 1, [])
    fontName = resolve_plot_font_name();
    f  = figure('Color','w','Position',[100 100 1000 700]);
    ax = axes(f);
    hold(ax,'on');

    lgd    = gobjects(0,1);
    lgdTxt = strings(0,1);

    ReVals = unique(Re_all(valid));
    cmap = lines(max(numel(ReVals), 1));
    markerStyles = {'o', 's', '^', 'd', 'v', '>', '<', 'p', 'h'};

    for r = 1:numel(ReVals)
        Rei   = ReVals(r);
        mask  = valid & Re_all == Rei;
        kD_v  = kD_all(mask);
        vf_v  = vf_all(mask);
        col = cmap(r, :);
        markerStyle = markerStyles{mod(r - 1, numel(markerStyles)) + 1};

        [kD_sorted, si] = sort(kD_v);
        vf_sorted = vf_v(si);

        hPts = plot(ax, kD_sorted, vf_sorted, [markerStyle '-'], ...
            'Color',           col, ...
            'LineWidth',       0.75, ...
            'MarkerSize',      12, ...
            'MarkerFaceColor', col, ...
            'MarkerEdgeColor', [0 0 0]);
        lgd(end+1,1)    = hPts; %#ok<AGROW>
        lgdTxt(end+1,1) = sprintf('Data, Re=%g', Rei); %#ok<AGROW>

        if numel(kD_sorted) >= 2
            p = polyfit(kD_sorted, log10(vf_sorted), 1);
            xFit = linspace(min(kD_sorted), max(kD_sorted), 200).';
            yFit = 10.^(p(2) + p(1)*xFit);
            hFit = plot(ax, xFit, yFit, '--', 'LineWidth', 1.8, 'Color', col);
            lgd(end+1,1)    = hFit; %#ok<AGROW>
            lgdTxt(end+1,1) = sprintf('Fit, Re=%g', Rei); %#ok<AGROW>
        end
    end

    xlabel(ax, '$k/d$',                        'Interpreter','latex');
    ylabel(ax, 'Mean void fraction $\bar{\alpha}$ (\%)', 'Interpreter','latex');
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

    outBase = fullfile(figDir, "VoidFraction_vs_kD_" + theme);
    save_fig_dual_safe(f, outBase, plotOpts);
    if ~isfield(plotOpts,'keepFiguresOpen') || ~plotOpts.keepFiguresOpen
        close(f);
    end
end
fprintf('Saved void fraction plot to: %s\n', figDir);
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
