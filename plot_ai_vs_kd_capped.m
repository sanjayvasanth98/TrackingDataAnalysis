function plot_ai_vs_kd_capped(allLoc, figDir, plotOpts)
% PLOT_AI_VS_KD_CAPPED  A/I plot with activation count capped at maxActivationsPerCase.
%
%   For each case: A_capped = min(nActivated, maxActivationsPerCase).
%   Plots A_capped / nInjected vs k/d per Re group.
%   Mirrors the style of plot_ai_vs_kdh_re but without CI bars.

if nargin < 3 || ~isfield(plotOpts, 'themes') || isempty(plotOpts.themes)
    plotOpts.themes = "normal";
end

maxAct = 250;
if isfield(plotOpts, 'maxActivationsPerCase') && isfinite(plotOpts.maxActivationsPerCase)
    maxAct = plotOpts.maxActivationsPerCase;
end

nCases = numel(allLoc.caseName);
if nCases == 0
    warning('No cases in allLoc. Skipping capped A/I plot.');
    return;
end

if ~isfield(allLoc, 'nActivated') || ~isfield(allLoc, 'nInjected')
    warning('allLoc missing nActivated/nInjected fields. Skipping capped A/I plot.');
    return;
end

ReVals = unique(allLoc.Re(:));
if isempty(ReVals)
    warning('No Reynolds numbers found. Skipping capped A/I plot.');
    return;
end

for theme = reshape(plotOpts.themes, 1, [])
    fontName = resolve_plot_font_name();
    f = figure('Color', 'w', 'Position', [100 100 1000 700]);
    ax = axes(f);
    hold(ax, 'on');

    cmap = lines(max(numel(ReVals), 1));
    markerStyles = {'o', 's', '^', 'd', 'v', '>', '<', 'p', 'h'};
    lgd    = gobjects(0,1);
    lgdTxt = strings(0,1);

    for r = 1:numel(ReVals)
        Rei    = ReVals(r);
        idxRe  = find(allLoc.Re == Rei);
        col = cmap(r, :);
        markerStyle = markerStyles{mod(r - 1, numel(markerStyles)) + 1};

        kD_vals  = allLoc.kD(idxRe);
        nAct     = allLoc.nActivated(idxRe);
        nInj     = allLoc.nInjected(idxRe);

        A_capped = min(nAct, maxAct);
        AI_capped = A_capped ./ max(nInj, 1);

        valid = isfinite(kD_vals) & isfinite(AI_capped) & AI_capped > 0;
        kD_v  = kD_vals(valid);
        AI_v  = AI_capped(valid);

        if isempty(kD_v), continue; end

        [kD_sorted, si] = sort(kD_v);
        AI_sorted = AI_v(si);

        hPts = plot(ax, kD_sorted, AI_sorted, [markerStyle '-'], ...
            'Color', col, ...
            'LineWidth', 0.75, ...
            'LineStyle', '-', ...
            'MarkerSize', 12, ...
            'MarkerFaceColor', col, ...
            'MarkerEdgeColor', [0 0 0]);
        lgd(end+1,1)    = hPts; %#ok<AGROW>
        lgdTxt(end+1,1) = sprintf('Data (cap=%d), Re=%g', maxAct, Rei); %#ok<AGROW>

        if numel(kD_sorted) >= 2
            p = polyfit(kD_sorted, log10(AI_sorted), 1);
            xFit = linspace(min(kD_sorted), max(kD_sorted), 200).';
            yFit = 10.^(p(2) + p(1)*xFit);
            hFit = plot(ax, xFit, yFit, '--', 'LineWidth', 1.8, 'Color', col);
            lgd(end+1,1)    = hFit; %#ok<AGROW>
            lgdTxt(end+1,1) = sprintf('Fit (cap=%d), Re=%g', maxAct, Rei); %#ok<AGROW>
        end
    end

    xlabel(ax, '$k/d$', 'Interpreter', 'latex');
    ylabel(ax, sprintf('A_{cap}/I (cap = %d)', maxAct), 'Interpreter', 'tex');
    title(ax, '');
    set(ax, 'YScale', 'linear', 'FontName', fontName);
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

    outBase = fullfile(figDir, sprintf('AI_capped%d_vs_kD_%s', maxAct, theme));
    save_fig_dual_safe(f, outBase, plotOpts);
    if ~isfield(plotOpts, 'keepFiguresOpen') || ~plotOpts.keepFiguresOpen
        close(f);
    end
end

fprintf('Saved capped A/I plot (cap=%d) to: %s\n', maxAct, figDir);
end


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
