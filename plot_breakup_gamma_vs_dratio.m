function plot_breakup_gamma_vs_dratio(breakupData, figDir, plotOpts)
% PLOT_BREAKUP_GAMMA_VS_DRATIO
%   Scatter of gamma = (x_child - x_parent)/d_roughness vs d_child/d_parent.
%   One scatter point per child-parent pair, 50% transparent markers.
%   Per-case coloured mean line (log-spaced bins, no markers).
%   Zero-line reference at gamma=0.
%   Square axes.  X plotted as log10(dRatio) with 10^n tick labels.
%   Y centred at 0 with symmetric limits.
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

% Pool all data to determine sensible axis ranges.
allDRatio = [];
allGamma  = [];
for ci = 1:nCases
    ev = breakupData(ci).events;
    if isempty(ev), continue; end
    allDRatio = [allDRatio; [ev.dRatio].']; %#ok<AGROW>
    allGamma  = [allGamma;  [ev.gamma].'];  %#ok<AGROW>
end
if isempty(allDRatio)
    warning('No breakup events found across all cases.');
    return;
end

% X-axis in log10 space (linear axis, labels as 10^n).
logDR = log10(max(allDRatio, 1e-6));       % guard against log(0)
xMinExp = floor(min(logDR));               % e.g. -2
xMaxExp = ceil(max(logDR));                % e.g.  0
if xMaxExp <= xMinExp, xMaxExp = xMinExp + 1; end
xLimLog = [xMinExp, xMaxExp];

% Tick positions and labels: 10^-2, 10^-1, 10^0 ...
xTicks   = xMinExp:xMaxExp;
xTickLbl = arrayfun(@(e) sprintf('10^{%d}', e), xTicks, 'UniformOutput', false);

% Y-axis: symmetric about 0.
yAbsMax = ceil(max(abs(allGamma)));
if yAbsMax == 0, yAbsMax = 1; end
yLimPlot = [-yAbsMax, yAbsMax];

% Per-case colour (using lines colourmap).
cmap = lines(max(nCases, 1));

% Number of bins for the mean line (equally spaced in log10 space).
nBins = 8;

for theme = reshape(plotOpts.themes, 1, [])
    f = figure('Color', 'w', 'Position', [100 100 700 700]);
    ax = axes(f);
    hold(ax, 'on');

    lgd    = gobjects(0,1);
    lgdTxt = strings(0,1);

    for ci = 1:nCases
        ev = breakupData(ci).events;
        if isempty(ev), continue; end

        dRatio   = [ev.dRatio].';
        gamma    = [ev.gamma].';
        logDRci  = log10(max(dRatio, 1e-6));
        col      = cmap(ci, :);

        % --- Scatter (50% transparent) ---
        scatter(ax, logDRci, gamma, 40, col, 'filled', ...
            'MarkerFaceAlpha', 0.50, ...
            'MarkerEdgeColor', 'none', ...
            'HandleVisibility', 'off');

        % --- Binned mean line (equally spaced in log10 space, no markers) ---
        binEdges  = linspace(xLimLog(1), xLimLog(2), nBins + 1);
        binCentre = 0.5 * (binEdges(1:end-1) + binEdges(2:end));
        meanGamma = nan(nBins, 1);
        for b = 1:nBins
            inBin = logDRci >= binEdges(b) & logDRci < binEdges(b+1);
            if sum(inBin) >= 1
                meanGamma(b) = mean(gamma(inBin));
            end
        end
        hasData = isfinite(meanGamma);
        if sum(hasData) >= 2
            hLine = plot(ax, binCentre(hasData), meanGamma(hasData), '-', ...
                'Color', col, ...
                'LineWidth', 2.0);
            lgd(end+1,1)    = hLine; %#ok<AGROW>
            lgdTxt(end+1,1) = sprintf('k/d = %.2f', breakupData(ci).kD); %#ok<AGROW>
        end
    end

    % Zero reference line.
    plot(ax, xLimLog, [0 0], '--', 'Color', [0.5 0.5 0.5], ...
        'LineWidth', 1.2, 'HandleVisibility', 'off');

    xlim(ax, xLimLog);
    ylim(ax, yLimPlot);
    set(ax, 'XTick', xTicks, 'XTickLabel', xTickLbl, 'TickLabelInterpreter', 'tex');
    pbaspect(ax, [1 1 1]);

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
