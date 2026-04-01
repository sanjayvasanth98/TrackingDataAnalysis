function plot_breakup_gamma_vs_dratio(breakupData, figDir, plotOpts, arLabel)
% PLOT_BREAKUP_GAMMA_VS_DRATIO
%   Scatter of gamma = (x_child - x_parent)/d_roughness vs d_child/d_parent.
%   One scatter point per child-parent pair, 50% transparent markers.
%   Per-case coloured mean line (linear bins, no markers).
%   Zero-line reference at gamma=0.
%   Linear x-axis.  Y centred at 0 with symmetric limits.
%
%   breakupData: struct array, one element per case, with fields:
%     .caseName  string
%     .kD        double
%     .events    struct array from analyze_breakup_events
%
%   arLabel (optional): string appended to output filename, e.g. 'AR1p5'.
%     Also shown in the title for identification.

if nargin < 3 || ~isfield(plotOpts, 'themes') || isempty(plotOpts.themes)
    plotOpts.themes = "normal";
end
if nargin < 4 || isempty(arLabel)
    arLabel = "";
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

% X-axis: linear scale with tight limits around the data.
xPad = 0.05 * (max(allDRatio) - min(allDRatio));
if xPad == 0, xPad = 0.5; end
xLim = [max(0, min(allDRatio) - xPad), max(allDRatio) + xPad];

% Y-axis: symmetric about 0.
yAbsMax = ceil(max(abs(allGamma)));
if yAbsMax == 0, yAbsMax = 1; end
yLimPlot = [-yAbsMax, yAbsMax];

% Per-case colour (using lines colourmap).
cmap = lines(max(nCases, 1));

% Number of bins for the mean line (equally spaced in log space).
nBins = 8;

for theme = reshape(plotOpts.themes, 1, [])
    f = figure('Color', 'w', 'Position', [100 100 1100 700]);
    ax = axes(f);
    hold(ax, 'on');

    lgd    = gobjects(0,1);
    lgdTxt = strings(0,1);

    for ci = 1:nCases
        ev = breakupData(ci).events;
        if isempty(ev), continue; end

        dRatio   = [ev.dRatio].';
        gamma    = [ev.gamma].';
        col      = cmap(ci, :);

        % --- Scatter (50% transparent) — used for legend ---
        hPt = scatter(ax, dRatio, gamma, 40, col, 'filled', ...
            'MarkerFaceAlpha', 0.50, ...
            'MarkerEdgeColor', 'none');
        lgd(end+1,1)    = hPt; %#ok<AGROW>
        lgdTxt(end+1,1) = sprintf('k/d = %.2f', breakupData(ci).kD); %#ok<AGROW>

        % --- Binned mean line (equally spaced in linear space, no markers) ---
        binEdges  = linspace(xLim(1), xLim(2), nBins + 1);
        binCentre = 0.5 * (binEdges(1:end-1) + binEdges(2:end));
        meanGamma = nan(nBins, 1);
        for b = 1:nBins
            inBin = dRatio >= binEdges(b) & dRatio < binEdges(b+1);
            if sum(inBin) >= 1
                meanGamma(b) = mean(gamma(inBin));
            end
        end
        hasData = isfinite(meanGamma);
        if sum(hasData) >= 2
            plot(ax, binCentre(hasData), meanGamma(hasData), '-', ...
                'Color', col, ...
                'LineWidth', 2.0, ...
                'HandleVisibility', 'off');
        end
    end

    % Zero reference line.
    plot(ax, xLim, [0 0], '--', 'Color', [0.5 0.5 0.5], ...
        'LineWidth', 1.2, 'HandleVisibility', 'off');

    xlim(ax, xLim);
    ylim(ax, yLimPlot);

    xlabel(ax, '$d_\mathrm{child}/d_\mathrm{parent}$', 'Interpreter', 'latex');
    ylabel(ax, '$\gamma = (x_\mathrm{child} - x_\mathrm{parent})\,/\,d$', ...
        'Interpreter', 'latex');
    grid(ax, 'off');
    box(ax, 'on');

    if ~isempty(lgd)
        leg = legend(ax, lgd, cellstr(lgdTxt), ...
            'Location', 'northeast', 'Box', 'off');
    else
        leg = [];
    end

    apply_plot_theme(ax, char(theme));
    style_legend_for_theme(leg, char(theme));

    if arLabel ~= ""
        outBase = fullfile(figDir, "Breakup_gamma_vs_dRatio_" + arLabel + "_" + theme);
    else
        outBase = fullfile(figDir, "Breakup_gamma_vs_dRatio_" + theme);
    end
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
