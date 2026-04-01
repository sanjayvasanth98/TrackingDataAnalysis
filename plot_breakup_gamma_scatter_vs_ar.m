function plot_breakup_gamma_scatter_vs_ar(breakupData, figDir, plotOpts, arLabel, matDir)
% PLOT_BREAKUP_GAMMA_SCATTER_VS_AR
%   Scatter of gamma = (x_child - x_parent)/d vs parent aspect ratio.
%   One colour per k/d case, 50% transparent markers, no edges.
%
%   breakupData : struct array with .caseName, .kD, .events
%   figDir      : output directory for figures
%   plotOpts    : standard plotting options struct
%   arLabel     : (optional) AR threshold tag, e.g. 'AR1p5'
%   matDir      : (optional) directory to save .mat plot data

if nargin < 3 || ~isfield(plotOpts, 'themes') || isempty(plotOpts.themes)
    plotOpts.themes = "normal";
end
if nargin < 4 || isempty(arLabel), arLabel = ""; end
if nargin < 5, matDir = ""; end

nCases = numel(breakupData);
if nCases == 0, warning('No breakup data.'); return; end

% ---- Collect all events with case metadata ----
allAR    = [];
allGamma = [];
allKD    = [];
caseIdx  = [];
for ci = 1:nCases
    ev = breakupData(ci).events;
    if isempty(ev), continue; end
    nEv = numel(ev);
    allAR    = [allAR;    [ev.parentAR].'];  %#ok<AGROW>
    allGamma = [allGamma; [ev.gamma].'];     %#ok<AGROW>
    allKD    = [allKD;    repmat(breakupData(ci).kD, nEv, 1)]; %#ok<AGROW>
    caseIdx  = [caseIdx;  repmat(ci, nEv, 1)];                %#ok<AGROW>
end

if isempty(allAR)
    warning('No breakup events found across all cases.');
    return;
end

% ---- Unique cases for colouring ----
casesPresent = unique(caseIdx, 'stable');
nPresent     = numel(casesPresent);
cmap = lines(max(nPresent, 1));

% ---- Save .mat ----
if matDir ~= ""
    plotData = struct();
    plotData.parentAR  = allAR;
    plotData.gamma     = allGamma;
    plotData.kD        = allKD;
    plotData.caseIdx   = caseIdx;
    plotData.arLabel   = arLabel;
    caseNames = strings(nCases, 1);
    caseKDs   = nan(nCases, 1);
    for ci = 1:nCases
        caseNames(ci) = string(breakupData(ci).caseName);
        caseKDs(ci)   = breakupData(ci).kD;
    end
    plotData.caseNames = caseNames;
    plotData.caseKDs   = caseKDs;
    if arLabel ~= ""
        matFile = fullfile(matDir, sprintf("breakup_gamma_vs_ar_%s.mat", arLabel));
    else
        matFile = fullfile(matDir, "breakup_gamma_vs_ar.mat");
    end
    save(matFile, 'plotData');
end

% ---- Axis ranges ----
xMin = floor(min(allAR));
xMax = ceil(max(allAR));
if xMax <= xMin, xMax = xMin + 1; end

yAbsMax = ceil(max(abs(allGamma)));
if yAbsMax == 0, yAbsMax = 1; end

% ---- Plot per theme ----
for theme = reshape(plotOpts.themes, 1, [])
    fontName = resolve_plot_font_name();
    f  = figure('Color', 'w', 'Position', [100 100 1000 700]);
    ax = axes(f); %#ok<LAXES>
    hold(ax, 'on');

    lgd    = gobjects(0, 1);
    lgdTxt = strings(0, 1);

    for p = 1:nPresent
        ci  = casesPresent(p);
        sel = caseIdx == ci;
        col = cmap(p, :);
        nSel = sum(sel);
        if nSel == 0, continue; end

        hPt = scatter(ax, allAR(sel), allGamma(sel), 50, col, 'filled', ...
            'MarkerFaceAlpha', 0.50, ...
            'MarkerEdgeColor', 'none');
        lgd(end+1, 1)    = hPt; %#ok<AGROW>
        lgdTxt(end+1, 1) = sprintf('k/d = %.2f  (n = %d)', ...
            breakupData(ci).kD, nSel); %#ok<AGROW>
    end

    % ---- Zero reference ----
    plot(ax, [xMin, xMax], [0 0], ':', ...
        'Color', [0.55 0.55 0.55], 'LineWidth', 1.0, ...
        'HandleVisibility', 'off');

    % ---- Axes ----
    xlim(ax, [xMin, xMax]);
    ylim(ax, [-yAbsMax, yAbsMax]);
    set(ax, 'FontName', fontName);

    xlabel(ax, 'Parent aspect ratio, $\mathrm{AR}_\mathrm{parent}$', ...
        'Interpreter', 'latex');
    ylabel(ax, '$\gamma = (x_\mathrm{child} - x_\mathrm{parent})\,/\,d$', ...
        'Interpreter', 'latex');
    if arLabel ~= ""
        title(ax, sprintf('Parent AR $\\geq$ %s', ...
            strrep(char(arLabel), 'AR', '')), 'Interpreter', 'latex');
    else
        title(ax, '');
    end
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

    if arLabel ~= ""
        outBase = fullfile(figDir, "Breakup_gamma_scatter_AR_" + arLabel + "_" + theme);
    else
        outBase = fullfile(figDir, "Breakup_gamma_scatter_AR_" + theme);
    end
    save_fig_dual_safe(f, outBase, plotOpts);
    if ~isfield(plotOpts, 'keepFiguresOpen') || ~plotOpts.keepFiguresOpen
        close(f);
    end
end

fprintf('Saved breakup scatter (gamma vs AR) to: %s\n', figDir);
end


% =========================================================================
function style_legend_for_theme(leg, theme)
if isempty(leg) || ~isgraphics(leg), return; end
if strcmp(theme, 'poster')
    leg.TextColor = [1 1 1];
    leg.Color     = [0 0 0];
else
    leg.TextColor = [0 0 0];
    leg.Color     = 'none';
end
end
