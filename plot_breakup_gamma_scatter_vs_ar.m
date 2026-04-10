function plot_breakup_gamma_scatter_vs_ar(breakupData, figDir, plotOpts, arLabel, matDir)
% PLOT_BREAKUP_GAMMA_SCATTER_VS_AR
%   Scatter of gamma = (x_child - x_parent)/d vs parent aspect ratio.
%   One colour per k/d case, lightly transparent markers, no edges.
%   A single density-based contour is added per case to show the dominant
%   cluster location and spread.
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

markerAlpha = get_opt_value(plotOpts, 'breakupARMarkerAlpha', 0.22);
densityMass = get_opt_value(plotOpts, 'breakupARDensityMass', 0.60);
densityGridSize = round(get_opt_value(plotOpts, 'breakupARDensityGridSize', 80));
densityGridSize = max(densityGridSize, 40);
clusterMinCount = round(get_opt_value(plotOpts, 'breakupARClusterMinCount', 5));
clusterMinCount = max(clusterMinCount, 3);

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

validAll = isfinite(allAR) & isfinite(allGamma) & isfinite(caseIdx);
allAR = allAR(validAll);
allGamma = allGamma(validAll);
allKD = allKD(validAll);
caseIdx = caseIdx(validAll);
if isempty(allAR)
    warning('No finite breakup events found across all cases.');
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
xPad = 0.05 * (max(allAR) - min(allAR));
if ~isfinite(xPad) || xPad <= 0, xPad = 0.5; end
xLim = [max(0, min(allAR) - xPad), max(allAR) + xPad];
if xLim(2) <= xLim(1), xLim(2) = xLim(1) + 1; end

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

        draw_density_contour(ax, allAR(sel), allGamma(sel), col, densityMass, densityGridSize, clusterMinCount);

        hPt = scatter(ax, allAR(sel), allGamma(sel), 50, col, 'filled', ...
            'MarkerFaceAlpha', markerAlpha, ...
            'MarkerEdgeColor', 'none');
        lgd(end+1, 1)    = hPt; %#ok<AGROW>
        lgdTxt(end+1, 1) = sprintf('k/d = %.2f  (n = %d)', ...
            breakupData(ci).kD, nSel); %#ok<AGROW>
    end

    % ---- Zero reference ----
    plot(ax, xLim, [0 0], ':', ...
        'Color', [0.55 0.55 0.55], 'LineWidth', 1.0, ...
        'HandleVisibility', 'off');

    % ---- Axes ----
    xlim(ax, xLim);
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
            'Location', 'best', 'NumColumns', 1, 'Box', 'on');
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
function draw_density_contour(ax, xVals, yVals, col, targetMass, gridSize, minCount)
valid = isfinite(xVals) & isfinite(yVals);
xVals = xVals(valid);
yVals = yVals(valid);
if numel(xVals) < minCount
    return;
end

xVals = xVals(:);
yVals = yVals(:);
n = numel(xVals);

xMin = min(xVals);
xMax = max(xVals);
yMin = min(yVals);
yMax = max(yVals);

xRange = xMax - xMin;
yRange = yMax - yMin;
if ~isfinite(xRange) || xRange < 0, xRange = 0; end
if ~isfinite(yRange) || yRange < 0, yRange = 0; end

xPad = max(0.10 * max(xRange, eps), 0.20 * estimate_bandwidth(xVals));
yPad = max(0.10 * max(yRange, eps), 0.20 * estimate_bandwidth(yVals));

xGrid = linspace(xMin - xPad, xMax + xPad, gridSize);
yGrid = linspace(yMin - yPad, yMax + yPad, gridSize);
[X, Y] = meshgrid(xGrid, yGrid);

hx = estimate_bandwidth(xVals);
hy = estimate_bandwidth(yVals);
if ~isfinite(hx) || hx <= 0 || ~isfinite(hy) || hy <= 0
    return;
end

Z = zeros(size(X));
for ii = 1:n
    ZX = ((X - xVals(ii)) ./ hx) .^ 2;
    ZY = ((Y - yVals(ii)) ./ hy) .^ 2;
    Z = Z + exp(-0.5 * (ZX + ZY));
end
Z = Z ./ (n * 2 * pi * hx * hy);

level = density_level_for_mass(Z, targetMass);
if ~isfinite(level) || level <= 0
    return;
end

[~, hContour] = contour(ax, xGrid, yGrid, Z, [level level], ...
    'LineColor', col, ...
    'LineWidth', 1.8, ...
    'HandleVisibility', 'off');
if isgraphics(hContour)
    hContour.HandleVisibility = 'off';
end
end


% =========================================================================
function bw = estimate_bandwidth(vals)
vals = vals(isfinite(vals));
n = numel(vals);
if n <= 1
    bw = eps;
    return;
end

sigma = std(vals);
iqrVal = local_percentile(vals, 0.75) - local_percentile(vals, 0.25);
sigmaRobust = iqrVal / 1.349;

if isfinite(sigmaRobust) && sigmaRobust > 0
    sigmaUse = min(sigma, sigmaRobust);
else
    sigmaUse = sigma;
end

if ~isfinite(sigmaUse) || sigmaUse <= 0
    sigmaUse = max(max(vals) - min(vals), eps) / 4;
end

bw = 1.06 * sigmaUse * n^(-1/5);
if ~isfinite(bw) || bw <= 0
    bw = max(max(vals) - min(vals), eps) / 20;
end
end


% =========================================================================
function level = density_level_for_mass(Z, targetMass)
targetMass = min(max(targetMass, 0.05), 0.95);
zVals = Z(isfinite(Z) & Z > 0);
if isempty(zVals)
    level = nan;
    return;
end

zVals = sort(zVals, 'descend');
cumMass = cumsum(zVals) / sum(zVals);
idx = find(cumMass >= targetMass, 1, 'first');
if isempty(idx)
    idx = numel(zVals);
end
level = zVals(idx);
end


% =========================================================================
function pctVal = local_percentile(vals, p)
vals = vals(isfinite(vals));
if isempty(vals)
    pctVal = nan;
    return;
end

vals = sort(vals(:));
p = min(max(p, 0), 1);
if numel(vals) == 1
    pctVal = vals;
    return;
end

idx = 1 + p * (numel(vals) - 1);
idxLo = floor(idx);
idxHi = ceil(idx);
w = idx - idxLo;
pctVal = (1 - w) * vals(idxLo) + w * vals(idxHi);
end


% =========================================================================
function val = get_opt_value(plotOpts, fieldName, defaultVal)
if nargin < 1 || ~isstruct(plotOpts) || ~isfield(plotOpts, fieldName) || isempty(plotOpts.(fieldName))
    val = defaultVal;
else
    val = plotOpts.(fieldName);
end
end


% =========================================================================
function style_legend_for_theme(leg, theme)
if isempty(leg) || ~isgraphics(leg), return; end
if strcmp(theme, 'poster')
    leg.TextColor = [1 1 1];
    leg.Color     = [0 0 0];
    leg.EdgeColor = [0.70 0.70 0.70];
else
    leg.TextColor = [0 0 0];
    leg.Color     = [1 1 1];
    leg.EdgeColor = [0.80 0.80 0.80];
end
end
