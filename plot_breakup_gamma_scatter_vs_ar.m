function plot_breakup_gamma_scatter_vs_ar(breakupData, figDir, plotOpts, arLabel, matDir)
% PLOT_BREAKUP_GAMMA_SCATTER_VS_AR
%   Scatter of gamma = (x_child - x_parent)/d vs parent aspect ratio.
%   One colour per k/d case, lightly transparent markers, no edges.
%   A PDF panel is added above the scatter to show the marginal density of
%   parent AR values for each case. Binned dashed mean trend lines are
%   drawn after the scatter markers so they stay visually on top.
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

markerSize = round(get_opt_value(plotOpts, 'breakupARMarkerSize', 30));
if ~isfinite(markerSize) || markerSize <= 0
    markerSize = 30;
end
markerSize = max(markerSize, 5);
markerAlpha = get_opt_value(plotOpts, 'breakupARMarkerAlpha', 0.40);
if ~isfinite(markerAlpha)
    markerAlpha = 0.40;
end
markerAlpha = min(max(markerAlpha, 0), 1);
trendMaxBins = round(get_opt_value(plotOpts, 'breakupARTrendMaxBins', 12));
trendMaxBins = max(trendMaxBins, 2);
trendMinCount = round(get_opt_value(plotOpts, 'breakupARTrendMinCount', 5));
trendMinCount = max(trendMinCount, 2);

xLim = get_opt_value(plotOpts, 'breakupARXLim', [0, 4]);
yLim = get_opt_value(plotOpts, 'breakupARGammaYLim', [-0.5, 0.5]);
xLim = sanitize_limit_pair(xLim, [0, 4]);
yLim = sanitize_limit_pair(yLim, [-0.5, 0.5]);

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

% ---- Plot per theme ----
for theme = reshape(plotOpts.themes, 1, [])
    fontName = resolve_plot_font_name();
    f = figure('Color', 'w', 'Position', [100 100 1000 820]);
    mainL = 0.09;
    mainR = 0.06;
    mainB = 0.12;
    mainW = 1 - mainL - mainR;
    pdfH = 0.17;
    pdfGap = 0.02;
    mainH = 0.56;
    pdfB = mainB + mainH + pdfGap;

    axPdf = axes(f, 'Position', [mainL, pdfB, mainW, pdfH]); %#ok<LAXES>
    axMain = axes(f, 'Position', [mainL, mainB, mainW, mainH]); %#ok<LAXES>
    hold(axPdf, 'on');
    hold(axMain, 'on');
    trendX = cell(nPresent, 1);
    trendY = cell(nPresent, 1);
    trendKeep = false(nPresent, 1);
    trendHandles = gobjects(0, 1);
    pdfYMax = 0;

    lgd    = gobjects(0, 1);
    lgdTxt = strings(0, 1);

    for p = 1:nPresent
        ci  = casesPresent(p);
        sel = caseIdx == ci;
        col = cmap(p, :);
        xVals = allAR(sel);
        yVals = allGamma(sel);
        scatterMask = xVals >= xLim(1) & xVals <= xLim(2) & ...
                      yVals >= yLim(1) & yVals <= yLim(2);
        pdfMask = xVals >= xLim(1) & xVals <= xLim(2);

        nVisible = nnz(scatterMask);
        if nVisible > 0
            hPt = scatter(axMain, xVals(scatterMask), yVals(scatterMask), markerSize, col, 'filled', ...
                'MarkerFaceAlpha', markerAlpha, ...
                'MarkerEdgeColor', 'none', ...
                'Clipping', 'on');
            lgd(end+1, 1)    = hPt; %#ok<AGROW>
            lgdTxt(end+1, 1) = sprintf('k/d = %.2f  (n = %d)', ...
                breakupData(ci).kD, nVisible); %#ok<AGROW>
        end

        xPdf = xVals(pdfMask);
        if numel(xPdf) >= 2
            xi = linspace(xLim(1), xLim(2), 300);
            fhat = estimate_pdf_density(xPdf, xi);
            if any(isfinite(fhat))
                plot(axPdf, xi, fhat, '-', ...
                    'Color', col, ...
                    'LineWidth', 2.0, ...
                    'HandleVisibility', 'off', ...
                    'Clipping', 'on');
                pdfYMax = max(pdfYMax, max(fhat(isfinite(fhat))));
            end
        end

        [xTrend, yTrend] = build_binned_mean_trend( ...
            xVals(scatterMask), yVals(scatterMask), trendMaxBins, trendMinCount, xLim);
        if numel(xTrend) >= 2
            trendX{p} = xTrend;
            trendY{p} = yTrend;
            trendKeep(p) = true;
        end
    end

    % ---- Zero reference ----
    plot(axMain, xLim, [0 0], ':', ...
        'Color', [0.55 0.55 0.55], 'LineWidth', 1.0, ...
        'HandleVisibility', 'off');

    % ---- Mean trend lines ----
    for p = 1:nPresent
        if ~trendKeep(p)
            continue;
        end
        hTrend = plot(axMain, trendX{p}, trendY{p}, '--', ...
            'Color', cmap(p, :), ...
            'LineWidth', 2.4, ...
            'HandleVisibility', 'off', ...
            'Clipping', 'on');
        trendHandles(end+1, 1) = hTrend; %#ok<AGROW>
    end

    % ---- Axes ----
    xlim(axPdf, xLim);
    if pdfYMax > 0
        ylim(axPdf, [0, pdfYMax * 1.12]);
    else
        ylim(axPdf, [0, 1]);
    end
    set(axPdf, 'FontName', fontName);
    ylabel(axPdf, 'PDF', 'Interpreter', 'latex');
    title(axPdf, '');
    grid(axPdf, 'off');
    box(axPdf, 'on');

    xlim(axMain, xLim);
    ylim(axMain, yLim);
    set(axMain, 'FontName', fontName);

    xlabel(axMain, 'Parent aspect ratio, $\mathrm{AR}_\mathrm{parent}$', ...
        'Interpreter', 'latex');
    ylabel(axMain, '$\gamma = (x_\mathrm{child} - x_\mathrm{parent})\,/\,d$', ...
        'Interpreter', 'latex');
    if arLabel ~= ""
        title(axMain, sprintf('Parent AR $\\geq$ %s', ...
            strrep(char(arLabel), 'AR', '')), 'Interpreter', 'latex');
    else
        title(axMain, '');
    end
    grid(axMain, 'off');
    box(axMain, 'on');

    if ~isempty(lgd)
        leg = legend(axMain, lgd, cellstr(lgdTxt), ...
            'Location', 'best', 'NumColumns', 1, 'Box', 'on');
    else
        leg = [];
    end

    apply_plot_theme(axPdf, char(theme));
    apply_plot_theme(axMain, char(theme));
    set(axPdf, 'XTickLabel', []);

    for h = reshape(trendHandles, 1, [])
        try
            uistack(h, 'top');
        catch
        end
    end

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
function [xTrend, yTrend] = build_binned_mean_trend(xVals, yVals, maxTrendBins, minBinCount, xLim)
valid = isfinite(xVals) & isfinite(yVals);
xVals = xVals(valid);
yVals = yVals(valid);

if numel(xVals) < minBinCount
    xTrend = nan(0, 1);
    yTrend = nan(0, 1);
    return;
end

nBins = min(maxTrendBins, floor(numel(xVals) / minBinCount));
if nBins < 1
    xTrend = nan(0, 1);
    yTrend = nan(0, 1);
    return;
end

x0 = xLim(1);
x1 = xLim(2);
if ~isfinite(x0) || ~isfinite(x1) || x1 <= x0
    x0 = min(xVals);
    x1 = max(xVals);
end
if ~isfinite(x0) || ~isfinite(x1) || x1 <= x0
    xTrend = nan(0, 1);
    yTrend = nan(0, 1);
    return;
end

edges = linspace(x0, x1, nBins + 1);
xTrend = nan(0, 1);
yTrend = nan(0, 1);
for bi = 1:nBins
    if bi < nBins
        inBin = xVals >= edges(bi) & xVals < edges(bi + 1);
    else
        inBin = xVals >= edges(bi) & xVals <= edges(bi + 1);
    end
    if nnz(inBin) < minBinCount
        continue;
    end
    xTrend(end+1, 1) = mean(xVals(inBin)); %#ok<AGROW>
    yTrend(end+1, 1) = mean(yVals(inBin)); %#ok<AGROW>
end
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
function pair = sanitize_limit_pair(pair, fallbackPair)
if nargin < 2 || isempty(fallbackPair)
    fallbackPair = [0, 1];
end

pair = double(pair(:).');
if numel(pair) < 2 || any(~isfinite(pair(1:2)))
    pair = fallbackPair;
else
    pair = sort(pair(1:2));
end

if pair(2) <= pair(1)
    pair = fallbackPair;
end
end


% =========================================================================
function fhat = estimate_pdf_density(x, xi)
x = x(:);
x = x(isfinite(x));
xi = xi(:).';
n = numel(x);

fhat = zeros(size(xi));
if n < 2 || isempty(xi)
    return;
end

if exist('ksdensity', 'file') == 2
    try
        fhat = ksdensity(x, xi);
        fhat = fhat(:).';
        return;
    catch
    end
end

h = silverman_bandwidth(x);
if ~(isfinite(h) && h > 0)
    xSpan = max(xi) - min(xi);
    if ~(isfinite(xSpan) && xSpan > 0)
        xSpan = max(x) - min(x);
    end
    h = max(xSpan / 50, sqrt(eps));
end

normConst = 1 / (n * h * sqrt(2*pi));
blockSize = 5000;
for s = 1:blockSize:n
    e = min(s + blockSize - 1, n);
    xb = x(s:e);
    u = bsxfun(@minus, xi, xb) / h;
    fhat = fhat + sum(exp(-0.5 * (u .^ 2)), 1);
end
fhat = normConst * fhat;
end


% =========================================================================
function h = silverman_bandwidth(x)
n = numel(x);
h = NaN;
if n < 2
    return;
end

sx = std(x, 0);
iqrx = iqr_linear(x);
scale = min(sx, iqrx / 1.34);
if ~(isfinite(scale) && scale > 0)
    scale = max(sx, iqrx / 1.34);
end
if ~(isfinite(scale) && scale > 0)
    return;
end

h = 0.9 * scale * (n ^ (-1/5));
end


% =========================================================================
function q = percentile_linear(x, p)
q = NaN;
if isempty(x) || ~isfinite(p)
    return;
end

x = x(:);
x = x(isfinite(x));
if isempty(x)
    return;
end

x = sort(x);
n = numel(x);
if n == 1
    q = x(1);
    return;
end

p = min(100, max(0, p));
idx = 1 + (n - 1) * (p / 100);
i0 = floor(idx);
i1 = ceil(idx);
if i0 == i1
    q = x(i0);
else
    frac = idx - i0;
    q = x(i0) + frac * (x(i1) - x(i0));
end
end


% =========================================================================
function w = iqr_linear(x)
q75 = percentile_linear(x, 75);
q25 = percentile_linear(x, 25);
w = q75 - q25;
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
