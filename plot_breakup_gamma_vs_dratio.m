function plot_breakup_gamma_vs_dratio(breakupData, figDir, plotOpts, arLabel)
% PLOT_BREAKUP_GAMMA_VS_DRATIO
%   Scatter of gamma = (x_child - x_parent)/d_roughness vs d_child/d_parent.
%   One scatter point per child-parent pair, lightly transparent markers.
%   Per-case coloured robust trend line using quantile bins and medians.
%   The x-axis is clipped at a high percentile of dRatio and hidden
%   outliers are reported in a text file.
%   Zero-line reference at gamma=0. Y centred at 0 with symmetric limits.
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

clipPercentile = get_opt_value(plotOpts, 'breakupDRatioClipPercentile', 99);
markerAlpha    = get_opt_value(plotOpts, 'breakupDRatioMarkerAlpha', 0.22);
maxTrendBins   = round(get_opt_value(plotOpts, 'breakupDRatioTrendMaxBins', 8));
minBinCount    = round(get_opt_value(plotOpts, 'breakupDRatioTrendMinCount', 5));
maxTrendBins   = max(maxTrendBins, 1);
minBinCount    = max(minBinCount, 1);

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

validAll = isfinite(allDRatio) & isfinite(allGamma);
allDRatio = allDRatio(validAll);
allGamma  = allGamma(validAll);
if isempty(allDRatio)
    warning('No finite breakup events found across all cases.');
    return;
end

% X-axis: clip at a robust percentile so a few extreme outliers do not
% dominate the whole view. Outliers are still counted and reported.
xUpper = local_percentile(allDRatio, clipPercentile / 100);
if ~isfinite(xUpper)
    xUpper = max(allDRatio);
end
visibleAll = allDRatio <= xUpper;
if ~any(visibleAll)
    visibleAll = true(size(allDRatio));
    xUpper = max(allDRatio);
end
nTotalPts   = numel(allDRatio);
nClippedPts = sum(~visibleAll);

xVisible = allDRatio(visibleAll);
xMinVisible = min(xVisible);
xRangeVisible = xUpper - xMinVisible;
xPad = 0.04 * xRangeVisible;
if ~isfinite(xPad) || xPad <= 0
    xPad = 0.25;
end
xLim = [max(0, xMinVisible - xPad), xUpper];

% Y-axis: symmetric about 0.
yGammaVisible = allGamma(visibleAll);
yAbsMax = ceil(max(abs(yGammaVisible)));
if yAbsMax == 0, yAbsMax = 1; end
yLimPlot = [-yAbsMax, yAbsMax];

% Per-case colour (using lines colourmap).
cmap = lines(max(nCases, 1));

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
        validCase = isfinite(dRatio) & isfinite(gamma);
        dRatio = dRatio(validCase);
        gamma  = gamma(validCase);
        if isempty(dRatio), continue; end

        visibleCase = dRatio <= xUpper;
        col      = cmap(ci, :);

        if any(visibleCase)
            % --- Scatter (light transparency) — used for legend ---
            hPt = scatter(ax, dRatio(visibleCase), gamma(visibleCase), 34, col, 'filled', ...
                'MarkerFaceAlpha', markerAlpha, ...
                'MarkerEdgeColor', 'none');
            lgd(end+1,1)    = hPt; %#ok<AGROW>
            lgdTxt(end+1,1) = sprintf('k/d = %.2f', breakupData(ci).kD); %#ok<AGROW>
        end

        % --- Robust trend line: quantile bins + median gamma + min count ---
        [xTrend, gammaTrend] = build_quantile_trend( ...
            dRatio(visibleCase), gamma(visibleCase), maxTrendBins, minBinCount);
        if numel(xTrend) >= 2
            plot(ax, xTrend, gammaTrend, '-', ...
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

write_outlier_report(figDir, breakupData, xUpper, clipPercentile, arLabel);

fprintf('Saved breakup plot to: %s\n', figDir);
end


% =========================================================================
function [xTrend, yTrend] = build_quantile_trend(xVals, yVals, maxTrendBins, minBinCount)
valid = isfinite(xVals) & isfinite(yVals);
xVals = xVals(valid);
yVals = yVals(valid);

if numel(xVals) < 2 * minBinCount
    xTrend = nan(0,1);
    yTrend = nan(0,1);
    return;
end

nBins = min(maxTrendBins, floor(numel(xVals) / minBinCount));
if nBins < 1
    xTrend = nan(0,1);
    yTrend = nan(0,1);
    return;
end

qEdges = local_percentile(xVals, linspace(0, 1, nBins + 1));
qEdges = unique(qEdges(:).', 'stable');
if numel(qEdges) < 2
    xTrend = nan(0,1);
    yTrend = nan(0,1);
    return;
end

xTrend = nan(numel(qEdges) - 1, 1);
yTrend = nan(numel(qEdges) - 1, 1);
keep = false(numel(qEdges) - 1, 1);

for b = 1:(numel(qEdges) - 1)
    if b < numel(qEdges) - 1
        inBin = xVals >= qEdges(b) & xVals < qEdges(b+1);
    else
        inBin = xVals >= qEdges(b) & xVals <= qEdges(b+1);
    end
    if sum(inBin) < minBinCount
        continue;
    end
    xTrend(b) = median(xVals(inBin));
    yTrend(b) = median(yVals(inBin));
    keep(b) = true;
end

xTrend = xTrend(keep);
yTrend = yTrend(keep);
end


% =========================================================================
function pctVals = local_percentile(xVals, probs)
xVals = xVals(isfinite(xVals));
if isempty(xVals)
    pctVals = nan(size(probs));
    return;
end

xVals = sort(xVals(:));
probs = min(max(probs, 0), 1);
if numel(xVals) == 1
    pctVals = repmat(xVals, size(probs));
    return;
end

idx = 1 + probs(:) * (numel(xVals) - 1);
idxLo = floor(idx);
idxHi = ceil(idx);
w = idx - idxLo;
pctVals = (1 - w) .* xVals(idxLo) + w .* xVals(idxHi);
pctVals = reshape(pctVals, size(probs));
end


% =========================================================================
function write_outlier_report(figDir, breakupData, xUpper, clipPercentile, arLabel)
if arLabel ~= ""
    reportFile = fullfile(figDir, "Breakup_gamma_vs_dRatio_" + arLabel + "_outliers.txt");
else
    reportFile = fullfile(figDir, "Breakup_gamma_vs_dRatio_outliers.txt");
end

fid = fopen(reportFile, 'w');
if fid < 0
    warning('Could not write breakup outlier report: %s', reportFile);
    return;
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

totalPts = 0;
totalClipped = 0;
fprintf(fid, 'Breakup gamma vs d_child/d_parent outlier report\n');
fprintf(fid, '================================================\n');
fprintf(fid, 'Visible x-limit: %.6g\n', xUpper);
fprintf(fid, 'Percentile used: %.2f\n\n', clipPercentile);
fprintf(fid, '%-20s %10s %12s\n', 'Case', 'TotalPts', 'ClippedPts');
fprintf(fid, '%s\n', repmat('-', 1, 46));

for ci = 1:numel(breakupData)
    ev = breakupData(ci).events;
    if isempty(ev)
        nPts = 0;
        nClip = 0;
    else
        dRatio = [ev.dRatio].';
        dRatio = dRatio(isfinite(dRatio));
        nPts = numel(dRatio);
        nClip = sum(dRatio > xUpper);
    end
    totalPts = totalPts + nPts;
    totalClipped = totalClipped + nClip;
    fprintf(fid, '%-20s %10d %12d\n', char(string(breakupData(ci).caseName)), nPts, nClip);
end

fprintf(fid, '%s\n', repmat('-', 1, 46));
fprintf(fid, '%-20s %10d %12d\n', 'TOTAL', totalPts, totalClipped);
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
