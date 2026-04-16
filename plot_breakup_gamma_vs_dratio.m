function plot_breakup_gamma_vs_dratio(breakupData, figDir, plotOpts, arLabel)
% PLOT_BREAKUP_GAMMA_VS_DRATIO
%   Scatter of gamma = (x_child - x_parent)/d_roughness vs d_child/d_parent.
%   One scatter point per child-parent pair, lightly transparent markers.
%   Per-case coloured dashed binned mean trend line.
%   The x-axis defaults to log scale with robust low/high clipping, and
%   extreme tails are reported in a text file.
%   The y-axis is linear and fixed to the published breakup view window.
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

axisScale      = lower(string(get_opt_value(plotOpts, 'breakupDRatioXScale', "log")));
fixedXLim      = get_opt_value(plotOpts, 'breakupDRatioXLim', []);
fixedYLim      = get_opt_value(plotOpts, 'breakupGammaYLim', []);
yClipPercentile = get_opt_value(plotOpts, 'breakupGammaYClipPercentile', [1, 99]);
clipPercentile = get_opt_value(plotOpts, 'breakupDRatioClipPercentile', 99.5);
lowPercentile  = get_opt_value(plotOpts, 'breakupDRatioClipLowPercentile', 0.5);
markerAlpha    = get_opt_value(plotOpts, 'breakupDRatioMarkerAlpha', 0.10);
markerSize     = get_opt_value(plotOpts, 'breakupDRatioMarkerSize', 9);
maxTrendBins   = round(get_opt_value(plotOpts, 'breakupDRatioTrendMaxBins', 12));
minBinCount    = round(get_opt_value(plotOpts, 'breakupDRatioTrendMinCount', 5));
if axisScale ~= "log" && axisScale ~= "linear"
    axisScale = "log";
end
maxTrendBins   = max(maxTrendBins, 1);
minBinCount    = max(minBinCount, 1);
useFixedXLim = isnumeric(fixedXLim) && numel(fixedXLim) >= 2 && all(isfinite(fixedXLim(1:2)));
if useFixedXLim
    fixedXLim = sort(double(fixedXLim(1:2)));
    if fixedXLim(1) >= fixedXLim(2) || (axisScale == "log" && fixedXLim(1) <= 0)
        warning('Ignoring invalid breakupDRatioXLim. Falling back to percentile x-limits.');
        useFixedXLim = false;
    end
end
useFixedYLim = isnumeric(fixedYLim) && numel(fixedYLim) >= 2 && all(isfinite(fixedYLim(1:2)));
if useFixedYLim
    fixedYLim = sort(double(fixedYLim(1:2)));
    if fixedYLim(1) >= fixedYLim(2)
        warning('Ignoring invalid breakupGammaYLim. Falling back to robust y-limits.');
        useFixedYLim = false;
    end
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

validAll = isfinite(allDRatio) & isfinite(allGamma);
allDRatio = allDRatio(validAll);
allGamma  = allGamma(validAll);
if isempty(allDRatio)
    warning('No finite breakup events found across all cases.');
    return;
end
if axisScale == "log"
    positiveMask = allDRatio > 0;
    if ~any(positiveMask)
        warning('No positive dRatio values found for log-scale breakup plot.');
        return;
    end
    allDRatio = allDRatio(positiveMask);
    allGamma  = allGamma(positiveMask);
end

% X-axis: use a fixed requested window when supplied; otherwise clip at
% robust percentiles so extreme tails do not dominate the view.
if useFixedXLim
    xLower = fixedXLim(1);
    xUpper = fixedXLim(2);
else
    xUpper = local_percentile(allDRatio, clipPercentile / 100);
    if ~isfinite(xUpper)
        xUpper = max(allDRatio);
    end
    if axisScale == "log"
        xLower = local_percentile(allDRatio, lowPercentile / 100);
        if ~isfinite(xLower) || xLower <= 0
            xLower = min(allDRatio(allDRatio > 0));
        end
    else
        xLower = min(allDRatio);
    end
end

visibleAll = allDRatio >= xLower & allDRatio <= xUpper;
if ~any(visibleAll) && ~useFixedXLim
    visibleAll = true(size(allDRatio));
    xLower = min(allDRatio);
    xUpper = max(allDRatio);
end
xVisible = allDRatio(visibleAll);
if useFixedXLim
    xLim = [xLower, xUpper];
elseif axisScale == "log"
    xLowerPlot = max(xLower, min(xVisible(xVisible > 0)));
    logLo = log10(xLowerPlot);
    logHi = log10(xUpper);
    logPad = 0.04 * max(logHi - logLo, eps);
    xLim = 10.^[logLo - logPad, logHi + logPad];
else
    xMinVisible = min(xVisible);
    xRangeVisible = xUpper - xMinVisible;
    xPad = 0.04 * xRangeVisible;
    if ~isfinite(xPad) || xPad <= 0
        xPad = 0.25;
    end
    xLim = [max(0, xMinVisible - xPad), xUpper];
end

% Y-axis: robust limits avoid letting a few outliers leave the panel empty.
yGammaVisible = allGamma(visibleAll);
if isempty(yGammaVisible)
    yGammaVisible = allGamma;
end
if useFixedYLim
    yLimPlot = fixedYLim;
else
    if ~isnumeric(yClipPercentile) || numel(yClipPercentile) < 2 || any(~isfinite(yClipPercentile(1:2)))
        yClipPercentile = [1, 99];
    end
    yClipPercentile = sort(double(yClipPercentile(1:2)));
    yClipPercentile = min(max(yClipPercentile, 0), 100);
    yLoHi = local_percentile(yGammaVisible, yClipPercentile ./ 100);
    yLo = min(yLoHi(1), 0);
    yHi = max(yLoHi(2), 0);
    yRange = yHi - yLo;
    if ~isfinite(yRange) || yRange <= 0
        yRange = max(abs(yGammaVisible));
        if ~isfinite(yRange) || yRange <= 0
            yRange = 1;
        end
        yLo = -0.1 * yRange;
        yHi = 0.1 * yRange;
    end
    yPad = 0.08 * (yHi - yLo);
    yLimPlot = [yLo - yPad, yHi + yPad];
end

% Keep the published breakup view fixed so the panel is consistent run-to-run.
yLimPlot = [-0.5, 0.5];

% Per-case colour (using lines colourmap).
cmap = lines(max(nCases, 1));

for theme = reshape(plotOpts.themes, 1, [])
    f = figure('Color', 'w', 'Position', [100 100 900 900]);
    ax = axes(f);
    hold(ax, 'on');
    if axisScale == "log"
        set(ax, 'XScale', 'log');
        ax.XMinorTick = 'on';
    end

    lgd    = gobjects(0,1);
    lgdTxt = strings(0,1);
    trendHandles = gobjects(0,1);
    trendX = cell(nCases, 1);
    trendY = cell(nCases, 1);
    trendKeep = false(nCases, 1);

    for ci = 1:nCases
        ev = breakupData(ci).events;
        if isempty(ev), continue; end

        dRatio   = [ev.dRatio].';
        gamma    = [ev.gamma].';
        validCase = isfinite(dRatio) & isfinite(gamma);
        dRatio = dRatio(validCase);
        gamma  = gamma(validCase);
        if isempty(dRatio), continue; end
        if axisScale == "log"
            keepPositive = dRatio > 0;
            dRatio = dRatio(keepPositive);
            gamma  = gamma(keepPositive);
            if isempty(dRatio), continue; end
        end

        visibleCase = dRatio >= xLower & dRatio <= xUpper & ...
                      gamma  >= yLimPlot(1) & gamma  <= yLimPlot(2);
        col      = cmap(ci, :);

        if any(visibleCase)
            % Scatter with low opacity so dense clusters read as density.
            hPt = scatter(ax, dRatio(visibleCase), gamma(visibleCase), markerSize, col, 'filled', ...
                'MarkerFaceAlpha', markerAlpha, ...
                'MarkerEdgeColor', 'none', ...
                'Clipping', 'on');
            lgd(end+1,1)    = hPt; %#ok<AGROW>
            lgdTxt(end+1,1) = sprintf('k/d = %.2f', breakupData(ci).kD); %#ok<AGROW>
        end

        % Binned mean trend line matched to the visible x-window.
        [xTrend, gammaTrend] = build_binned_mean_trend( ...
            dRatio(visibleCase), gamma(visibleCase), maxTrendBins, minBinCount, axisScale, xLim);
        if numel(xTrend) >= 2
            trendX{ci} = xTrend;
            trendY{ci} = gammaTrend;
            trendKeep(ci) = true;
        end
    end

    % Zero reference line.
    plot(ax, xLim, [0 0], '--', 'Color', [0.5 0.5 0.5], ...
        'LineWidth', 1.2, 'HandleVisibility', 'off', 'Clipping', 'on');

    % Draw all trend lines after the scatters so they sit on top.
    for ci = 1:nCases
        if ~trendKeep(ci)
            continue;
        end
        hTrend = plot(ax, trendX{ci}, trendY{ci}, '--', ...
            'Color', cmap(ci, :), ...
            'LineWidth', 2.2, ...
            'HandleVisibility', 'off', ...
            'Clipping', 'on');
        trendHandles(end+1,1) = hTrend; %#ok<AGROW>
    end

    xlim(ax, xLim);
    ylim(ax, yLimPlot);

    enhance_minor_ticks(ax, axisScale, xLim, yLimPlot);

    xlabel(ax, '$d_\mathrm{child}/d_\mathrm{parent}$', 'Interpreter', 'latex');
    ylabel(ax, '$\gamma = (x_\mathrm{child} - x_\mathrm{parent})\,/\,d$', ...
        'Interpreter', 'latex');
    grid(ax, 'off');
    box(ax, 'on');

    if ~isempty(lgd)
        leg = legend(ax, lgd, cellstr(lgdTxt), ...
            'Location', 'southoutside', ...
            'Orientation', 'horizontal', ...
            'Box', 'off');
        try
            leg.NumColumns = max(1, numel(lgd));
        catch
            % Older MATLAB releases ignore NumColumns; southoutside still works.
        end
    else
        leg = [];
    end

    apply_plot_theme(ax, char(theme));

    % Re-assert the stacking order after theming so mean lines stay above
    % the scatter markers even when the renderer/theme updates children.
    for h = reshape(trendHandles, 1, [])
        try
            uistack(h, 'top');
        catch
        end
    end

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

write_outlier_report(figDir, breakupData, xLower, xUpper, lowPercentile, clipPercentile, axisScale, useFixedXLim, arLabel);

fprintf('Saved breakup plot to: %s\n', figDir);
end


% =========================================================================
function enhance_minor_ticks(ax, axisScale, xLim, yLimPlot)
if ~isgraphics(ax)
    return;
end

try
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
catch
end

try
    drawnow;
catch
end

% Add denser minor tick values by subdividing between major ticks.
try
    xMajor = xticks(ax);
    xMajor = xMajor(isfinite(xMajor));
    if numel(xMajor) >= 2
        xMinor = [];
        for k = 1:(numel(xMajor) - 1)
            x0 = xMajor(k);
            x1 = xMajor(k + 1);
            if ~isfinite(x0) || ~isfinite(x1) || x1 <= x0
                continue;
            end
            if axisScale == "log"
                if x0 <= 0 || x1 <= 0
                    continue;
                end
                % Use more internal positions than MATLAB's default 2..9
                % minor ticks so the crowded left side reads more clearly.
                vals = logspace(log10(x0), log10(x1), 13);
            else
                vals = linspace(x0, x1, 8);
            end
            xMinor = [xMinor, vals(2:end-1)]; %#ok<AGROW>
        end
        xMinor = unique(xMinor);
        xMinor = xMinor(isfinite(xMinor) & xMinor > xLim(1) & xMinor < xLim(2));
        if ~isempty(xMinor)
            ax.XAxis.MinorTickValues = xMinor;
        end
    end
catch
    % If the MATLAB release does not support custom minor tick values,
    % the automatic minor ticks still remain enabled above.
end

try
    yMajor = yticks(ax);
    yMajor = yMajor(isfinite(yMajor));
    if numel(yMajor) >= 2
        yMinor = [];
        for k = 1:(numel(yMajor) - 1)
            y0 = yMajor(k);
            y1 = yMajor(k + 1);
            if ~isfinite(y0) || ~isfinite(y1) || y1 <= y0
                continue;
            end
            vals = linspace(y0, y1, 7);
            yMinor = [yMinor, vals(2:end-1)]; %#ok<AGROW>
        end
        yMinor = unique(yMinor);
        yMinor = yMinor(isfinite(yMinor) & yMinor > yLimPlot(1) & yMinor < yLimPlot(2));
        if ~isempty(yMinor)
            ax.YAxis.MinorTickValues = yMinor;
        end
    end
catch
    % Older releases may not expose settable minor tick values on rulers.
end
end


% =========================================================================
function [xTrend, yTrend] = build_binned_mean_trend(xVals, yVals, maxTrendBins, minBinCount, axisScale, xLim)
valid = isfinite(xVals) & isfinite(yVals);
xVals = xVals(valid);
yVals = yVals(valid);
if axisScale == "log"
    positiveMask = xVals > 0;
    xVals = xVals(positiveMask);
    yVals = yVals(positiveMask);
end

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

if axisScale == "log"
    x0 = max(xLim(1), min(xVals));
    x1 = min(xLim(2), max(xVals));
    if ~isfinite(x0) || ~isfinite(x1) || x0 <= 0 || x0 >= x1
        x0 = min(xVals);
        x1 = max(xVals);
    end
    qEdges = logspace(log10(x0), log10(x1), nBins + 1);
else
    x0 = max(xLim(1), min(xVals));
    x1 = min(xLim(2), max(xVals));
    if ~isfinite(x0) || ~isfinite(x1) || x0 >= x1
        x0 = min(xVals);
        x1 = max(xVals);
    end
    qEdges = linspace(x0, x1, nBins + 1);
end
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
    xTrend(b) = mean(xVals(inBin));
    yTrend(b) = mean(yVals(inBin));
    keep(b) = true;
end

xTrend = xTrend(keep);
yTrend = yTrend(keep);
end


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
function write_outlier_report(figDir, breakupData, xLower, xUpper, lowPercentile, clipPercentile, axisScale, useFixedXLim, arLabel)
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
totalClippedLow = 0;
totalClippedHigh = 0;
fprintf(fid, 'Breakup gamma vs d_child/d_parent outlier report\n');
fprintf(fid, '================================================\n');
fprintf(fid, 'X scale: %s\n', char(axisScale));
fprintf(fid, 'Visible x-window: %.6g to %.6g\n', xLower, xUpper);
if useFixedXLim
    fprintf(fid, 'Limit mode: fixed breakupDRatioXLim\n\n');
else
    fprintf(fid, 'Limit mode: percentile\n');
    fprintf(fid, 'Low percentile used: %.2f\n', lowPercentile);
    fprintf(fid, 'High percentile used: %.2f\n\n', clipPercentile);
end
fprintf(fid, '%-20s %10s %12s %13s\n', 'Case', 'TotalPts', 'ClippedLow', 'ClippedHigh');
fprintf(fid, '%s\n', repmat('-', 1, 62));

for ci = 1:numel(breakupData)
    ev = breakupData(ci).events;
    if isempty(ev)
        nPts = 0;
        nClipLow = 0;
        nClipHigh = 0;
    else
        dRatio = [ev.dRatio].';
        dRatio = dRatio(isfinite(dRatio));
        if axisScale == "log"
            dRatio = dRatio(dRatio > 0);
        end
        nPts = numel(dRatio);
        nClipLow = sum(dRatio < xLower);
        nClipHigh = sum(dRatio > xUpper);
    end
    totalPts = totalPts + nPts;
    totalClippedLow = totalClippedLow + nClipLow;
    totalClippedHigh = totalClippedHigh + nClipHigh;
    fprintf(fid, '%-20s %10d %12d %13d\n', char(string(breakupData(ci).caseName)), nPts, nClipLow, nClipHigh);
end

fprintf(fid, '%s\n', repmat('-', 1, 62));
fprintf(fid, '%-20s %10d %12d %13d\n', 'TOTAL', totalPts, totalClippedLow, totalClippedHigh);
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
