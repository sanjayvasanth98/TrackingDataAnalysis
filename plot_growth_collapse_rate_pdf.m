function plot_growth_collapse_rate_pdf(allGrowthCollapse, figDir, plotOpts, matDir)
%PLOT_GROWTH_COLLAPSE_RATE_PDF  PDF of normalized axial growth/collapse rate.
%
%   Plots signed dL/dt/U distributions:
%     collapse: negative dL/dt/U, red
%     growth:   positive dL/dt/U, blue
%
%   Line style indicates the aspect-ratio group.

if nargin < 3 || ~isfield(plotOpts, 'themes') || isempty(plotOpts.themes)
    plotOpts.themes = "normal";
end
if nargin < 4
    matDir = "";
end
if ~isfolder(figDir), mkdir(figDir); end

if isempty(allGrowthCollapse) || ~isfield(allGrowthCollapse, 'data') || isempty(allGrowthCollapse.data)
    warning('plot_growth_collapse_rate_pdf: no growth/collapse data. Skipping.');
    return;
end

ReVals = unique(allGrowthCollapse.Re(:));
ReVals = ReVals(isfinite(ReVals));
if isempty(ReVals)
    warning('plot_growth_collapse_rate_pdf: no finite Reynolds values. Skipping.');
    return;
end

for ri = 1:numel(ReVals)
    ReNow = ReVals(ri);
    caseIdx = find(allGrowthCollapse.Re == ReNow);
    if isempty(caseIdx), continue; end

    plotData = build_pdf_plot_data(allGrowthCollapse, caseIdx, ReNow);
    if isempty(plotData.group)
        continue;
    end

    if matDir ~= ""
        if ~isfolder(matDir), mkdir(matDir); end
        reTag = sanitize_tag(sprintf('Re_%g', ReNow));
        save(fullfile(matDir, "growth_collapse_rate_pdf_plot_data_" + reTag + ".mat"), 'plotData');
    end

    xAll = [];
    for gi = 1:numel(plotData.group)
        xAll = [xAll; plotData.group(gi).growthRate_over_U(:); plotData.group(gi).collapseRate_over_U(:)]; %#ok<AGROW>
    end
    xAll = xAll(isfinite(xAll));
    if isempty(xAll)
        continue;
    end

    xLimAbs = percentile_linear(abs(xAll), 99);
    if ~(isfinite(xLimAbs) && xLimAbs > 0)
        xLimAbs = max(abs(xAll));
    end
    if ~(isfinite(xLimAbs) && xLimAbs > 0)
        xLimAbs = 1;
    end
    xLimAbs = min(max(xLimAbs, 0.25), max(abs(xAll)));
    xiCollapse = linspace(-xLimAbs, 0, 300);
    xiGrowth = linspace(0, xLimAbs, 300);

    for theme = reshape(plotOpts.themes, 1, [])
        f = figure('Color', 'w', 'Position', [120 120 1050 650]);
        ax = axes(f);
        hold(ax, 'on');

        growthColor = [0.00 0.27 0.80];
        collapseColor = [0.82 0.08 0.08];
        lineStyles = {'-', '--', ':', '-.'};
        lgd = gobjects(0,1);
        lgdTxt = strings(0,1);
        yMax = 0;
        yPositive = nan(0,1);

        for gi = 1:numel(plotData.group)
            ls = lineStyles{mod(gi - 1, numel(lineStyles)) + 1};
            label = char(string(plotData.group(gi).label));

            cVals = plotData.group(gi).collapseRate_over_U(:);
            cVals = cVals(isfinite(cVals) & cVals < 0 & cVals >= -xLimAbs);
            if numel(cVals) >= 2
                cPdf = estimate_pdf_density(cVals, xiCollapse);
                hC = plot(ax, xiCollapse, cPdf, ...
                    'Color', collapseColor, 'LineStyle', ls, 'LineWidth', 2.0);
                yMax = max(yMax, max(cPdf));
                yPositive = [yPositive; cPdf(cPdf > 0).']; %#ok<AGROW>
                lgd(end+1,1) = hC; %#ok<AGROW>
                lgdTxt(end+1,1) = sprintf('Collapse %s', label); %#ok<AGROW>
            end

            gVals = plotData.group(gi).growthRate_over_U(:);
            gVals = gVals(isfinite(gVals) & gVals > 0 & gVals <= xLimAbs);
            if numel(gVals) >= 2
                gPdf = estimate_pdf_density(gVals, xiGrowth);
                hG = plot(ax, xiGrowth, gPdf, ...
                    'Color', growthColor, 'LineStyle', ls, 'LineWidth', 2.0);
                yMax = max(yMax, max(gPdf));
                yPositive = [yPositive; gPdf(gPdf > 0).']; %#ok<AGROW>
                lgd(end+1,1) = hG; %#ok<AGROW>
                lgdTxt(end+1,1) = sprintf('Growth %s', label); %#ok<AGROW>
            end
        end

        yPositive = yPositive(isfinite(yPositive) & yPositive > 0);
        if isempty(yPositive)
            warning('plot_growth_collapse_rate_pdf: no PDF curves with at least two samples for Re=%g.', ReNow);
            close(f);
            continue;
        end

        xlim(ax, [-xLimAbs, xLimAbs]);
        if ~isempty(yPositive) && yMax > 0
            yMin = max(min(yPositive) * 0.8, yMax * 1e-4);
            yUpper = yMax * 1.20;
            set(ax, 'YScale', 'log');
            ylim(ax, [yMin, yUpper]);
            plot(ax, [0 0], [yMin, yUpper], '-', ...
                'Color', [0.70 0.70 0.70], 'LineWidth', 1.0, ...
                'HandleVisibility', 'off');
        end

        xlabel(ax, 'Axial rate, $dL/dt\,/\,U$', 'Interpreter', 'latex');
        ylabel(ax, 'P.d.f.', 'Interpreter', 'latex');
        title(ax, '');
        grid(ax, 'off');
        box(ax, 'on');

        text(ax, 0.25, 0.05, 'Collapse rate', ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'Color', collapseColor);
        text(ax, 0.75, 0.05, 'Growth rate', ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'Color', growthColor);

        if ~isempty(lgd)
            leg = legend(ax, lgd, cellstr(lgdTxt), ...
                'Location', 'northeast', 'Box', 'on', 'NumColumns', 1);
        else
            leg = [];
        end

        apply_plot_theme(ax, char(theme));
        style_legend_for_theme(leg, char(theme));

        outBase = fullfile(figDir, sprintf('GrowthCollapseRatePDF_Re_%g_%s', ReNow, theme));
        save_fig_dual_safe(f, outBase, plotOpts);
        if ~isfield(plotOpts, 'keepFiguresOpen') || ~plotOpts.keepFiguresOpen
            close(f);
        end
    end
end

fprintf('Saved growth/collapse rate PDF plots to: %s\n', figDir);
end


% =========================================================================
function plotData = build_pdf_plot_data(allGrowthCollapse, caseIdx, ReNow)
firstData = allGrowthCollapse.data{caseIdx(1)};
edges = firstData.aspectRatioEdges(:).';
labels = firstData.aspectRatioLabels(:);
nGroups = max(numel(edges) - 1, 0);

plotData = struct();
plotData.Re = ReNow;
plotData.caseName = string(allGrowthCollapse.caseName(caseIdx));
plotData.kD = allGrowthCollapse.kD(caseIdx);
plotData.U_m_s = nan(numel(caseIdx), 1);
plotData.aspectRatioEdges = edges;
plotData.aspectRatioLabels = labels;
plotData.group = repmat(struct( ...
    'label', "", ...
    'arMin', NaN, ...
    'arMax', NaN, ...
    'growthRate_over_U', nan(0,1), ...
    'collapseRate_over_U', nan(0,1), ...
    'nGrowth', 0, ...
    'nCollapse', 0), nGroups, 1);

for gi = 1:nGroups
    plotData.group(gi).label = labels(gi);
    plotData.group(gi).arMin = edges(gi);
    plotData.group(gi).arMax = edges(gi + 1);
end

for ii = 1:numel(caseIdx)
    data = allGrowthCollapse.data{caseIdx(ii)};
    if isfield(data, 'U_m_s')
        plotData.U_m_s(ii) = data.U_m_s;
    end

    for gi = 1:nGroups
        gMask = segment_ar_mask(data.growthSegments, edges(gi), edges(gi + 1));
        cMask = segment_ar_mask(data.collapseSegments, edges(gi), edges(gi + 1));

        plotData.group(gi).growthRate_over_U = [plotData.group(gi).growthRate_over_U; ...
            data.growthSegments.rate_over_U(gMask)]; %#ok<AGROW>
        plotData.group(gi).collapseRate_over_U = [plotData.group(gi).collapseRate_over_U; ...
            data.collapseSegments.rate_over_U(cMask)]; %#ok<AGROW>
    end
end

for gi = 1:nGroups
    plotData.group(gi).growthRate_over_U = plotData.group(gi).growthRate_over_U(:);
    plotData.group(gi).collapseRate_over_U = plotData.group(gi).collapseRate_over_U(:);
    plotData.group(gi).nGrowth = sum(isfinite(plotData.group(gi).growthRate_over_U));
    plotData.group(gi).nCollapse = sum(isfinite(plotData.group(gi).collapseRate_over_U));
end
end


% =========================================================================
function mask = segment_ar_mask(seg, arMin, arMax)
if isempty(seg) || ~isfield(seg, 'rate_over_U') || isempty(seg.rate_over_U)
    mask = false(0,1);
    return;
end
if isinf(arMax)
    arMask = seg.aspectRatio >= arMin;
else
    arMask = seg.aspectRatio >= arMin & seg.aspectRatio < arMax;
end
mask = isfinite(seg.rate_over_U) & isfinite(seg.aspectRatio) & arMask;
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
if n < 2, return; end
sx = std(x, 0);
iqrx = percentile_linear(x, 75) - percentile_linear(x, 25);
scale = min(sx, iqrx / 1.34);
if ~(isfinite(scale) && scale > 0)
    scale = max(sx, iqrx / 1.34);
end
if ~(isfinite(scale) && scale > 0), return; end
h = 0.9 * scale * (n ^ (-1/5));
end


% =========================================================================
function q = percentile_linear(x, p)
q = NaN;
if isempty(x) || ~isfinite(p), return; end
x = x(:);
x = x(isfinite(x));
if isempty(x), return; end
x = sort(x);
n = numel(x);
if n == 1, q = x(1); return; end
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
function tag = sanitize_tag(raw)
tag = string(raw);
tag = regexprep(tag, '[^\w]+', '_');
tag = regexprep(tag, '_+$', '');
end


% =========================================================================
function style_legend_for_theme(leg, theme)
if isempty(leg) || ~isgraphics(leg), return; end
if strcmp(theme, 'poster')
    leg.TextColor = [1 1 1];
    leg.Color = [0 0 0];
    leg.EdgeColor = [0.70 0.70 0.70];
else
    leg.TextColor = [0 0 0];
    leg.Color = [1 1 1];
    leg.EdgeColor = [0.80 0.80 0.80];
end
end
