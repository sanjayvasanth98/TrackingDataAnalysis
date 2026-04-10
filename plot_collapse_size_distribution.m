function plot_collapse_size_distribution(allCollapse, figDir, plotOpts)
%PLOT_COLLAPSE_SIZE_DISTRIBUTION  KDE of peak equivalent diameter at collapse.
%
%   Overlaid density curves, one per case. X-axis: peak equivalent
%   diameter d_eq (um) computed from peakArea_px2. Mirrors the style of
%   plot_upstream_size_distribution_by_re.

if nargin < 3 || ~isfield(plotOpts,'themes') || isempty(plotOpts.themes)
    plotOpts.themes = "normal";
end

nCases = numel(allCollapse.caseName);
if nCases == 0
    warning('plot_collapse_size_distribution: no cases. Skipping.');
    return;
end

hasPixelSize = isfield(allCollapse, 'pixelSize') && numel(allCollapse.pixelSize) >= nCases;

for theme = reshape(plotOpts.themes, 1, [])
    f = figure('Color','w','Position',[120 120 1000 700]);
    ax = axes(f);
    hold(ax,'on');

    cmap = lines(max(nCases,1));
    lgd    = gobjects(0,1);
    lgdTxt = strings(0,1);
    yMaxPlot = 0;
    anyPlotted = false;

    % Pool all diameters first to determine shared x-axis limits
    pooledX = [];
    for ci = 1:nCases
        cd = allCollapse.data{ci};
        if isempty(cd) || cd.nQualified < 5, continue; end
        areas = cd.peakArea_px2(:);
        areas = areas(isfinite(areas) & areas > 0);
        if hasPixelSize
            d_eq = 2 * sqrt(areas / pi) * allCollapse.pixelSize(ci) * 1000; % um
        else
            d_eq = 2 * sqrt(areas / pi); % px (fallback)
        end
        pooledX = [pooledX; d_eq]; %#ok<AGROW>
    end

    if isempty(pooledX)
        warning('plot_collapse_size_distribution: no cases with >= 5 events. Skipping.');
        close(f);
        continue;
    end

    xLimUse = robust_size_xlim(pooledX);

    for ci = 1:nCases
        cd = allCollapse.data{ci};
        if isempty(cd) || cd.nQualified < 5, continue; end

        areas = cd.peakArea_px2(:);
        areas = areas(isfinite(areas) & areas > 0);
        if hasPixelSize
            x = 2 * sqrt(areas / pi) * allCollapse.pixelSize(ci) * 1000; % um
        else
            x = 2 * sqrt(areas / pi); % px
        end

        if numel(x) < 5, continue; end

        % Remove extreme outliers (1st-99th percentile)
        pLo = percentile_linear(x, 1);
        pHi = percentile_linear(x, 99);
        if isfinite(pLo) && isfinite(pHi) && pHi > pLo
            x = x(x >= pLo & x <= pHi);
        end
        if numel(x) < 5, continue; end

        xi = linspace(xLimUse(1), xLimUse(2), 300);
        fhat = estimate_pdf_density(x, xi);

        col = cmap(ci,:);
        hLine = plot(ax, xi, fhat, '-', 'LineWidth', 2.0, 'Color', col);
        yMaxPlot = max(yMaxPlot, max(fhat));
        anyPlotted = true;

        lgd(end+1,1)    = hLine; %#ok<AGROW>
        lgdTxt(end+1,1) = sprintf('k/d=%.4g (n=%d)', ...
            allCollapse.kD(ci), cd.nQualified); %#ok<AGROW>
    end

    if ~anyPlotted
        warning('plot_collapse_size_distribution: nothing to plot after filtering.');
        close(f);
        continue;
    end

    xlim(ax, xLimUse);
    if yMaxPlot > 0
        ylim(ax, [0, yMaxPlot * 1.10]);
    end

    if hasPixelSize
        xlabel(ax, 'Peak equivalent diameter at collapse, $d_{eq}\,(\mu\mathrm{m})$', ...
            'Interpreter','latex');
    else
        xlabel(ax, 'Peak equivalent diameter at collapse, $d_{eq}\,(\mathrm{px})$', ...
            'Interpreter','latex');
    end
    ylabel(ax, 'PDF', 'Interpreter','latex');
    title(ax, '');
    grid(ax, 'off');
    box(ax, 'on');
    set(ax, 'LooseInset', max(get(ax,'LooseInset'), get(ax,'TightInset')));

    if ~isempty(lgd)
        leg = legend(ax, lgd, cellstr(lgdTxt), 'Location','northeast','Box','off');
    else
        leg = [];
    end

    apply_plot_theme(ax, char(theme));
    style_legend_for_theme(leg, char(theme));

    outBase = fullfile(figDir, "CollapseSizeDist_" + theme);
    save_fig_dual_safe(f, outBase, plotOpts);
    if ~isfield(plotOpts,'keepFiguresOpen') || ~plotOpts.keepFiguresOpen
        close(f);
    end
end
fprintf('Saved collapse size distribution plot to: %s\n', figDir);
end


% =========================================================================
function xLim = robust_size_xlim(x)
if isempty(x)
    xLim = [0 30];
    return;
end
x = x(isfinite(x));
if isempty(x)
    xLim = [0 30];
    return;
end
xMin = min(x);
xHi  = percentile_linear(x, 99);
if ~(isfinite(xHi) && xHi > xMin)
    xHi = max(x);
end
if ~(isfinite(xHi) && xHi > xMin)
    xHi = xMin + max(abs(xMin) * 0.1, 1);
end
pad = max(0.05 * (xHi - xMin), 5);
xLo = max(0, xMin - pad);
xLim = [xLo, xHi + pad];
end


% =========================================================================
function fhat = estimate_pdf_density(x, xi)
x  = x(:);
x  = x(isfinite(x));
xi = xi(:).';
n  = numel(x);
fhat = zeros(size(xi));
if n < 2 || isempty(xi), return; end

if exist('ksdensity','file') == 2
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
sx   = std(x, 0);
iqrx = iqr_linear(x);
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
p   = min(100, max(0, p));
idx = 1 + (n - 1) * (p / 100);
i0  = floor(idx);
i1  = ceil(idx);
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
w   = q75 - q25;
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
