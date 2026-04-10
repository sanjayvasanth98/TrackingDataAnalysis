function plot_breakup_gamma_beeswarm_vs_kd(breakupData, figDir, plotOpts, arLabel, matDir)
% PLOT_BREAKUP_GAMMA_BEESWARM_VS_KD
%   Beeswarm plot of gamma = (x_child - x_parent)/d vs k/d.
%   One beeswarm column per unique k/d value.
%   Fully opaque markers with black edges. Median bar per group.
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

markerSize = get_opt_value(plotOpts, 'breakupKDMarkerSize', 40.5);

% ---- Pool events per unique k/d value ----
kD_raw   = [breakupData.kD];
kD_unique = sort(unique(kD_raw));
nGroups   = numel(kD_unique);

gammaByGroup = cell(nGroups, 1);
countByGroup = zeros(nGroups, 1);
for ci = 1:nCases
    ev = breakupData(ci).events;
    if isempty(ev), continue; end
    gIdx = find(kD_unique == breakupData(ci).kD, 1);
    gammaByGroup{gIdx} = [gammaByGroup{gIdx}; [ev.gamma].'];
end
for g = 1:nGroups
    countByGroup(g) = numel(gammaByGroup{g});
end

totalEvents = sum(countByGroup);
if totalEvents == 0, warning('No breakup events.'); return; end

% ---- Compute beeswarm offsets for each group ----
allGamma = vertcat(gammaByGroup{:});
yRange = max(allGamma) - min(allGamma);
if yRange == 0, yRange = 1; end
maxHalfWidth = 0.35;

offsetByGroup = cell(nGroups, 1);
medianByGroup = nan(nGroups, 1);
for g = 1:nGroups
    if countByGroup(g) == 0
        offsetByGroup{g} = [];
    else
        offsetByGroup{g} = beeswarm_offsets(gammaByGroup{g}, yRange, maxHalfWidth);
        medianByGroup(g) = median(gammaByGroup{g});
    end
end

% ---- Save .mat ----
if matDir ~= ""
    plotData = struct();
    plotData.kD_values     = kD_unique;
    plotData.gamma_groups  = gammaByGroup;
    plotData.n_per_group   = countByGroup;
    plotData.median_gamma  = medianByGroup;
    plotData.arLabel       = arLabel;
    if arLabel ~= ""
        matFile = fullfile(matDir, sprintf("breakup_gamma_vs_kd_%s.mat", arLabel));
    else
        matFile = fullfile(matDir, "breakup_gamma_vs_kd.mat");
    end
    save(matFile, 'plotData');
end

% ---- Colours ----
cmap = lines(max(nGroups, 1));

% ---- Y-axis limits: symmetric about 0 ----
yAbsMax = ceil(max(abs(allGamma)));
if yAbsMax == 0, yAbsMax = 1; end

% ---- Plot per theme ----
for theme = reshape(plotOpts.themes, 1, [])
    fontName = resolve_plot_font_name();
    f  = figure('Color', 'w', 'Position', [100 100 900 700]);
    ax = axes(f); %#ok<LAXES>
    hold(ax, 'on');

    lgd    = gobjects(0, 1);
    lgdTxt = strings(0, 1);

    for g = 1:nGroups
        if countByGroup(g) == 0, continue; end
        xPos = g + offsetByGroup{g};
        yPos = gammaByGroup{g};
        col  = cmap(g, :);

        hPt = scatter(ax, xPos, yPos, markerSize, col, 'filled', ...
            'MarkerFaceAlpha', 1.0, ...
            'MarkerEdgeColor', [0 0 0], ...
            'LineWidth', 1.2);
        lgd(end+1, 1)    = hPt; %#ok<AGROW>
        lgdTxt(end+1, 1) = sprintf('k/d = %.2f  (n = %d)', ...
            kD_unique(g), countByGroup(g)); %#ok<AGROW>
    end

    % ---- Median bars ----
    for g = 1:nGroups
        if countByGroup(g) < 2, continue; end
        plot(ax, g + [-0.30, 0.30], [medianByGroup(g), medianByGroup(g)], '-', ...
            'Color', [0.15 0.15 0.15], 'LineWidth', 2.5, ...
            'HandleVisibility', 'off');
    end

    % ---- Dashed line connecting medians ----
    hasMedian = isfinite(medianByGroup);
    if sum(hasMedian) >= 2
        gx = find(hasMedian);
        plot(ax, gx, medianByGroup(hasMedian), '--', ...
            'Color', [0.4 0.4 0.4], 'LineWidth', 1.2, ...
            'HandleVisibility', 'off');
    end

    % ---- Zero reference ----
    plot(ax, [0.5, nGroups + 0.5], [0 0], ':', ...
        'Color', [0.55 0.55 0.55], 'LineWidth', 1.0, ...
        'HandleVisibility', 'off');

    % ---- Axes ----
    xlim(ax, [0.5, nGroups + 0.5]);
    ylim(ax, [-yAbsMax, yAbsMax]);
    set(ax, 'XTick', 1:nGroups, ...
        'XTickLabel', arrayfun(@(v) sprintf('%.2f', v), kD_unique, ...
        'UniformOutput', false), ...
        'FontName', fontName);

    xlabel(ax, '$k/d$', 'Interpreter', 'latex');
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
        outBase = fullfile(figDir, "Breakup_gamma_beeswarm_kD_" + arLabel + "_" + theme);
    else
        outBase = fullfile(figDir, "Breakup_gamma_beeswarm_kD_" + theme);
    end
    save_fig_dual_safe(f, outBase, plotOpts);
    if ~isfield(plotOpts, 'keepFiguresOpen') || ~plotOpts.keepFiguresOpen
        close(f);
    end
end

fprintf('Saved breakup beeswarm (gamma vs k/d) to: %s\n', figDir);
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
function xOff = beeswarm_offsets(yVals, yRange, maxHalfWidth)
%BEESWARM_OFFSETS  Deterministic horizontal offsets for a beeswarm column.
%   Uses KDE to shape the swarm width and golden-ratio-based placement
%   for an even, non-overlapping distribution.
    n = numel(yVals);
    xOff = zeros(n, 1);
    if n <= 1, return; end

    % Silverman bandwidth for KDE
    s = std(yVals);
    if ~isfinite(s) || s <= 0, s = yRange / 4; end
    sortedY = sort(yVals);
    q25 = sortedY(max(1, round(n * 0.25)));
    q75 = sortedY(min(n, round(n * 0.75)));
    iqrVal = q75 - q25;
    bwBase = min(s, max(iqrVal / 1.34, eps));
    bw = 0.9 * bwBase * n^(-1/5);
    if ~isfinite(bw) || bw <= 0, bw = yRange / 10; end

    % KDE density at each point
    dens = zeros(n, 1);
    for i = 1:n
        dens(i) = mean(exp(-0.5 * ((yVals - yVals(i)) / bw).^2));
    end
    if max(dens) > 0
        dens = dens / max(dens);
    end

    % Golden-ratio deterministic placement (quasi-uniform, no randomness)
    [~, order] = sort(yVals);
    xSorted = zeros(n, 1);
    phi = (sqrt(5) - 1) / 2;  % golden ratio conjugate
    for ii = 1:n
        t = mod(ii * phi, 1) * 2 - 1;  % maps to [-1, 1] quasi-uniformly
        xSorted(ii) = t * dens(order(ii)) * maxHalfWidth;
    end
    xOff(order) = xSorted;
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
