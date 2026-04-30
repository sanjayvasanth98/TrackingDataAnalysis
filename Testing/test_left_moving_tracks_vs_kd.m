%% test_left_moving_tracks_vs_kd.m
% Standalone test plot for left-moving track count vs k/d.
%
% Preferred input:
%   resultsDir\plot_data_mat\activation_summary_by_case.mat
%
% This .mat file stores summaryRows, including nLeftMovingTracks. That
% column is the strict left-moving/primary track count used as the current
% denominator for A/I.

clear; clc;

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot);
addpath(fileparts(mfilename('fullpath')));

%% Paths
% Option 1: point this to your plot_data_mat folder from a completed run.
matDir = test_plotmat_location("activation_summary_by_case.mat");

% Option 2: or point directly to activation_summary_by_case.mat.
% If non-empty, this takes priority over matDir.
matFile = "";

outDir = fullfile(fileparts(mfilename('fullpath')), ...
    'test_outputs', 'LeftMovingTracksVsKD');
if ~isfolder(outDir), mkdir(outDir); end

%% Plot options
set(0, 'DefaultFigureVisible', 'on');
plotOpts = struct();
plotOpts.enableNormalTheme = true;
plotOpts.enablePosterTheme = false;
plotOpts.savePNG = true;
plotOpts.saveSVG = false;
plotOpts.themes = "normal";
plotOpts.keepFiguresOpen = true;

countField = "nLeftMovingTracks";
countErrorMode = "sqrtPoisson"; % "sqrtPoisson" uses +/- sqrt(N); "none" hides error bars

%% Load data
summaryRows = load_activation_summary_rows(matDir, matFile);
requiredCols = {'Case', 'Re', 'kD', char(countField)};
missingCols = requiredCols(~ismember(requiredCols, summaryRows.Properties.VariableNames));
if ~isempty(missingCols)
    error('Summary table is missing required column(s): %s', strjoin(missingCols, ', '));
end

summaryRows = sortrows(summaryRows, {'Re', 'kD'});
fprintf('Loaded %d cases for left-moving track count plot.\n', height(summaryRows));
disp(summaryRows(:, requiredCols));

plotData = summaryRows(:, requiredCols);
[errLow, errHigh] = count_error_bars(plotData.(char(countField)), countErrorMode);
plotData.nLeftMovingTracks_err_low = errLow;
plotData.nLeftMovingTracks_err_high = errHigh;
save(fullfile(outDir, 'left_moving_tracks_vs_kd_plot_data.mat'), ...
    'plotData', 'countField', 'countErrorMode');
write_table_csv_compat(plotData, fullfile(outDir, 'left_moving_tracks_vs_kd_plot_data.csv'));

%% Generate plot
fprintf('\nGenerating left-moving tracks vs k/d plot...\n');
plot_left_moving_tracks_vs_kd(summaryRows, countField, countErrorMode, outDir, plotOpts);
fprintf('Done. Output in: %s\n', outDir);


% =========================================================================
function summaryRows = load_activation_summary_rows(matDir, matFile)
matDir = string(matDir);
matFile = string(matFile);

if strlength(strtrim(matFile)) == 0
    matFile = fullfile(matDir, "activation_summary_by_case.mat");
end

if ~isfile(matFile)
    error(['Could not find activation summary MAT file:\n  %s\n' ...
        'Set matDir to a folder containing activation_summary_by_case.mat, ' ...
        'or set matFile directly.'], matFile);
end

S = load(matFile);
if ~isfield(S, 'summaryRows')
    error('MAT file does not contain summaryRows: %s', matFile);
end
summaryRows = S.summaryRows;
if ~istable(summaryRows)
    error('summaryRows in %s is not a MATLAB table.', matFile);
end
fprintf('Loaded: %s\n', matFile);
end


% =========================================================================
function plot_left_moving_tracks_vs_kd(summaryRows, countField, countErrorMode, outDir, plotOpts)
countField = char(countField);
ReVals = unique(summaryRows.Re(:));
ReVals = ReVals(isfinite(ReVals));
if isempty(ReVals)
    warning('No finite Reynolds values found. Skipping plot.');
    return;
end

for theme = reshape(plotOpts.themes, 1, [])
    fontName = resolve_plot_font_name();
    f = figure('Color', 'w', 'Position', [120 120 980 680]);
    ax = axes(f);
    hold(ax, 'on');

    cmap = lines(max(numel(ReVals), 1));
    markerStyles = {'o', 's', '^', 'd', 'v', '>', '<', 'p', 'h'};
    lgd = gobjects(0, 1);
    lgdTxt = strings(0, 1);

    for r = 1:numel(ReVals)
        Rei = ReVals(r);
        sub = summaryRows(summaryRows.Re == Rei, :);
        sub = sortrows(sub, 'kD');

        x = sub.kD(:);
        y = sub.(countField)(:);
        [errLow, errHigh] = count_error_bars(y, countErrorMode);
        valid = isfinite(x) & isfinite(y);
        x = x(valid);
        y = y(valid);
        errLow = errLow(valid);
        errHigh = errHigh(valid);
        if isempty(x)
            continue;
        end

        col = cmap(r, :);
        markerStyle = markerStyles{mod(r - 1, numel(markerStyles)) + 1};
        if strcmpi(char(countErrorMode), 'none')
            h = plot(ax, x, y, '-', ...
                'Color', col, ...
                'LineWidth', 1.4, ...
                'Marker', markerStyle, ...
                'MarkerSize', 10, ...
                'MarkerFaceColor', col, ...
                'MarkerEdgeColor', [0 0 0]);
        else
            h = errorbar(ax, x, y, errLow, errHigh, markerStyle, ...
                'LineStyle', '-', ...
                'CapSize', 8, ...
                'Color', col, ...
                'LineWidth', 1.4, ...
                'MarkerSize', 10, ...
                'MarkerFaceColor', col, ...
                'MarkerEdgeColor', [0 0 0]);
        end
        set(h, ...
            'Color', col, ...
            'LineWidth', 1.4, ...
            'MarkerSize', 10, ...
            'MarkerFaceColor', col, ...
            'MarkerEdgeColor', [0 0 0]);
        lgd(end+1,1) = h; %#ok<AGROW>
        lgdTxt(end+1,1) = sprintf('Re=%g', Rei); %#ok<AGROW>
    end

    xlabel(ax, '$k/d$', 'Interpreter', 'latex');
    ylabel(ax, 'Number of left-moving tracks', 'Interpreter', 'latex');
    title(ax, 'Left-moving track count vs roughness', 'FontName', fontName);
    set(ax, 'FontName', fontName, 'YScale', 'linear');
    grid(ax, 'off');
    box(ax, 'on');

    yl = ylim(ax);
    if isfinite(yl(2)) && yl(2) > 0
        ylim(ax, [0, yl(2) * 1.08]);
    end

    if ~isempty(lgd)
        leg = legend(ax, lgd, cellstr(lgdTxt), ...
            'Location', 'northeast', 'Box', 'off');
    else
        leg = [];
    end

    apply_plot_theme(ax, char(theme));
    style_legend_for_theme(leg, char(theme));

    outBase = fullfile(outDir, "LeftMovingTracks_vs_kD_by_Re_" + theme);
    save_fig_dual_safe(f, outBase, plotOpts);
    if ~isfield(plotOpts, 'keepFiguresOpen') || ~plotOpts.keepFiguresOpen
        close(f);
    end
end
end


% =========================================================================
function [errLow, errHigh] = count_error_bars(counts, countErrorMode)
counts = double(counts(:));
errLow = zeros(size(counts));
errHigh = zeros(size(counts));

mode = lower(char(string(countErrorMode)));
switch mode
    case {'none', 'off'}
        return;
    case {'sqrtpoisson', 'sqrt', 'poisson'}
        sigma = sqrt(max(counts, 0));
        errLow = min(sigma, max(counts, 0));
        errHigh = sigma;
    otherwise
        warning('Unknown countErrorMode "%s"; using sqrtPoisson.', char(string(countErrorMode)));
        sigma = sqrt(max(counts, 0));
        errLow = min(sigma, max(counts, 0));
        errHigh = sigma;
end

errLow(~isfinite(errLow)) = 0;
errHigh(~isfinite(errHigh)) = 0;
end


% =========================================================================
function style_legend_for_theme(leg, theme)
if isempty(leg) || ~isgraphics(leg)
    return;
end
if strcmp(theme, 'poster')
    leg.TextColor = [1 1 1];
    leg.Color = [0 0 0];
else
    leg.TextColor = [0 0 0];
    leg.Color = 'none';
end
end
