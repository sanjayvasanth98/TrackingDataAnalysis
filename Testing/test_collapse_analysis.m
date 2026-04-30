%% test_collapse_analysis.m
% Standalone test for the collapse analysis pipeline.
% Loads saved .mat data from a completed main run and re-runs:
%   - plot_collapse_rate_vs_frame  (cumulative count vs time)
%   - plot_collapse_rate_vs_kd     (collapse rate vs k/d)
%   - plot_collapse_size_distribution (peak diameter KDE)
%   - write_collapse_analysis_csv
%
% Edit plotOpts below to experiment, then reflect
% any changes back to main_batch_trackmate_local.m.

clear; clc;

addpath(fileparts(mfilename('fullpath')));
addpath(fileparts(fileparts(mfilename('fullpath'))));

%% Paths
matDir   = test_plotmat_location("collapse_analysis_by_case.mat");
outDir   = fullfile(fileparts(mfilename('fullpath')), 'test_outputs', 'CollapseAnalysis');
csvOut   = fullfile(fileparts(mfilename('fullpath')), 'test_outputs', 'collapse_analysis.csv');
if ~isfolder(outDir), mkdir(outDir); end

%% Load allCollapse
S = load(fullfile(matDir, "collapse_analysis_by_case.mat"), 'allCollapse');
allCollapse = S.allCollapse;
nCases = numel(allCollapse.caseName);
fprintf('Loaded %d cases from collapse_analysis_by_case.mat\n\n', nCases);

% Backward compat: add pixelSize if missing from old .mat files
if ~isfield(allCollapse, 'pixelSize')
    allCollapse.pixelSize = repmat(0.00375009375, 1, nCases);
    fprintf('(pixelSize not in .mat — using default 0.00375 mm/px)\n\n');
end

%% Print per-case summary
fprintf('%-12s  %6s  %8s  %10s  %10s  %12s  %12s  %12s\n', ...
    'Case','Re','k/d','nTotal','nQualified','rate/frame','rate/ms','rate/s');
fprintf('%s\n', repmat('-',1,80));
for ci = 1:nCases
    cd = allCollapse.data{ci};
    if isempty(cd)
        fprintf('%-12s  (no data)\n', allCollapse.caseName{ci});
        continue;
    end
    fprintf('%-12s  %6g  %8.4g  %10d  %10d  %12.4g  %12.4g  %12.4g\n', ...
        allCollapse.caseName{ci}, ...
        allCollapse.Re(ci), ...
        allCollapse.kD(ci), ...
        cd.nTotal, ...
        cd.nQualified, ...
        cd.ratePerFrame, ...
        get_collapse_rate_per_ms(cd), ...
        cd.ratePerSec);
    if cd.nQualified > 0
        pa = cd.peakArea_px2(:);
        pa = pa(isfinite(pa) & pa > 0);
        if ~isempty(pa)
            fprintf('  Peak area (px2): mean=%.1f  median=%.1f  std=%.1f\n', ...
                mean(pa), median(pa), std(pa));
        end
    end
end
fprintf('\n');

%% Plot options — edit here to test changes
set(0, 'DefaultFigureVisible', 'on');
plotOpts = struct();
plotOpts.enableNormalTheme = true;
plotOpts.enablePosterTheme = false;
plotOpts.savePNG  = true;
plotOpts.saveSVG  = false;
plotOpts.themes   = "normal";
plotOpts.keepFiguresOpen = true;

%% Re-run plots
fprintf('Generating cumulative collapse count vs time plot...\n');
plot_collapse_rate_vs_frame(allCollapse, outDir, plotOpts);

fprintf('Generating collapse rate vs k/d plot...\n');
plot_collapse_rate_vs_kd(allCollapse, outDir, plotOpts);

fprintf('Generating collapse size distribution plot...\n');
plot_collapse_size_distribution(allCollapse, outDir, plotOpts);

fprintf('Writing CSV...\n');
write_collapse_analysis_csv(allCollapse, csvOut);

fprintf('\nDone. Output in:\n  %s\n  %s\n', outDir, csvOut);


%% ========================================================================
function ratePerMs = get_collapse_rate_per_ms(cd)
ratePerMs = NaN;
if isfield(cd, 'ratePerMs') && ~isempty(cd.ratePerMs) && isfinite(cd.ratePerMs)
    ratePerMs = cd.ratePerMs;
elseif isfield(cd, 'ratePerSec') && ~isempty(cd.ratePerSec) && isfinite(cd.ratePerSec)
    ratePerMs = cd.ratePerSec / 1000;
end
end
