%% test_growth_collapse_rate.m
% Standalone test for Figure-7-style growth/collapse axial rate PDFs.
% Loads saved .mat data from a completed main run and re-runs:
%   - plot_growth_collapse_rate_pdf
%
% Edit matDir below to point to either:
%   combined: resultsDir/plot_data_mat
%   chunk:    resultsDir/results individual/chunk_N/plot_data_mat

clear; clc;

%% Paths — edit matDir to point to your latest results run
matDir = "C:\Users\kbsanjayvasanth\Downloads\plot_data_mat";
outDir = fullfile(fileparts(mfilename('fullpath')), 'test_outputs', 'GrowthCollapseRate');
if ~isfolder(outDir), mkdir(outDir); end

addpath(fileparts(fileparts(mfilename('fullpath'))));

%% Load growth/collapse data
S = load(fullfile(matDir, "growth_collapse_rate_by_case.mat"));
if isfield(S, 'allGrowthCollapse')
    allGrowthCollapse = S.allGrowthCollapse;
elseif isfield(S, 'chunkGrowthCollapse')
    allGrowthCollapse = S.chunkGrowthCollapse;
else
    error('growth_collapse_rate_by_case.mat must contain allGrowthCollapse or chunkGrowthCollapse.');
end

nCases = numel(allGrowthCollapse.data);
fprintf('Loaded %d cases from growth_collapse_rate_by_case.mat\n\n', nCases);

%% Print per-case summary
fprintf('%-12s  %6s  %8s  %10s  %10s\n', ...
    'Case', 'Re', 'k/d', 'N growth', 'N collapse');
fprintf('%s\n', repmat('-', 1, 56));
for ci = 1:nCases
    data = allGrowthCollapse.data{ci};
    nGrowth = 0;
    nCollapse = 0;
    if isfield(data, 'growthSegments')
        nGrowth = numel(data.growthSegments.rate_over_U);
    end
    if isfield(data, 'collapseSegments')
        nCollapse = numel(data.collapseSegments.rate_over_U);
    end
    fprintf('%-12s  %6g  %8.4g  %10d  %10d\n', ...
        char(string(allGrowthCollapse.caseName(ci))), ...
        allGrowthCollapse.Re(ci), ...
        allGrowthCollapse.kD(ci), ...
        nGrowth, nCollapse);
end
fprintf('\n');

%% Plot options — edit here to test changes
set(0, 'DefaultFigureVisible', 'on');
plotOpts = struct();
plotOpts.enableNormalTheme = true;
plotOpts.enablePosterTheme = false;
plotOpts.savePNG = true;
plotOpts.saveSVG = false;
plotOpts.themes = "normal";
plotOpts.keepFiguresOpen = true;

%% Re-run plot
fprintf('Generating growth/collapse axial rate PDF plot...\n');
plot_growth_collapse_rate_pdf(allGrowthCollapse, outDir, plotOpts, outDir);

fprintf('\nDone. Output in:\n  %s\n', outDir);
