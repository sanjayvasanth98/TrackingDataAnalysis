%% test_plot_upstream_size_dist.m
% Standalone test for plot_upstream_size_distribution_by_re — upstream
% moving microbubble equivalent diameter PDF plot.
% Loads saved .mat data and calls the plot function with the same
% arguments as main_batch_trackmate_local.m.
%
% Edit plotOpts below to experiment, then reflect changes back to main code.

clear; clc;

%% Paths
matDir = "C:\Users\kbsanjayvasanth\Downloads\plot_data_mat";
outDir = fullfile(fileparts(mfilename('fullpath')), 'test_outputs', 'UpstreamSizeDistributions');
if ~isfolder(outDir), mkdir(outDir); end

addpath(fileparts(fileparts(mfilename('fullpath'))));

%% Load data
S = load(fullfile(matDir, "upstream_size_distribution_by_case.mat"), 'allSize');
allSize = S.allSize;
nCases = numel(allSize.caseName);
fprintf('Loaded %d cases from upstream_size_distribution_by_case.mat\n', nCases);
for i = 1:nCases
    nSamples = 0;
    if ~isempty(allSize.size_eqd{i})
        nSamples = numel(allSize.size_eqd{i});
    end
    fprintf('  %s: Re=%g, k/d=%.4g, %d eqd samples\n', ...
        allSize.caseName(i), allSize.Re(i), allSize.kD(i), nSamples);
end

%% Plot options (mirror main code — edit here to test changes)
set(0, 'DefaultFigureVisible', 'on');
plotOpts = struct();
plotOpts.enableNormalTheme = true;
plotOpts.enablePosterTheme = false;
plotOpts.savePNG = true;
plotOpts.saveSVG = false;
plotOpts.themes = "normal";
plotOpts.keepFiguresOpen = true;
plotOpts.upstreamSizeXLim_um = [];  % auto-compute; set [lo hi] to override

% binSize_phys passed as 3rd arg (same as main code)
binSize_phys = 0.02;

%% Generate plot
fprintf('\nGenerating upstream size distribution plot...\n');
plot_upstream_size_distribution_by_re(allSize, outDir, binSize_phys, plotOpts);
fprintf('Done. Output in: %s\n', outDir);
