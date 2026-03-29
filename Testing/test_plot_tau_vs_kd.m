%% test_plot_tau_vs_kd.m
% Standalone test for plot_tau_vs_kdh_re — Tau vs k/d plot.
% Loads saved .mat data and calls the plot function with the same
% arguments as main_batch_trackmate_local.m.
%
% Edit plotOpts below to experiment, then reflect changes back to main code.

clear; clc;

%% Paths
matDir = "E:\March Re 90,000 inception data\Processed images\results\results 28\plot_data_mat";
outDir = fullfile(fileparts(mfilename('fullpath')), 'test_outputs');
if ~isfolder(outDir), mkdir(outDir); end

addpath(fileparts(fileparts(mfilename('fullpath'))));

%% Load data
S = load(fullfile(matDir, "activation_summary_by_case.mat"), 'summaryRows');
summaryRows = S.summaryRows;
fprintf('Loaded %d cases from activation_summary_by_case.mat\n', height(summaryRows));
disp(summaryRows(:, {'Case','Re','kD','tau_mean','tau_std','tau_sem','nTau'}));

%% Plot options (mirror main code — edit here to test changes)
set(0, 'DefaultFigureVisible', 'on');
plotOpts = struct();
plotOpts.enableNormalTheme = true;
plotOpts.enablePosterTheme = false;
plotOpts.savePNG = true;
plotOpts.saveSVG = false;
plotOpts.themes = "normal";
plotOpts.keepFiguresOpen = true;

%% Generate plot
figDir = outDir;

fprintf('\nGenerating Tau vs k/d plot...\n');
plot_tau_vs_kdh_re(summaryRows, figDir, plotOpts);
fprintf('Done. Output in: %s\n', outDir);
