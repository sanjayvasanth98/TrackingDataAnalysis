%% test_plot_ai_vs_kd.m
% Standalone test for plot_ai_vs_kdh_re — A/I vs k/d plot.
% Loads saved .mat data and calls the plot function with the same
% arguments as main_batch_trackmate_local.m.
%
% Edit plotOpts below to experiment, then reflect changes back to main code.

clear; clc;

addpath(fileparts(mfilename('fullpath')));
addpath(fileparts(fileparts(mfilename('fullpath'))));

%% Paths
matDir = test_plotmat_location("activation_summary_by_case.mat");
outDir = fullfile(fileparts(mfilename('fullpath')), 'test_outputs');
if ~isfolder(outDir), mkdir(outDir); end

%% Load data
S = load(fullfile(matDir, "activation_summary_by_case.mat"), 'summaryRows');
summaryRows = S.summaryRows;
fprintf('Loaded %d cases from activation_summary_by_case.mat\n', height(summaryRows));
disp(summaryRows(:, {'Case','Re','kD','A_over_I','tau_mean','tau_sem'}));

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
fitTxtFile = fullfile(outDir, "fit_AI_vs_kD_by_Re.txt");

fprintf('\nGenerating A/I vs k/d plot...\n');
plot_ai_vs_kdh_re(summaryRows, figDir, fitTxtFile, plotOpts);
fprintf('Done. Output in: %s\n', outDir);
