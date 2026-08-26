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
outDir = fullfile(fileparts(mfilename('fullpath')), 'test_outputs', 'AI_vs_kD');
if ~isfolder(outDir), mkdir(outDir); end

%% Load data
S = load(fullfile(matDir, "activation_summary_by_case.mat"), 'summaryRows');
summaryRows = S.summaryRows;
fprintf('Loaded %d cases from activation_summary_by_case.mat\n', height(summaryRows));
displayCols = {'Case','Re','kD','A_over_I','A_over_I_absolute','tau_mean','tau_sem'};
displayCols = displayCols(ismember(displayCols, summaryRows.Properties.VariableNames));
disp(summaryRows(:, displayCols));

%% Plot options (mirror main code — edit here to test changes)
set(0, 'DefaultFigureVisible', 'on');
plotOpts = struct();
plotOpts.enableNormalTheme = true;
plotOpts.enablePosterTheme = false;
plotOpts.savePNG = true;
plotOpts.saveSVG = false;
plotOpts.themes = "normal";
plotOpts.keepFiguresOpen = true;
plotOpts.aiVsKDYScale = "linear";

%% Generate plot
figDir = outDir;

fprintf('\nGenerating A/I vs k/d plot...\n');
plotOpts.aiMetricField = "A_over_I";
plotOpts.aiMetricOutputStem = "A_over_I_vs_kD_by_Re";
fitTxtFile = fullfile(outDir, "fit_A_over_I_vs_kD_by_Re.txt");
plot_ai_vs_kdh_re(summaryRows, figDir, fitTxtFile, plotOpts);

if ismember('A_over_I_absolute', summaryRows.Properties.VariableNames)
    fprintf('\nGenerating A/I absolute vs k/d plot...\n');
    plotOpts.aiMetricField = "A_over_I_absolute";
    plotOpts.aiMetricOutputStem = "A_over_I_absolute_vs_kD_by_Re";
    fitTxtFile = fullfile(outDir, "fit_A_over_I_absolute_vs_kD_by_Re.txt");
    plot_ai_vs_kdh_re(summaryRows, figDir, fitTxtFile, plotOpts);
else
    warning('summaryRows does not contain A_over_I_absolute. Skipping absolute A/I plot.');
end
fprintf('Done. Output in: %s\n', outDir);
