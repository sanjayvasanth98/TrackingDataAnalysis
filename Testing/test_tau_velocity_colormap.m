%% test_tau_velocity_colormap.m
% Standalone test for residence-time plots coloured by activated velocity.
% Loads saved .mat data from a completed main run and re-runs:
%   - plot_tau_velocity_colormap
%
% Edit matDir below to point to either:
%   combined: resultsDir/plot_data_mat
%   chunk:    resultsDir/results individual/chunk_N/plot_data_mat

clear; clc;

%% Paths - edit matDir to point to your latest results run
matDir = "C:\Users\kbsanjayvasanth\Downloads\plot_data_mat";
outDir = fullfile(fileparts(mfilename('fullpath')), 'test_outputs', 'TauVelocityColormap');
if ~isfolder(outDir), mkdir(outDir); end

addpath(fileparts(fileparts(mfilename('fullpath'))));

%% Load activation summary data
summaryFile = fullfile(matDir, "activation_summary_by_case.mat");
S = load(summaryFile);
if isfield(S, 'summaryRows')
    summaryTable = S.summaryRows;
elseif isfield(S, 'chunkSR')
    summaryTable = S.chunkSR;
else
    error('activation_summary_by_case.mat must contain summaryRows or chunkSR.');
end

required = {'Case','Re','kD','tau_mean','tau_sem','activatedVelocity_mean_m_s'};
missing = required(~ismember(required, summaryTable.Properties.VariableNames));
if ~isempty(missing)
    error(['Missing columns needed for velocity-coloured tau plot: %s\n', ...
        'Re-run main_batch_trackmate_local.m or main_batch_trackmate_arc.m so the new velocity metrics are saved.'], ...
        strjoin(missing, ', '));
end

fprintf('Loaded %d cases from:\n  %s\n\n', height(summaryTable), summaryFile);

%% Print per-case summary
fprintf('%-12s  %8s  %8s  %12s  %12s\n', ...
    'Case', 'Re', 'k/d', 'tau_mean', 'v_mean_m_s');
fprintf('%s\n', repmat('-', 1, 62));
for ci = 1:height(summaryTable)
    fprintf('%-12s  %8g  %8.4g  %12.4g  %12.4g\n', ...
        char(string(summaryTable.Case(ci))), ...
        summaryTable.Re(ci), ...
        summaryTable.kD(ci), ...
        summaryTable.tau_mean(ci), ...
        summaryTable.activatedVelocity_mean_m_s(ci));
end
fprintf('\n');

%% Plot options - edit here to test changes
set(0, 'DefaultFigureVisible', 'on');
plotOpts = struct();
plotOpts.enableNormalTheme = true;
plotOpts.enablePosterTheme = false;
plotOpts.savePNG = true;
plotOpts.saveSVG = false;
plotOpts.themes = "normal";
plotOpts.keepFiguresOpen = true;

%% Re-run plot
fprintf('Generating velocity-coloured residence-time plot...\n');
plot_tau_velocity_colormap(summaryTable, outDir, plotOpts, outDir);

fprintf('\nDone. Output in:\n  %s\n', outDir);
