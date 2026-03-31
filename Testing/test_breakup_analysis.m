%% test_breakup_analysis.m
% Standalone test for the breakup analysis pipeline.
% Loads saved .mat data from a completed main run and re-runs:
%   - plot_breakup_gamma_vs_dratio
%   - write_breakup_analysis_xlsx
%
% Edit matDir below to point to the latest results run.

clear; clc;

%% Paths — edit matDir to point to your latest results run
matDir  = "E:\March Re 90,000 inception data\Processed images\results\results 30 local\plot_data_mat";
outDir  = fullfile(fileparts(mfilename('fullpath')), 'test_outputs', 'BreakupAnalysis');
xlsxOut = fullfile(fileparts(mfilename('fullpath')), 'test_outputs', 'breakup_events.xlsx');
if ~isfolder(outDir), mkdir(outDir); end

addpath(fileparts(fileparts(mfilename('fullpath'))));

%% Load allBreakup
S = load(fullfile(matDir, "breakup_analysis_by_case.mat"), 'allBreakup');
allBreakup = S.allBreakup;
nCases = numel(allBreakup);
fprintf('Loaded %d cases from breakup_analysis_by_case.mat\n\n', nCases);

%% Print per-case summary
fprintf('%-12s  %6s  %10s\n', 'Case', 'k/d', 'N events');
fprintf('%s\n', repmat('-', 1, 34));
for ci = 1:nCases
    nEv = numel(allBreakup(ci).events);
    fprintf('%-12s  %6.2f  %10d\n', ...
        allBreakup(ci).caseName, allBreakup(ci).kD, nEv);
end
fprintf('\n');

%% Plot options — edit here to test changes
set(0, 'DefaultFigureVisible', 'on');
plotOpts = struct();
plotOpts.enableNormalTheme  = true;
plotOpts.enablePosterTheme  = false;
plotOpts.savePNG            = true;
plotOpts.saveSVG            = false;
plotOpts.themes             = "normal";
plotOpts.keepFiguresOpen    = true;

%% Re-run plot
fprintf('Generating breakup gamma vs dRatio plot...\n');
plot_breakup_gamma_vs_dratio(allBreakup, outDir, plotOpts);

%% Write XLSX
fprintf('Writing XLSX...\n');
write_breakup_analysis_xlsx(allBreakup, xlsxOut);

fprintf('\nDone. Output in:\n  %s\n  %s\n', outDir, xlsxOut);
