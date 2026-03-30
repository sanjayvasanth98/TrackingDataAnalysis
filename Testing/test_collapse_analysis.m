%% test_collapse_analysis.m
% Standalone test for the collapse-frequency analysis pipeline.
% Loads saved .mat data from a completed main run and re-runs:
%   - plot_collapse_rate_vs_frame
%   - plot_collapse_power_spectrum
%   - write_collapse_analysis_csv
%
% Edit plotOpts and collapseOpts below to experiment, then reflect
% any changes back to main_batch_trackmate_local.m.

clear; clc;

%% Paths — edit matDir to point to your latest results run
matDir   = "E:\March Re 90,000 inception data\Processed images\results\results 27\plot_data_mat";
outDir   = fullfile(fileparts(mfilename('fullpath')), 'test_outputs', 'CollapseAnalysis');
csvOut   = fullfile(fileparts(mfilename('fullpath')), 'test_outputs', 'collapse_analysis.csv');
if ~isfolder(outDir), mkdir(outDir); end

addpath(fileparts(fileparts(mfilename('fullpath'))));

%% Load allCollapse
S = load(fullfile(matDir, "collapse_analysis_by_case.mat"), 'allCollapse');
allCollapse = S.allCollapse;
nCases = numel(allCollapse.caseName);
fprintf('Loaded %d cases from collapse_analysis_by_case.mat\n\n', nCases);

%% Print per-case summary
fprintf('%-12s  %6s  %8s  %10s  %10s  %12s  %12s\n', ...
    'Case','Re','k/d','nTotal','nQualified','rate/frame','rate/s');
fprintf('%s\n', repmat('-',1,80));
for ci = 1:nCases
    cd = allCollapse.data{ci};
    if isempty(cd)
        fprintf('%-12s  (no data)\n', allCollapse.caseName{ci});
        continue;
    end
    fprintf('%-12s  %6g  %8.4g  %10d  %10d  %12.4g  %12.4g\n', ...
        allCollapse.caseName{ci}, ...
        allCollapse.Re(ci), ...
        allCollapse.kD(ci), ...
        cd.nTotal, ...
        cd.nQualified, ...
        cd.ratePerFrame, ...
        cd.ratePerSec);
    if ~isempty(cd.domFreqs_Hz)
        fprintf('  Dominant frequencies (Hz): ');
        fprintf('%.1f  ', cd.domFreqs_Hz);
        fprintf('\n');
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
fprintf('Generating collapse rate vs frame plot...\n');
plot_collapse_rate_vs_frame(allCollapse, outDir, plotOpts);

fprintf('Generating collapse power spectrum plot...\n');
plot_collapse_power_spectrum(allCollapse, outDir, plotOpts);

fprintf('Writing CSV...\n');
write_collapse_analysis_csv(allCollapse, csvOut);

fprintf('\nDone. Output in:\n  %s\n  %s\n', outDir, csvOut);
