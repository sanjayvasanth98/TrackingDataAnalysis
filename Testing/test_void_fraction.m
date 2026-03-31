%% test_void_fraction.m
% Standalone test for the void fraction analysis pipeline.
% Loads saved .mat data from a completed main run and re-runs:
%   - plot_void_fraction_vs_kd
%
% Edit matDir and plotOpts below to experiment, then reflect any changes
% back to main_batch_trackmate_local.m.

clear; clc;

%% Paths — edit matDir to point to your latest results run
matDir = "E:\March Re 90,000 inception data\Processed images\results\results 30 local\plot_data_mat";
outDir = fullfile(fileparts(mfilename('fullpath')), 'test_outputs', 'VoidFraction');
if ~isfolder(outDir), mkdir(outDir); end

addpath(fileparts(fileparts(mfilename('fullpath'))));

%% Load allVoidFrac
S = load(fullfile(matDir, "void_fraction_by_case.mat"), 'allVoidFrac');
allVoidFrac = S.allVoidFrac;
nCases = numel(allVoidFrac.caseName);
fprintf('Loaded %d cases from void_fraction_by_case.mat\n\n', nCases);

%% Print per-case summary
fprintf('%-12s  %6s  %8s  %10s  %8s  %8s  %12s  %12s  %12s\n', ...
    'Case','Re','k/d','FOV(px)','masked%','nFrames','alpha_mean%','alpha_std%','alpha_med%');
fprintf('%s\n', repmat('-',1,100));
for ci = 1:nCases
    d = allVoidFrac.data{ci};
    if isempty(d) || ~isfinite(d.vfMean)
        fprintf('%-12s  (no data)\n', allVoidFrac.caseName{ci});
        continue;
    end
    maskedPct = d.maskedArea_px2 / d.fovArea_px2 * 100;
    fprintf('%-12s  %6g  %8.4g  %4dx%-4d  %7.2f%%  %8d  %12.4g  %12.4g  %12.4g\n', ...
        allVoidFrac.caseName{ci}, ...
        allVoidFrac.Re(ci), ...
        allVoidFrac.kD(ci), ...
        d.fovWidth_px, d.fovHeight_px, ...
        maskedPct, ...
        d.nFrames, ...
        d.vfMean   * 100, ...
        d.vfStd    * 100, ...
        d.vfMedian * 100);
end
fprintf('\n');

%% Plot options
set(0, 'DefaultFigureVisible', 'on');
plotOpts = struct();
plotOpts.enableNormalTheme  = true;
plotOpts.enablePosterTheme  = false;
plotOpts.savePNG            = true;
plotOpts.saveSVG            = false;
plotOpts.themes             = "normal";
plotOpts.keepFiguresOpen    = true;

%% Re-run plot
fprintf('Generating void fraction vs k/d plot...\n');
plot_void_fraction_vs_kd(allVoidFrac, outDir, plotOpts);

fprintf('\nDone. Output in: %s\n', outDir);
