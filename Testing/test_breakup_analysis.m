%% test_breakup_analysis.m
% Standalone test for the breakup analysis pipeline.
% Loads saved .mat data from a completed main run and re-runs:
%   - breakup plots (saved as square-format outputs only)
%   - write_breakup_analysis_xlsx
%
% Edit matDir below to point to the latest results run.

clear; clc;

%% Paths — edit matDir to point to your latest results run
matDir  = "C:\Users\kbsanjayvasanth\Downloads\plot_data_mat";
outDir  = fullfile(fileparts(mfilename('fullpath')), 'test_outputs', 'BreakupAnalysis');
xlsxOut = fullfile(fileparts(mfilename('fullpath')), 'test_outputs', 'breakup_events.xlsx');
if ~isfolder(outDir), mkdir(outDir); end

addpath(fileparts(fileparts(mfilename('fullpath'))));

%% Load allBreakup
S = load(fullfile(matDir, "breakup_analysis_by_case.mat"), 'allBreakup');
allBreakup = S.allBreakup;
nCases = numel(allBreakup);
fprintf('Loaded %d cases from breakup_analysis_by_case.mat\n\n', nCases);

%% Print per-case summary (all events, lowest AR)
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
plotOpts.savePNG            = false; % suppress default rectangular saves
plotOpts.saveSVG            = false; % square copies are saved manually below
plotOpts.themes             = "normal";
plotOpts.keepFiguresOpen    = true;
plotOpts.breakupDRatioXScale = "log";
plotOpts.breakupDRatioXLim = [0.11, 1.4];
plotOpts.breakupDRatioClipLowPercentile = 0.5;
plotOpts.breakupDRatioClipPercentile = 99.5;
plotOpts.breakupGammaYLim = [];
plotOpts.breakupGammaYClipPercentile = [1, 99];
plotOpts.breakupDRatioMarkerSize = 30;
plotOpts.breakupDRatioMarkerAlpha = 0.40;
plotOpts.breakupDRatioTrendMaxBins = 12;
plotOpts.breakupDRatioTrendMinCount = 5;

%% Global gamma-vs-dRatio plot (no AR filtering)
nTotalAll = 0;
for ci = 1:nCases
    nTotalAll = nTotalAll + numel(allBreakup(ci).events);
end
fprintf('Generating breakup plot for all events (%d events total)...\n', nTotalAll);
plot_breakup_gamma_vs_dratio(allBreakup, outDir, plotOpts);
save_square_only(gcf, outDir, "Breakup_gamma_vs_dRatio_normal", plotOpts);

fprintf('Generating gamma vs parent AR plot for all events (%d events total)...\n', nTotalAll);
plot_breakup_gamma_scatter_vs_ar(allBreakup, outDir, plotOpts, "", outDir);
save_square_only(gcf, outDir, "Breakup_gamma_scatter_AR_normal", plotOpts);

fprintf('Generating gamma vs k/d beeswarm for all events (%d events total)...\n', nTotalAll);
plot_breakup_gamma_beeswarm_vs_kd(allBreakup, outDir, plotOpts, "", outDir);
save_square_only(gcf, outDir, "Breakup_gamma_beeswarm_kD_normal", plotOpts);

%% AR thresholds — same as main batch for the other breakup plots
arThresholds = [1.5, 2.0, 3.0, 4.0];

%% Re-run the AR-filtered breakup plots
for arThr = arThresholds
    arTag = sprintf('AR%s', strrep(sprintf('%.1f', arThr), '.', 'p'));
    filteredBreakup = filter_breakup_by_ar(allBreakup, arThr);

    nTotal = 0;
    for ci = 1:numel(filteredBreakup)
        nTotal = nTotal + numel(filteredBreakup(ci).events);
    end
    fprintf('Generating AR-filtered breakup plots for %s (%d events total)...\n', arTag, nTotal);

    plot_breakup_gamma_beeswarm_vs_kd(filteredBreakup, outDir, plotOpts, arTag, outDir);
    save_square_only(gcf, outDir, "Breakup_gamma_beeswarm_kD_" + arTag + "_normal", plotOpts);
end

%% Write XLSX (full dataset)
fprintf('Writing XLSX...\n');
write_breakup_analysis_xlsx(allBreakup, xlsxOut);

fprintf('\nDone. Output in:\n  %s\n  %s\n', outDir, xlsxOut);


%% ========================================================================
function save_square_only(fig, outDir, baseName, plotOpts)
%SAVE_SQUARE_ONLY  Resize an open figure to square and save only that version.
    origPos = fig.Position;
    fig.Position = [100 100 700 700];
    drawnow;
    squarePlotOpts = plotOpts;
    if ~isfield(squarePlotOpts, 'savePNG') || ~squarePlotOpts.savePNG
        squarePlotOpts.savePNG = true;
    end
    save_fig_dual_safe(fig, fullfile(outDir, baseName + "_square"), squarePlotOpts);
    fig.Position = origPos;
end
