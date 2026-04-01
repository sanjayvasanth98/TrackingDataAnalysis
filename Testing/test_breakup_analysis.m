%% test_breakup_analysis.m
% Standalone test for the breakup analysis pipeline.
% Loads saved .mat data from a completed main run and re-runs:
%   - plot_breakup_gamma_vs_dratio  (per AR threshold: 1.5, 2.0, 4.0)
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
plotOpts.savePNG            = true;
plotOpts.saveSVG            = false;
plotOpts.themes             = "normal";
plotOpts.keepFiguresOpen    = true;

%% AR thresholds — same as main batch
arThresholds = [1.5, 2.0, 3.0, 4.0];

%% Re-run plots for each AR threshold
for arThr = arThresholds
    arTag = sprintf('AR%s', strrep(sprintf('%.1f', arThr), '.', 'p'));
    filteredBreakup = filter_breakup_by_ar(allBreakup, arThr);

    nTotal = 0;
    for ci = 1:numel(filteredBreakup)
        nTotal = nTotal + numel(filteredBreakup(ci).events);
    end
    fprintf('Generating breakup plot for %s (%d events total)...\n', arTag, nTotal);

    plot_breakup_gamma_vs_dratio(filteredBreakup, outDir, plotOpts, arTag);
    save_square_copy(gcf, outDir, "Breakup_gamma_vs_dRatio_" + arTag + "_normal", plotOpts);

    plot_breakup_gamma_beeswarm_vs_kd(filteredBreakup, outDir, plotOpts, arTag, outDir);
    save_square_copy(gcf, outDir, "Breakup_gamma_beeswarm_kD_" + arTag + "_normal", plotOpts);

    plot_breakup_gamma_scatter_vs_ar(filteredBreakup, outDir, plotOpts, arTag, outDir);
    save_square_copy(gcf, outDir, "Breakup_gamma_scatter_AR_" + arTag + "_normal", plotOpts);
end

%% Write XLSX (full dataset)
fprintf('Writing XLSX...\n');
write_breakup_analysis_xlsx(allBreakup, xlsxOut);

fprintf('\nDone. Output in:\n  %s\n  %s\n', outDir, xlsxOut);


%% ========================================================================
function save_square_copy(fig, outDir, baseName, plotOpts)
%SAVE_SQUARE_COPY  Resize an open figure to square and save a second copy.
%   The rectangle version is already saved by the plot function.
%   This creates an additional _square version.
    origPos = fig.Position;
    fig.Position = [100 100 700 700];
    drawnow;
    save_fig_dual_safe(fig, fullfile(outDir, baseName + "_square"), plotOpts);
    fig.Position = origPos;   % restore rectangle for display
end
