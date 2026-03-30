%% test_breakup_analysis.m
% Runs breakup event detection across all cases and produces the
% gamma vs d_child/d_parent scatter + mean-line plot.
%
% Adjust the parameters block if needed before running.

clear; clc;

%% ---- Cases ---------------------------------------------------------------
cases(1).name      = "5um";
cases(1).kD        = 0.030;
cases(1).xmlFile   = "E:\March Re 90,000 inception data\Processed images\Smooth variation 2\test_nofilter_smoothvar2_48lit.xml";
cases(1).pixelSize = 0.00375009375;

cases(2).name      = "12um";
cases(2).kD        = 0.080;
cases(2).xmlFile   = "E:\March Re 90,000 inception data\Processed images\P10S100\P10S100_48lit.xml";
cases(2).pixelSize = 0.00375009375;

cases(3).name      = "20um";
cases(3).kD        = 0.141;
cases(3).xmlFile   = "E:\March Re 90,000 inception data\Processed images\P10S70\P10S70_48lit.xml";
cases(3).pixelSize = 0.00375009375;

cases(4).name      = "30um";
cases(4).kD        = 0.267;
cases(4).xmlFile   = "E:\March Re 90,000 inception data\Processed images\P10S50\P10S50_48lit.xml";
cases(4).pixelSize = 0.00375009375;

cases(5).name      = "53um";
cases(5).kD        = 0.444;
cases(5).xmlFile   = "E:\March Re 90,000 inception data\Processed images\P10S30\P10S30_48lit.xml";
cases(5).pixelSize = 0.00375009375;

cases(6).name      = "80um";
cases(6).kD        = 0.720;
cases(6).xmlFile   = "E:\March Re 90,000 inception data\Processed images\P10S20\P10S20_48lit.xml";
cases(6).pixelSize = 0.00375009375;

%% ---- Parameters ----------------------------------------------------------
aspectRatioMin       = 4.0;     % parent ELLIPSE_ASPECTRATIO threshold
childAreaMin_px2     = 100.0;   % minimum child area (px^2) — excludes microbubbles
dRoughnessSpacing_mm = 0.384;   % roughness spacing d (mm)

%% ---- ROI mask (unwanted area) -------------------------------------------
roiFile = "E:\March Re 90,000 inception data\Processed images\results\2000 frames\Smooth variation 2000\ROI_throat.mat";
roiData = [];
if isfile(roiFile)
    R = load(roiFile);
    roiData = struct();
    roiData.unwantedTrackMask = R.ROI_throat.unwantedTrackMask;
    roiData.maskPixelSize     = cases(1).pixelSize;
    fprintf('Loaded ROI mask: %d unwanted pixels.\n', sum(roiData.unwantedTrackMask(:)));
else
    warning('ROI file not found — mask filtering disabled.');
end

%% ---- Output directory ---------------------------------------------------
outDir = "E:\March Re 90,000 inception data\Processed images\results\2000 frames\Smooth variation 2000\Figures_PNG_SVG\Breakup";
if ~isfolder(outDir), mkdir(outDir); end

%% ---- Run detection per case ---------------------------------------------
nCases = numel(cases);
breakupData(nCases) = struct('caseName', "", 'kD', 0, 'events', []);

for ci = 1:nCases
    fprintf('\n[%d/%d] Case: %s (k/d = %.4g)\n', ci, nCases, cases(ci).name, cases(ci).kD);
    breakupData(ci).caseName = cases(ci).name;
    breakupData(ci).kD       = cases(ci).kD;
    breakupData(ci).events   = analyze_breakup_events( ...
        cases(ci).xmlFile, cases(ci).pixelSize, ...
        'roiData',                roiData, ...
        'aspectRatioMin',         aspectRatioMin, ...
        'childAreaMin_px2',       childAreaMin_px2, ...
        'dRoughnessSpacing_mm',   dRoughnessSpacing_mm);
end

%% ---- Summary table -------------------------------------------------------
fprintf('\n%-10s  %6s  %10s\n', 'Case', 'k/d', 'N points');
fprintf('%s\n', repmat('-', 1, 32));
for ci = 1:nCases
    fprintf('%-10s  %6.4g  %10d\n', ...
        breakupData(ci).caseName, breakupData(ci).kD, numel(breakupData(ci).events));
end

%% ---- Validation XLSX (one sheet per case) --------------------------------
xlsxFile = fullfile(outDir, "breakup_events_validation.xlsx");
headers = {'Case','kD','TrackID','ParentSpotID','ParentFrame', ...
    'ParentArea_px2','ParentAR','nChildren', ...
    'ChildSpotID','ChildFrame','ChildArea_px2','gamma','dRatio'};

% Delete existing file so old sheets don't persist.
if isfile(xlsxFile), delete(xlsxFile); end

for ci = 1:nCases
    ev = breakupData(ci).events;
    sheetName = char(breakupData(ci).caseName);  % e.g. "12um"

    if isempty(ev)
        % Write header-only sheet so every case is represented.
        writecell(headers, xlsxFile, 'Sheet', sheetName, 'Range', 'A1');
        continue;
    end

    nEv = numel(ev);
    rows = cell(nEv, numel(headers));
    for ei = 1:nEv
        rows(ei,:) = { ...
            char(breakupData(ci).caseName), breakupData(ci).kD, ...
            ev(ei).trackID,    ev(ei).parentID,   ev(ei).parentFrame, ...
            ev(ei).parentArea, ev(ei).parentAR,   ev(ei).nChildren, ...
            ev(ei).childID,    ev(ei).childFrame,  ev(ei).childArea, ...
            ev(ei).gamma,      ev(ei).dRatio};
    end

    T = cell2table(rows, 'VariableNames', headers);
    writetable(T, xlsxFile, 'Sheet', sheetName);
end

fprintf('Validation XLSX saved: %s\n', xlsxFile);

%% ---- Plot ----------------------------------------------------------------
plotOpts.themes          = ["normal", "poster"];
plotOpts.keepFiguresOpen = false;
plotOpts.savePng         = true;
plotOpts.saveSvg         = true;

plot_breakup_gamma_vs_dratio(breakupData, outDir, plotOpts);

fprintf('\nDone. Figures saved to:\n  %s\n', outDir);
