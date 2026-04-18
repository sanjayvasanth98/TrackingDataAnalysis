%% test_plot_inception_locations.m
% Standalone test for plot_inception_locations_by_re — inception location
% scatter + marginal histograms + 2D KDE contour overlays.
% Loads saved .mat data and calls the plot function with the same
% arguments as main_batch_trackmate_local.m.
%
% Edit plotOpts below to experiment, then reflect changes back to main code.

clear; clc;

%% Paths
matDir = "C:\Users\kbsanjayvasanth\Downloads\plot_data_mat";
outDir = fullfile(fileparts(mfilename('fullpath')), 'test_outputs', 'InceptionLocations');
if ~isfolder(outDir), mkdir(outDir); end

addpath(fileparts(fileparts(mfilename('fullpath'))));

%% Load data
S = load(fullfile(matDir, "inception_locations_by_case.mat"), 'allLoc');
allLoc = S.allLoc;
nCases = numel(allLoc.caseName);
fprintf('Loaded %d cases from inception_locations_by_case.mat\n', nCases);
for i = 1:nCases
    nPts = 0;
    if ~isempty(allLoc.inception2x_xy{i})
        nPts = size(allLoc.inception2x_xy{i}, 1);
    end
    fprintf('  %s: Re=%g, k/d=%.4g, pixelSize=%.6g, %d activation points\n', ...
        allLoc.caseName(i), allLoc.Re(i), allLoc.kD(i), allLoc.pixelSize(i), nPts);
end

%% Load ROI data (wall mask, unwanted track mask, throat location)
roiFile = "E:\March Re 90,000 inception data\Processed images\results\2000 frames\Smooth variation 2000\ROI_throat.mat";
R = load(roiFile);
fprintf('Loaded ROI data from: %s\n', roiFile);
fprintf('  Throat at px: (%.1f, %.1f)\n', R.x_throat, R.y_throat);
fprintf('  Wall mask: %d pixels\n', sum(R.wallMask(:)));
fprintf('  Unwanted track mask: %d pixels\n', sum(R.unwantedTrackMask(:)));

roiData = struct();
roiData.wallMask = R.wallMask;
roiData.unwantedTrackMask = R.unwantedTrackMask;
roiData.throat_xy_px = R.ROI_throat.throat_xy_px;
roiData.maskPixelSize = median(allLoc.pixelSize);  % same calibration

%% Plot options (mirror main code — edit here to test changes)
set(0, 'DefaultFigureVisible', 'on');
plotOpts = struct();
plotOpts.enableNormalTheme = true;
plotOpts.enablePosterTheme = false;
plotOpts.savePNG = true;
plotOpts.saveSVG = false;
plotOpts.themes = "normal";
plotOpts.keepFiguresOpen = true;
plotOpts.inceptionImageSize_px = [1280 320];
plotOpts.inceptionXLim_mm = [0 4.8];
plotOpts.inceptionYLim_mm = [0 1.2];
plotOpts.roiData = roiData;

%% Generate plot
fprintf('\nGenerating inception location plots...\n');
plot_inception_locations_by_re(allLoc, outDir, plotOpts);
fprintf('Done. Output in: %s\n', outDir);
