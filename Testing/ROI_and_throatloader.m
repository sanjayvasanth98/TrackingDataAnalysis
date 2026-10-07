%% ------------------------------------------------------------
% Interactive Throat / Wall / Unwanted-Area Selection Tool
% Saves all selected data into ROI_throat.mat
% ------------------------------------------------------------
clear; clc; close all;
fprintf('--- Throat / Wall / Unwanted-Area Selection Tool ---\n');

%% --- Step 1: Select file (video OR image) ---
[fileName, filePath] = uigetfile({ ...
    '*.avi;*.mp4;*.mov;*.mkv;*.png;*.jpg;*.jpeg;*.tif;*.bmp', ...
    'Video or Image Files (*.avi, *.mp4, *.mov, *.mkv, *.png, *.jpg, *.jpeg, *.tif, *.bmp)'}, ...
    'Select a Video or Image File');

if isequal(fileName, 0)
    error('No file selected. Exiting.');
end

fullPath = fullfile(filePath, fileName);
[~, baseName, ext] = fileparts(fileName);
fprintf('Selected file: %s\n', fullPath);

%% --- Step 2: Create save folder ---
saveFolder = fullfile(filePath, baseName);
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
    fprintf('Created folder: %s\n', saveFolder);
else
    fprintf('Folder already exists: %s\n', saveFolder);
end

saveMatFile = fullfile(saveFolder, 'ROI_throat.mat');
summaryPlotFile = fullfile(saveFolder, [baseName '_selection_summary.png']);

%% --- Step 3: Load first frame or image ---
isVideo = ismember(lower(ext), {'.avi','.mp4','.mov','.mkv'});

if isVideo
    fprintf('File detected as VIDEO.\n');
    v = VideoReader(fullPath);
    rawFrame = readFrame(v);
else
    fprintf('File detected as IMAGE.\n');
    rawFrame = imread(fullPath);
end

if ndims(rawFrame) == 3
    frame1 = rgb2gray(rawFrame);
else
    frame1 = rawFrame;
end

frameSize = size(frame1);

%% --- Step 4: Select throat location ---
fprintf('\nStep 1/3: Select throat location.\n');
hFig1 = figure('Name', 'Step 1/3: Throat Location', 'NumberTitle', 'off', 'Color', 'w');
imshow(frame1, []);
axis on;
xlabel('x (px)');
ylabel('y (px)');
title({'Step 1/3: Click throat location', 'Click once, then press Enter'});

[x_throat, y_throat] = ginput(1);
if isempty(x_throat) || isempty(y_throat)
    error('No throat location selected. Exiting.');
end

hold on;
plot(x_throat, y_throat, 'r+', 'MarkerSize', 12, 'LineWidth', 2);
text(x_throat + 10, y_throat, 'Throat', 'Color', 'r', 'FontWeight', 'bold');
fprintf('Throat selected at (x, y) = (%.3f, %.3f) px\n', x_throat, y_throat);

%% --- Step 5: Select wall region ---
fprintf('\nStep 2/3: Draw wall region polygon.\n');
fprintf('Double-click to close polygon. Press Enter without drawing to skip.\n');

hFig2 = figure('Name', 'Step 2/3: Wall Region', 'NumberTitle', 'off', 'Color', 'w');
imshow(frame1, []);
axis on;
xlabel('x (px)');
ylabel('y (px)');
title({'Step 2/3: Draw wall region polygon', 'Double-click to close. Press Enter to skip'});

wallMask = roipoly;
if isempty(wallMask)
    wallMask = false(frameSize);
    fprintf('No wall region selected.\n');
else
    fprintf('Wall region selected.\n');
end

%% --- Step 6: Select unwanted-track area ---
fprintf('\nStep 3/3: Draw unwanted-track area polygon.\n');
fprintf('Double-click to close polygon. Press Enter without drawing to skip.\n');

hFig3 = figure('Name', 'Step 3/3: Unwanted Track Area', 'NumberTitle', 'off', 'Color', 'w');
imshow(frame1, []);
axis on;
xlabel('x (px)');
ylabel('y (px)');
title({'Step 3/3: Draw unwanted-track area polygon', 'Double-click to close. Press Enter to skip'});

unwantedTrackMask = roipoly;
if isempty(unwantedTrackMask)
    unwantedTrackMask = false(frameSize);
    fprintf('No unwanted-track area selected.\n');
else
    fprintf('Unwanted-track area selected.\n');
end

%% --- Step 7: Save everything into ROI_throat.mat ---
ROI_throat = struct();
ROI_throat.sourceFile = fullPath;
ROI_throat.fileName = fileName;
ROI_throat.baseName = baseName;
ROI_throat.isVideo = isVideo;
ROI_throat.frameSize = frameSize;
ROI_throat.x_throat = x_throat;
ROI_throat.y_throat = y_throat;
ROI_throat.throat_xy_px = [x_throat, y_throat];
ROI_throat.wallMask = wallMask;
ROI_throat.unwantedTrackMask = unwantedTrackMask;

save(saveMatFile, ...
    'ROI_throat', ...
    'x_throat', 'y_throat', ...
    'wallMask', 'unwantedTrackMask', ...
    'frame1', 'frameSize', ...
    'fullPath', 'fileName', 'baseName', 'isVideo');

fprintf('\nSaved ROI/throat data to: %s\n', saveMatFile);

%% --- Step 8: Show and save summary overlay ---
hFig4 = figure('Name', 'Selection Summary', 'NumberTitle', 'off', 'Color', 'w');
imshow(frame1, []);
axis on;
xlabel('x (px)');
ylabel('y (px)');
hold on;

hThroat = plot(x_throat, y_throat, 'r+', 'MarkerSize', 12, 'LineWidth', 2);

legendHandles = hThroat;
legendLabels = {'Throat'};

if any(wallMask(:))
    draw_cross_hatched_mask(wallMask, [1 1 0], 18, 0.8);
    contour(double(wallMask), [0.5 0.5], 'Color', [1 1 0], 'LineWidth', 1.5);
    hWall = line(nan, nan, 'Color', [1 1 0], 'LineWidth', 1.5);
    legendHandles(end+1) = hWall; %#ok<AGROW>
    legendLabels{end+1} = 'Wall region'; %#ok<AGROW>
end

if any(unwantedTrackMask(:))
    draw_cross_hatched_mask(unwantedTrackMask, [0 1 1], 18, 0.8);
    contour(double(unwantedTrackMask), [0.5 0.5], 'Color', [0 1 1], 'LineWidth', 1.5);
    hUnwanted = line(nan, nan, 'Color', [0 1 1], 'LineWidth', 1.5);
    legendHandles(end+1) = hUnwanted; %#ok<AGROW>
    legendLabels{end+1} = 'Unwanted track area'; %#ok<AGROW>
end

lgd = legend(legendHandles, legendLabels, 'Location', 'northeast');
lgd.Box = 'on';

exportgraphics(hFig4, summaryPlotFile, 'Resolution', 300);
fprintf('Saved selection summary plot: %s\n', summaryPlotFile);

fprintf('\n--- Selection Complete ---\n');
fprintf('Throat location saved.\n');
fprintf('Wall region saved.\n');
fprintf('Unwanted-track area saved.\n');
fprintf('Output MAT file: %s\n', saveMatFile);
fprintf('Summary plot saved: %s\n', summaryPlotFile);

function draw_cross_hatched_mask(mask, colorRGB, spacingPx, lineWidth)
[nRows, nCols] = size(mask);

plot_hatch_segments(mask, colorRGB, spacingPx, lineWidth, +1, nRows, nCols);
plot_hatch_segments(mask, colorRGB, spacingPx, lineWidth, -1, nRows, nCols);
end

function plot_hatch_segments(mask, colorRGB, spacingPx, lineWidth, slopeSign, nRows, nCols)
if slopeSign > 0
    offsets = -nRows:spacingPx:nCols;
    for offset = offsets
        rowIdx = 1:nRows;
        colIdx = rowIdx + offset;

        valid = colIdx >= 1 & colIdx <= nCols;
        rowIdx = rowIdx(valid);
        colIdx = colIdx(valid);

        if isempty(rowIdx)
            continue;
        end

        inside = mask(sub2ind([nRows, nCols], rowIdx, colIdx));
        draw_inside_segments(colIdx, rowIdx, inside, colorRGB, lineWidth);
    end
else
    sums = 1:spacingPx:(nRows + nCols);
    for s = sums
        rowIdx = 1:nRows;
        colIdx = s - rowIdx;

        valid = colIdx >= 1 & colIdx <= nCols;
        rowIdx = rowIdx(valid);
        colIdx = colIdx(valid);

        if isempty(rowIdx)
            continue;
        end

        inside = mask(sub2ind([nRows, nCols], rowIdx, colIdx));
        draw_inside_segments(colIdx, rowIdx, inside, colorRGB, lineWidth);
    end
end
end

function draw_inside_segments(xVals, yVals, insideMask, colorRGB, lineWidth)
insideMask = logical(insideMask(:)');
xVals = xVals(:)';
yVals = yVals(:)';

transitions = diff([false, insideMask, false]);
starts = find(transitions == 1);
stops = find(transitions == -1) - 1;

for k = 1:numel(starts)
    idx = starts(k):stops(k);
    line(xVals(idx), yVals(idx), 'Color', colorRGB, 'LineWidth', lineWidth);
end
end
