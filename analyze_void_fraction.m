function vfData = analyze_void_fraction(out, opts)
%ANALYZE_VOID_FRACTION  Compute per-frame 2-D void fraction from TrackMate data.
%
% ===========================================================================
% VOID FRACTION DEFINITION
% ===========================================================================
% The 2-D projected void fraction at frame f is:
%
%   alpha(f) = SUM_i[ AREA_i(f) ] / fluidArea_px2
%
% where:
%   - The sum is over ALL detected spots at frame f in out.spots (no
%     flow-direction or activation filter — every detected bubble blob
%     contributes its pixel area).
%   - fluidArea_px2 is the fluid-accessible area of the FOV in pixels:
%
%       fluidArea_px2 = fovArea_px2 - maskedArea_px2
%
%     fovArea_px2   = (ceil(max X) + 1) x (ceil(max Y) + 1), derived from
%                     spot position columns (can be overridden via opts).
%     maskedArea_px2 = number of true pixels in roiUnwantedMask, scaled by
%                     (roiPixelSize / cameraPixelSize)^2 when the mask was
%                     created at a different resolution.
%
%   Subtracting the wall / unwanted-track mask area gives the physically
%   correct denominator: the fraction of the FLOW DOMAIN occupied by gas,
%   not the fraction of all pixels including opaque walls.
%
% The frame grid is the uniform integer range [frameMin : frameMax] taken
% from out.spots.FRAME.  Frames with no detected spots contribute alpha = 0.
%
% Aggregate metrics:
%   vfMean   = mean(alpha)  over all frames in [frameMin, frameMax]
%   vfStd    = std(alpha)   over the same range
%   vfMedian = median(alpha)
%
% WHY ALL SPOTS (NOT JUST TRACKED)?
%   The spots table reflects every blob the detector found in each frame,
%   including bubbles too short-lived to form a track.  Using tracks alone
%   would underestimate the true area fraction when bubbles appear and
%   disappear within a single frame or are rejected by the tracker.
%
% Inputs
%   out   struct from TrackMate parser, requires:
%           .spots  table with columns: ID, FRAME, AREA, X, Y
%   opts  struct (all fields optional):
%     .fovWidth_px      override FOV width  in pixels (default: auto from X)
%     .fovHeight_px     override FOV height in pixels (default: auto from Y)
%     .roiUnwantedMask  logical 2-D array — pixels to exclude from denominator
%                       (same mask used by collapse / activation analyses)
%     .roiPixelSize     mm/px of the mask image (default: same as camera)
%     .cameraPixelSize  mm/px of camera, used to scale mask area
%                       (default: same as roiPixelSize → scale factor 1)
%
% Output: vfData struct with fields:
%   fovWidth_px, fovHeight_px, fovArea_px2
%   maskedArea_px2       pixels excluded from denominator (wall + unwanted)
%   fluidArea_px2        effective denominator = fovArea - maskedArea
%   frameAxis            (nFrames x 1) uniform integer frame indices
%   vfPerFrame           (nFrames x 1) void fraction per frame [0,1]
%   vfMean               scalar mean void fraction [0,1]
%   vfStd                scalar std  void fraction [0,1]
%   vfMedian             scalar median             [0,1]
%   nFrames              total frames in grid
%   nFramesWithSpots     frames where >= 1 spot was detected
% ===========================================================================

if nargin < 2, opts = struct(); end

spots = out.spots;

% ---- Validate spots table -----------------------------------------------
requiredCols = {'FRAME','AREA','X','Y'};
if ~istable(spots)
    warning('analyze_void_fraction: out.spots is not a table. Returning empty.');
    vfData = empty_vf_data();
    return;
end
missingCols = requiredCols(~ismember(requiredCols, spots.Properties.VariableNames));
if ~isempty(missingCols)
    warning('analyze_void_fraction: out.spots missing columns: %s. Returning empty.', ...
        strjoin(missingCols, ', '));
    vfData = empty_vf_data();
    return;
end

sFrame = double(spots.FRAME);
sArea  = double(spots.AREA);
sX     = double(spots.X);
sY     = double(spots.Y);

% Remove rows with non-finite values
validRows = isfinite(sFrame) & isfinite(sArea) & sArea >= 0 & isfinite(sX) & isfinite(sY);
sFrame = sFrame(validRows);
sArea  = sArea(validRows);
sX     = sX(validRows);
sY     = sY(validRows);

if isempty(sFrame)
    warning('analyze_void_fraction: no valid spots rows. Returning empty.');
    vfData = empty_vf_data();
    return;
end

% ---- FOV bounding box (in camera pixels) --------------------------------
if isfield(opts,'fovWidth_px') && isfinite(opts.fovWidth_px) && opts.fovWidth_px > 0
    fovW = opts.fovWidth_px;
else
    fovW = ceil(max(sX)) + 1;
end
if isfield(opts,'fovHeight_px') && isfinite(opts.fovHeight_px) && opts.fovHeight_px > 0
    fovH = opts.fovHeight_px;
else
    fovH = ceil(max(sY)) + 1;
end
fovArea = fovW * fovH;

% ---- Masked (wall + unwanted-track) area --------------------------------
maskedArea = 0;
if isfield(opts,'roiUnwantedMask') && ~isempty(opts.roiUnwantedMask) && islogical(opts.roiUnwantedMask)
    maskPx = sum(opts.roiUnwantedMask(:));

    % Scale mask pixels to camera pixels if mask was made at a different resolution
    scaleFactor = 1;
    if isfield(opts,'roiPixelSize') && isfield(opts,'cameraPixelSize') && ...
            isfinite(opts.roiPixelSize) && opts.roiPixelSize > 0 && ...
            isfinite(opts.cameraPixelSize) && opts.cameraPixelSize > 0
        scaleFactor = (opts.roiPixelSize / opts.cameraPixelSize)^2;
    end
    maskedArea = maskPx * scaleFactor;
end

fluidArea = fovArea - maskedArea;
if fluidArea <= 0
    warning('analyze_void_fraction: fluidArea <= 0 after mask subtraction. Using full FOV.');
    fluidArea = fovArea;
    maskedArea = 0;
end

% ---- Frame grid ---------------------------------------------------------
frameMin  = min(sFrame);
frameMax  = max(sFrame);
frameAxis = (frameMin : frameMax).';
nFrames   = numel(frameAxis);

% ---- Per-frame area sum (accumarray) ------------------------------------
frameIdx = round(sFrame - frameMin) + 1;
frameIdx = max(1, min(frameIdx, nFrames));  % safety clamp

areaSum = accumarray(frameIdx, sArea, [nFrames, 1], @sum, 0);

% ---- Void fraction ------------------------------------------------------
vfPerFrame = areaSum / fluidArea;
vfPerFrame = min(max(vfPerFrame, 0), 1);  % clamp numerical overshoot

vfMean   = mean(vfPerFrame);
vfStd    = std(vfPerFrame);
vfMedian = median(vfPerFrame);

nFramesWithSpots = sum(areaSum > 0);

% ---- Pack output --------------------------------------------------------
vfData = struct();
vfData.fovWidth_px       = fovW;
vfData.fovHeight_px      = fovH;
vfData.fovArea_px2       = fovArea;
vfData.maskedArea_px2    = maskedArea;
vfData.fluidArea_px2     = fluidArea;
vfData.frameAxis         = frameAxis;
vfData.vfPerFrame        = vfPerFrame;
vfData.vfMean            = vfMean;
vfData.vfStd             = vfStd;
vfData.vfMedian          = vfMedian;
vfData.nFrames           = nFrames;
vfData.nFramesWithSpots  = nFramesWithSpots;

fprintf('  Void fraction: mean=%.4g%% std=%.4g%% median=%.4g%% | fluid=%d px2 (FOV=%dx%d, masked=%d px) | frames=%d (%d with spots)\n', ...
    vfMean*100, vfStd*100, vfMedian*100, round(fluidArea), fovW, fovH, round(maskedArea), nFrames, nFramesWithSpots);
end


% =========================================================================
function vfData = empty_vf_data()
vfData = struct('fovWidth_px',NaN,'fovHeight_px',NaN,'fovArea_px2',NaN, ...
    'maskedArea_px2',NaN,'fluidArea_px2',NaN, ...
    'frameAxis',[],'vfPerFrame',[],'vfMean',NaN,'vfStd',NaN,'vfMedian',NaN, ...
    'nFrames',0,'nFramesWithSpots',0);
end
