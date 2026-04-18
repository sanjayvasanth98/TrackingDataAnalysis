function collapseData = analyze_collapse_events(out, pixelSize, dt, opts)
%ANALYZE_COLLAPSE_EVENTS  Detect bubble collapse events from TrackMate data.
%
% ===========================================================================
% COLLAPSE SIGNAL DEFINITION
% ===========================================================================
% Each tracked-bubble trajectory contributes at most ONE collapse event,
% defined as the FINAL FRAME of the track, subject to ALL of:
%
%   1. Basic validity: >= minTrackSpots spots with finite positions.
%      No flow-direction filter — all tracks in the field of view are
%      included (unlike the activation analysis which only uses left-moving
%      counterflow tracks).
%
%   2. Sufficient peak size: max(AREA over track) >= collapseMinPeakArea_px2.
%      Ensures the bubble genuinely grew to a detectable cavitation size,
%      filtering out micro-fragments and noise tracks.
%
%   3. Active shrinkage at track end: final AREA <=
%      collapseTruncationFactor × peakArea.
%      A bubble still near maximum size at its final frame most likely
%      left the camera field of view or was cut off by the recording
%      boundary — NOT a physical collapse. Default factor = 0.6, i.e. the
%      bubble must have shrunk to <= 60 % of peak before disappearing.
%
%   4. Not in unwanted ROI: the last tracked position must be outside the
%      masked-out region (e.g., throat wall area), consistent with the
%      rest of the analysis pipeline.
%
% PHYSICAL RATIONALE
%   In TrackMate tracking, a track ends when the bubble's image blob
%   becomes undetectable.  For criterion 3 to pass, the blob must have
%   actively shrunk before disappearing — consistent with physical
%   collapse back below the detection threshold, whereas a still-large
%   blob vanishing from one frame to the next points to the bubble
%   leaving the frame.
%
% OUTPUT SIGNAL
%   collapseCount(f) = number of qualified collapse events at frame f, on a
%   uniform integer grid [frameMin : frameMax] derived from out.spots.
%   This discrete time series is used for rate estimation and FFT analysis.
% ===========================================================================
%
% Inputs
%   out       struct from TrackMate parser: .trajectories (struct array),
%             .spots (table with at least ID, FRAME, AREA columns)
%   pixelSize mm per pixel (used for ROI mask lookup)
%   dt        seconds per frame
%   opts      struct (all fields optional, defaults applied):
%     .minTrackSpots            (default 3)
%     .collapseMinPeakArea_px2  (default 50)
%     .collapseTruncationFactor (default 0.6)
%     .fftNDomFreqs             (default 5) number of dominant peaks
%     .roiUnwantedMask          (default []) logical 2-D array
%     .roiPixelSize             (default 0)  mm/px of the mask image
%
% Output: collapseData struct — fields listed in the packing block below.

if nargin < 4, opts = struct(); end
opts = apply_collapse_defaults(opts);

traj  = out.trajectories;
spots = out.spots;

% ---- Build SpotID -> AREA lookup ----------------------------------------
hasArea = istable(spots) && ismember('AREA', spots.Properties.VariableNames);
areaById = containers.Map('KeyType','double','ValueType','double');
if hasArea
    for si = 1:height(spots)
        sid = double(spots.ID(si));
        if isfinite(sid)
            areaById(sid) = double(spots.AREA(si));
        end
    end
end

% ---- Full frame range from spots table ----------------------------------
frameMin = 0; frameMax = 0;
if istable(spots) && ismember('FRAME', spots.Properties.VariableNames)
    sFrames = double(spots.FRAME);
    sFrames = sFrames(isfinite(sFrames));
    if ~isempty(sFrames)
        frameMin = min(sFrames);
        frameMax = max(sFrames);
    end
end

% ---- Process each trajectory -------------------------------------------
nTraj      = numel(traj);
cfFrame    = nan(nTraj,1);
cfX        = nan(nTraj,1);
cfY        = nan(nTraj,1);
cfTrackId  = nan(nTraj,1);
cfPeakArea = nan(nTraj,1);
cfLastArea = nan(nTraj,1);

nQualified   = 0;
nTruncated   = 0;
nROIRejected = 0;

for k = 1:nTraj
    tk = traj(k);

    if ~isfield(tk,'x_phys') || ~isfield(tk,'y_phys') || ~isfield(tk,'frame')
        continue;
    end
    xv = double(tk.x_phys(:));
    yv = double(tk.y_phys(:));
    fv = double(tk.frame(:));
    n  = numel(xv);

    % Basic validity
    if n < opts.minTrackSpots, continue; end
    if numel(yv) ~= n || numel(fv) ~= n, continue; end
    if any(~isfinite(xv)) || any(~isfinite(yv)), continue; end

    % Area time series via spot IDs
    areaVals = nan(n,1);
    if hasArea && isfield(tk,'spotIds')
        sids = tk.spotIds(:);
        for ii = 1:min(numel(sids), n)
            sid = double(sids(ii));
            if isfinite(sid) && isKey(areaById, sid)
                areaVals(ii) = areaById(sid);
            end
        end
    end

    validA = areaVals(isfinite(areaVals) & areaVals > 0);
    if isempty(validA), continue; end

    peakA = max(validA);
    if peakA < opts.collapseMinPeakArea_px2, continue; end

    % Last valid area value
    lastA_all = areaVals(isfinite(areaVals) & areaVals >= 0);
    if isempty(lastA_all)
        lastA = validA(end);
    else
        lastA = lastA_all(end);
    end

    % Criterion 3: active shrinkage gate
    if lastA > opts.collapseTruncationFactor * peakA
        nTruncated = nTruncated + 1;
        continue;
    end

    % Criterion 4: ROI gate on final position
    xLast = xv(end);
    yLast = yv(end);
    if ~isempty(opts.roiUnwantedMask) && isfinite(opts.roiPixelSize) && opts.roiPixelSize > 0
        if pt_in_roi(xLast, yLast, opts.roiUnwantedMask, opts.roiPixelSize)
            nROIRejected = nROIRejected + 1;
            continue;
        end
    end

    % Qualified — record
    fLast = fv(end);
    if ~isfinite(fLast), fLast = frameMax; end

    nQualified = nQualified + 1;
    cfFrame(nQualified)    = fLast;
    cfX(nQualified)        = xLast;
    cfY(nQualified)        = yLast;
    if isfield(tk, 'TRACK_ID') && isfinite(tk.TRACK_ID)
        cfTrackId(nQualified) = double(tk.TRACK_ID);
    end
    cfPeakArea(nQualified) = peakA;
    cfLastArea(nQualified) = lastA;
end

% Trim pre-allocated arrays to actual count
cfFrame    = cfFrame(1:nQualified);
cfX        = cfX(1:nQualified);
cfY        = cfY(1:nQualified);
cfTrackId  = cfTrackId(1:nQualified);
cfPeakArea = cfPeakArea(1:nQualified);
cfLastArea = cfLastArea(1:nQualified);

% ---- Uniform frame grid -------------------------------------------------
if nQualified > 0
    frameMin = min(frameMin, min(cfFrame));
    frameMax = max(frameMax, max(cfFrame));
end
frameAxis     = (frameMin : frameMax).';
nFrames       = numel(frameAxis);
collapseCount = zeros(nFrames, 1);
for k = 1:nQualified
    idx = round(cfFrame(k)) - round(frameMin) + 1;
    if idx >= 1 && idx <= nFrames
        collapseCount(idx) = collapseCount(idx) + 1;
    end
end

% ---- Rate summaries -----------------------------------------------------
ratePerFrame = nQualified / max(nFrames, 1);
ratePerSec   = ratePerFrame / max(dt, eps);

% ---- FFT ----------------------------------------------------------------
[fftFreq_Hz, fftPower, domFreqs_Hz, domFreqPowers] = ...
    compute_collapse_fft(collapseCount, dt, opts.fftNDomFreqs);

fprintf('  Collapse: %d qualified | %d truncated | %d ROI-rejected | %.4g/frame | %.4g/s\n', ...
    nQualified, nTruncated, nROIRejected, ratePerFrame, ratePerSec);

% ---- Pack output --------------------------------------------------------
collapseData.collapseFrame   = cfFrame;
collapseData.collapseX_mm    = cfX;
collapseData.collapseY_mm    = cfY;
collapseData.collapseTrackId  = cfTrackId;
collapseData.peakArea_px2    = cfPeakArea;
collapseData.lastArea_px2    = cfLastArea;
collapseData.frameMin        = frameMin;
collapseData.frameMax        = frameMax;
collapseData.frameAxis       = frameAxis;
collapseData.collapseCount   = collapseCount;
collapseData.nTotal          = nTraj;
collapseData.nQualified      = nQualified;
collapseData.nTruncated      = nTruncated;
collapseData.nROIRejected    = nROIRejected;
collapseData.ratePerFrame    = ratePerFrame;
collapseData.ratePerSec      = ratePerSec;
collapseData.dt              = dt;
collapseData.fftFreq_Hz      = fftFreq_Hz;
collapseData.fftPower        = fftPower;
collapseData.domFreqs_Hz     = domFreqs_Hz;
collapseData.domFreqPowers   = domFreqPowers;
end


% =========================================================================
function [freqHz, power, domFreqs, domPowers] = compute_collapse_fft(signal, dt, nDom)
%COMPUTE_COLLAPSE_FFT  Single-sided FFT power spectrum of collapse-count signal.
%   Hann-windowed, zero-mean (DC removed).
signal  = double(signal(:));
N       = numel(signal);
freqHz  = nan(0,1); power = nan(0,1);
domFreqs = nan(0,1); domPowers = nan(0,1);

if N < 8 || ~isfinite(dt) || dt <= 0
    return;
end

% Hann window + DC removal
w   = 0.5 * (1 - cos(2*pi*(0:N-1).' / (N-1)));
sig = (signal - mean(signal)) .* w;

% FFT normalized by window RMS so amplitude is independent of window choice
wRMS  = sqrt(mean(w.^2));
F     = fft(sig) / (N * wRMS);
nPos  = floor(N/2) + 1;

% Single-sided amplitude (double non-DC, non-Nyquist bins)
A        = abs(F(1:nPos)) * 2;
A(1)     = abs(F(1));           % DC
if mod(N,2) == 0, A(end) = A(end)/2; end  % Nyquist (not doubled)

power  = A .^ 2;
freqHz = (0:nPos-1).' / (N * dt);

[domFreqs, domPowers] = find_dominant_freqs(freqHz, power, nDom);
end


% =========================================================================
function [freqVals, powerVals] = find_dominant_freqs(freqHz, power, nTop)
%FIND_DOMINANT_FREQS  Top-N local-maxima peaks in the power spectrum.
%   Skips DC (index 1). Uses 5-bin moving-average smoothing for peak
%   detection but reports unsmoothed power at the peak frequency.
power  = power(:);
freqHz = freqHz(:);
n      = numel(power);
freqVals  = nan(0,1);
powerVals = nan(0,1);
if n < 3, return; end

% Smooth for peak detection
smoothP = conv(power, ones(5,1)/5, 'same');

% Local maxima, skip DC (index 1)
isPeak = false(n,1);
for i = 2:(n-1)
    if smoothP(i) >= smoothP(i-1) && smoothP(i) >= smoothP(i+1) && power(i) > 0
        isPeak(i) = true;
    end
end
if n >= 2 && power(n) > power(n-1)
    isPeak(n) = true;
end

pidx = find(isPeak);
if isempty(pidx), return; end

[~, ord]  = sort(power(pidx), 'descend');
nActual   = min(nTop, numel(pidx));
topIdx    = pidx(ord(1:nActual));
freqVals  = freqHz(topIdx);
powerVals = power(topIdx);
end


% =========================================================================
function result = pt_in_roi(x_mm, y_mm, mask, pixelSize)
[nRows, nCols] = size(mask);
c = round(x_mm / pixelSize);
r = round(y_mm / pixelSize);
result = isfinite(c) && isfinite(r) && ...
    c >= 1 && c <= nCols && r >= 1 && r <= nRows && mask(r,c);
end


% =========================================================================
function opts = apply_collapse_defaults(opts)
if ~isfield(opts,'minTrackSpots'),            opts.minTrackSpots = 3;    end
if ~isfield(opts,'collapseMinPeakArea_px2'),  opts.collapseMinPeakArea_px2 = 50; end
if ~isfield(opts,'collapseTruncationFactor'), opts.collapseTruncationFactor = 0.6; end
if ~isfield(opts,'fftNDomFreqs'),             opts.fftNDomFreqs = 5;     end
if ~isfield(opts,'roiUnwantedMask'),          opts.roiUnwantedMask = []; end
if ~isfield(opts,'roiPixelSize'),             opts.roiPixelSize = 0;     end
end
