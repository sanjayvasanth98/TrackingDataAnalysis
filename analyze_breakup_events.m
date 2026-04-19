function events = analyze_breakup_events(xmlFile, pixelSize, varargin)
% ANALYZE_BREAKUP_EVENTS  Extract elongated-parent breakup events from TrackMate XML.
%
%   Refined logic:
%     1. Parse ALL tracks (no direction filter).
%     2. Only consider tracks with NUMBER_SPLITS > 0.
%     3. Within each split track, find split nodes: spots that appear as
%        SPOT_SOURCE_ID in 2+ edges (one parent -> multiple children).
%     4. Parent filter: ELLIPSE_ASPECTRATIO > aspectRatioMin (default 4).
%     5. Child filter: AREA > childAreaMin_px2 (default 100 px^2).
%        Microbubble children (<= 100 px^2) are discarded.
%     6. Valid event: at least 2 children pass the area filter.
%        (Ensures genuine breakup, not just shedding one microbubble.)
%     7. Apply unwanted area ROI mask to both parent and each child.
%     8. For each valid child: record gamma, dRatio as one data point.
%
%   Returns struct array, one element per valid child-parent pair:
%     .gamma        (x_child - x_parent) / dRoughnessSpacing_mm
%     .dRatio       d_child / d_parent  (equivalent diameters)
%     .parentAR     parent ELLIPSE_ASPECTRATIO
%     .parentArea   parent AREA (px^2)
%     .childArea    child  AREA (px^2)
%     .x_parent_mm, .y_parent_mm  (mm, image coords)
%     .x_child_mm,  .y_child_mm
%     .parentFrame  frame number of the parent spot
%     .childFrame   frame number of the child spot
%     .parentID     TrackMate spot ID of parent
%     .childID      TrackMate spot ID of child
%     .nChildren    total number of valid (non-microbubble) children at this split
%     .trackID      TrackMate TRACK_ID containing this split

p = inputParser;
addParameter(p, 'roiData',                 []);
addParameter(p, 'aspectRatioMin',          4.0);
addParameter(p, 'childAreaMin_px2',        100.0);
addParameter(p, 'dRoughnessSpacing_mm',    0.384);
addParameter(p, 'maxTracks',               Inf);
parse(p, varargin{:});

roiData    = p.Results.roiData;
arMin      = p.Results.aspectRatioMin;
childMin   = p.Results.childAreaMin_px2;
dRough     = p.Results.dRoughnessSpacing_mm;
maxTracks  = p.Results.maxTracks;

% Pre-allocate output as empty struct array with correct fields.
emptyEv = struct('gamma',[], 'dRatio',[], 'parentAR',[], ...
    'parentArea',[], 'childArea',[], ...
    'x_parent_mm',[], 'y_parent_mm',[], ...
    'x_child_mm',[], 'y_child_mm',[], ...
    'parentFrame',[], 'childFrame',[], ...
    'parentID',[], 'childID',[], ...
    'nChildren',[], 'trackID',[]);
events = emptyEv([]);

fprintf('  Parsing XML: %s\n', xmlFile);
xDoc = xmlread(xmlFile);

% ---- Walk tracks with splits --------------------------------------------
allTracksEl = xDoc.getElementsByTagName('AllTracks');
if allTracksEl.getLength == 0
    warning('No AllTracks element found in XML: %s', xmlFile);
    return;
end
trackList = allTracksEl.item(0).getElementsByTagName('Track');
nTracks = trackList.getLength;
nTracksToInspect = nTracks;
if isfinite(maxTracks)
    nTracksToInspect = min(nTracks, max(0, floor(maxTracks)));
end

% ---- Build spot lookup (ID -> spot data) --------------------------------
% In smoke-test mode, restrict the map to spots referenced by the inspected
% tracks so the breakup pass follows the same small-track limit as the main
% TrackMate parser.
keepSpotIds = [];
restrictSpotMap = false;
if nTracksToInspect < nTracks
    keepSpotIds = collect_track_spot_ids(trackList, nTracksToInspect);
    restrictSpotMap = true;
end
fprintf('  Building spot map...\n');
spotMap = build_spot_map(xDoc, keepSpotIds, restrictSpotMap);
fprintf('  Loaded %d spots.\n', spotMap.Count);

nSplitTracks = 0;
nRawEvents   = 0;

for ti = 0:nTracksToInspect-1
    track = trackList.item(ti);
    nSplits = str2double(track.getAttribute('NUMBER_SPLITS'));
    if ~(isfinite(nSplits) && nSplits > 0), continue; end
    nSplitTracks = nSplitTracks + 1;
    trackID = str2double(track.getAttribute('TRACK_ID'));

    % Build source -> [target1, target2, ...] map for this track.
    edgeEls = track.getElementsByTagName('Edge');
    srcToTgts = containers.Map('KeyType', 'double', 'ValueType', 'any');
    for ei = 0:edgeEls.getLength-1
        e = edgeEls.item(ei);
        src = str2double(e.getAttribute('SPOT_SOURCE_ID'));
        tgt = str2double(e.getAttribute('SPOT_TARGET_ID'));
        if isnan(src) || isnan(tgt), continue; end
        if srcToTgts.isKey(src)
            srcToTgts(src) = [srcToTgts(src), tgt];
        else
            srcToTgts(src) = tgt;
        end
    end

    % Find split nodes (sources with 2+ targets).
    srcIds = keys(srcToTgts);
    for si = 1:numel(srcIds)
        src    = srcIds{si};
        targets = srcToTgts(src);
        if numel(targets) < 2, continue; end   % not a split node

        % ---- Parent spot checks ----------------------------------------
        if ~spotMap.isKey(src), continue; end
        pSpot = spotMap(src);
        if pSpot.aspectRatio < arMin, continue; end   % not elongated enough

        % Parent ROI check.
        xP = pSpot.x * pixelSize;
        yP = pSpot.y * pixelSize;
        if ~isempty(roiData) && pt_in_roi(xP, yP, roiData), continue; end

        dParent = sqrt(4 * pSpot.area / pi) * pixelSize;  % mm

        % ---- Child spot checks -----------------------------------------
        validTgts = [];
        for tgt = targets
            if ~spotMap.isKey(tgt), continue; end
            cSpot = spotMap(tgt);
            if cSpot.area >= childMin
                validTgts(end+1) = tgt; %#ok<AGROW>
            end
        end
        % Need >= 2 valid (non-microbubble) children for a genuine breakup.
        if numel(validTgts) < 2, continue; end

        nRawEvents = nRawEvents + 1;
        nValidChildren = numel(validTgts);

        % ---- Record one point per valid child ---------------------------
        for tgt = validTgts
            cSpot = spotMap(tgt);
            xC = cSpot.x * pixelSize;
            yC = cSpot.y * pixelSize;
            if ~isempty(roiData) && pt_in_roi(xC, yC, roiData), continue; end

            dChild = sqrt(4 * cSpot.area / pi) * pixelSize;  % mm

            ev = emptyEv;
            ev.gamma        = (xC - xP) / dRough;
            ev.dRatio       = dChild / dParent;
            ev.parentAR     = pSpot.aspectRatio;
            ev.parentArea   = pSpot.area;
            ev.childArea    = cSpot.area;
            ev.x_parent_mm  = xP;
            ev.y_parent_mm  = yP;
            ev.x_child_mm   = xC;
            ev.y_child_mm   = yC;
            ev.parentFrame  = pSpot.frame;
            ev.childFrame   = cSpot.frame;
            ev.parentID     = src;
            ev.childID      = tgt;
            ev.nChildren    = nValidChildren;
            ev.trackID      = trackID;
            events(end+1) = ev; %#ok<AGROW>
        end
    end
end

fprintf('  Split tracks: %d | Valid breakup events: %d | Child data points: %d\n', ...
    nSplitTracks, nRawEvents, numel(events));
end


% =========================================================================
function spotIds = collect_track_spot_ids(trackList, nTracksToInspect)
spotIds = nan(0,1);
for ti = 0:nTracksToInspect-1
    track = trackList.item(ti);
    edgeEls = track.getElementsByTagName('Edge');
    for ei = 0:edgeEls.getLength-1
        e = edgeEls.item(ei);
        src = str2double(e.getAttribute('SPOT_SOURCE_ID'));
        tgt = str2double(e.getAttribute('SPOT_TARGET_ID'));
        if isfinite(src), spotIds(end+1,1) = src; end %#ok<AGROW>
        if isfinite(tgt), spotIds(end+1,1) = tgt; end %#ok<AGROW>
    end
end
spotIds = unique(spotIds(isfinite(spotIds)));
end


% =========================================================================
function spotMap = build_spot_map(xDoc, keepSpotIds, restrictSpotMap)
% Returns containers.Map: spot_ID (double) -> struct(x,y,area,aspectRatio,...)
spotMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
if nargin < 2
    keepSpotIds = [];
end
if nargin < 3
    restrictSpotMap = false;
end
keepAll = ~restrictSpotMap;
if ~keepAll
    keepMap = containers.Map('KeyType', 'double', 'ValueType', 'logical');
    for ii = 1:numel(keepSpotIds)
        keepMap(double(keepSpotIds(ii))) = true;
    end
else
    keepMap = containers.Map('KeyType', 'double', 'ValueType', 'logical');
end

allSpotsEl = xDoc.getElementsByTagName('AllSpots');
if allSpotsEl.getLength == 0, return; end

frames = allSpotsEl.item(0).getElementsByTagName('SpotsInFrame');
for fi = 0:frames.getLength-1
    spotEls = frames.item(fi).getElementsByTagName('Spot');
    for si = 0:spotEls.getLength-1
        s  = spotEls.item(si);
        id = str2double(s.getAttribute('ID'));
        if isnan(id), continue; end
        if ~keepAll && ~isKey(keepMap, id), continue; end
        sp.x           = str2double(s.getAttribute('POSITION_X'));
        sp.y           = str2double(s.getAttribute('POSITION_Y'));
        sp.frame       = str2double(s.getAttribute('FRAME'));
        sp.area        = str2double(s.getAttribute('AREA'));
        sp.aspectRatio = str2double(s.getAttribute('ELLIPSE_ASPECTRATIO'));
        sp.majorAxis   = str2double(s.getAttribute('ELLIPSE_MAJOR'));
        sp.minorAxis   = str2double(s.getAttribute('ELLIPSE_MINOR'));
        spotMap(id) = sp;
    end
end
end


% =========================================================================
function result = pt_in_roi(x_mm, y_mm_image, roiData)
% Returns true if the point falls inside the unwanted track mask.
mask = roiData.unwantedTrackMask;
ps   = roiData.maskPixelSize;
[nRows, nCols] = size(mask);
c = round(x_mm      / ps);
r = round(y_mm_image / ps);
result = isfinite(c) && isfinite(r) && ...
    c >= 1 && c <= nCols && r >= 1 && r <= nRows && mask(r, c);
end
