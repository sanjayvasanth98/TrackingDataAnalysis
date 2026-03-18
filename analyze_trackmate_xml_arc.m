function out = analyze_trackmate_xml_arc(xmlFile, opts)
%ANALYZE_TRACKMATE_XML_ARC  Parse TrackMate XML (v7+) and build MATLAB tables/trajectories.
% Cluster-safe version: NO plots unless opts.makePlots=true.

arguments
    xmlFile (1,:) char
    opts.pixelSize (1,1) double = 1.0
    opts.dt        (1,1) double = 1.0
    opts.maxTracks (1,1) double = Inf
    opts.parseTrackedSpotsOnly (1,1) logical = true
    opts.parseFilteredTracksOnly (1,1) logical = true
    opts.parseROI  (1,1) logical = false
    opts.verbose   (1,1) logical = true
    opts.makePlots (1,1) logical = false
end

xmlFile = sanitize_xml_path(xmlFile);

if ~isfile(xmlFile)
    error("XML file not found: %s", xmlFile);
end

if opts.verbose
    fprintf("Reading XML: %s\n", xmlFile);
end

doc = read_trackmate_xml_document(xmlFile);

%% -------------------------
%  Parse Tracks + Edges
% --------------------------
trackNodes = doc.getElementsByTagName('Track');
nTracksAll = double(trackNodes.getLength);
[filteredTrackIds, hasFilteredTracks] = parse_filtered_track_ids(doc);
useFilteredTrackSet = opts.parseFilteredTracksOnly && hasFilteredTracks;

if opts.parseFilteredTracksOnly && ~hasFilteredTracks
    warning('FilteredTracks block not found in %s. Falling back to all tracks.', xmlFile);
end
if useFilteredTrackSet && isempty(filteredTrackIds)
    warning('FilteredTracks block in %s is empty. Falling back to all tracks.', xmlFile);
    useFilteredTrackSet = false;
end

maxTracks = max(0, floor(opts.maxTracks));
if opts.verbose
    if isfinite(maxTracks)
        maxTrackMsg = sprintf('%d', maxTracks);
    else
        maxTrackMsg = 'all';
    end
    if useFilteredTrackSet
        fprintf("Found %d Track nodes; FilteredTracks has %d IDs; parsing up to %s filtered tracks.\n", ...
            nTracksAll, numel(filteredTrackIds), maxTrackMsg);
    else
        fprintf("Found %d Track nodes; parsing up to %s tracks.\n", nTracksAll, maxTrackMsg);
    end
end

TRACK_ID     = nan(0,1);
TRACK_INDEX  = nan(0,1);
N_SPOTS      = nan(0,1);
N_GAPS       = nan(0,1);
N_SPLITS     = nan(0,1);
N_MERGES     = nan(0,1);
N_COMPLEX    = nan(0,1);
TRACK_START  = nan(0,1);
TRACK_STOP   = nan(0,1);
DURATION     = nan(0,1);
DISP         = nan(0,1);
MEAN_SPEED   = nan(0,1);
MEAN_QUAL    = nan(0,1);
TOTAL_DIST   = nan(0,1);

edgeTrackId = nan(0,1);
edgeSource  = nan(0,1);
edgeTarget  = nan(0,1);
edgeCost    = nan(0,1);
edgeSpeed   = nan(0,1);
edgeDisp    = nan(0,1);
edgeTime    = nan(0,1);

filteredTrackMap = [];
if useFilteredTrackSet
    filteredTrackMap = containers.Map('KeyType', 'double', 'ValueType', 'logical');
    for ii = 1:numel(filteredTrackIds)
        filteredTrackMap(filteredTrackIds(ii)) = true;
    end
end

for ti = 1:nTracksAll
    if isfinite(maxTracks) && numel(TRACK_ID) >= maxTracks
        break;
    end

    tnode = trackNodes.item(ti-1);
    trackId = getAttrNum(tnode,'TRACK_ID');

    if useFilteredTrackSet
        if ~isKey(filteredTrackMap, trackId)
            continue;
        end
    end

    TRACK_ID(end+1,1)    = trackId; %#ok<AGROW>
    TRACK_INDEX(end+1,1) = getAttrNum(tnode,'TRACK_INDEX'); %#ok<AGROW>
    N_SPOTS(end+1,1)     = getAttrNum(tnode,'NUMBER_SPOTS'); %#ok<AGROW>
    N_GAPS(end+1,1)      = getAttrNum(tnode,'NUMBER_GAPS'); %#ok<AGROW>
    N_SPLITS(end+1,1)    = getAttrNum(tnode,'NUMBER_SPLITS'); %#ok<AGROW>
    N_MERGES(end+1,1)    = getAttrNum(tnode,'NUMBER_MERGES'); %#ok<AGROW>
    N_COMPLEX(end+1,1)   = getAttrNum(tnode,'NUMBER_COMPLEX'); %#ok<AGROW>
    TRACK_START(end+1,1) = getAttrNum(tnode,'TRACK_START'); %#ok<AGROW>
    TRACK_STOP(end+1,1)  = getAttrNum(tnode,'TRACK_STOP'); %#ok<AGROW>
    DURATION(end+1,1)    = getAttrNum(tnode,'TRACK_DURATION'); %#ok<AGROW>
    DISP(end+1,1)        = getAttrNum(tnode,'TRACK_DISPLACEMENT'); %#ok<AGROW>
    MEAN_SPEED(end+1,1)  = getAttrNum(tnode,'TRACK_MEAN_SPEED'); %#ok<AGROW>
    MEAN_QUAL(end+1,1)   = getAttrNum(tnode,'TRACK_MEAN_QUALITY'); %#ok<AGROW>
    TOTAL_DIST(end+1,1)  = getAttrNum(tnode,'TOTAL_DISTANCE_TRAVELED'); %#ok<AGROW>

    eNodes = tnode.getElementsByTagName('Edge');
    ne = eNodes.getLength;
    if ne == 0, continue; end

    edgeTrackId = [edgeTrackId; repmat(trackId, ne, 1)]; %#ok<AGROW>

    tmpSrc = nan(ne,1); tmpTgt = nan(ne,1);
    tmpCost = nan(ne,1); tmpSp = nan(ne,1); tmpDisp = nan(ne,1); tmpTime = nan(ne,1);

    for ei = 1:ne
        enode = eNodes.item(ei-1);
        tmpSrc(ei)  = getAttrNum(enode,'SPOT_SOURCE_ID');
        tmpTgt(ei)  = getAttrNum(enode,'SPOT_TARGET_ID');
        tmpCost(ei) = getAttrNum(enode,'LINK_COST');
        tmpSp(ei)   = getAttrNum(enode,'SPEED');
        tmpDisp(ei) = getAttrNum(enode,'DISPLACEMENT');
        tmpTime(ei) = getAttrNum(enode,'EDGE_TIME');
    end

    edgeSource = [edgeSource; tmpSrc]; %#ok<AGROW>
    edgeTarget = [edgeTarget; tmpTgt]; %#ok<AGROW>
    edgeCost   = [edgeCost;   tmpCost]; %#ok<AGROW>
    edgeSpeed  = [edgeSpeed;  tmpSp]; %#ok<AGROW>
    edgeDisp   = [edgeDisp;   tmpDisp]; %#ok<AGROW>
    edgeTime   = [edgeTime;   tmpTime]; %#ok<AGROW>
end

tracks = table(TRACK_ID, TRACK_INDEX, N_SPOTS, N_GAPS, N_SPLITS, N_MERGES, N_COMPLEX, ...
    TRACK_START, TRACK_STOP, DURATION, DISP, MEAN_SPEED, MEAN_QUAL, TOTAL_DIST, ...
    'VariableNames', {'TRACK_ID','TRACK_INDEX','NUMBER_SPOTS','NUMBER_GAPS','NUMBER_SPLITS', ...
    'NUMBER_MERGES','NUMBER_COMPLEX','TRACK_START','TRACK_STOP','TRACK_DURATION', ...
    'TRACK_DISPLACEMENT','TRACK_MEAN_SPEED','TRACK_MEAN_QUALITY','TOTAL_DISTANCE_TRAVELED'});

edges = table(edgeTrackId, edgeSource, edgeTarget, edgeCost, edgeSpeed, edgeDisp, edgeTime, ...
    'VariableNames', {'TRACK_ID','SPOT_SOURCE_ID','SPOT_TARGET_ID','LINK_COST','SPEED','DISPLACEMENT','EDGE_TIME'});

if opts.verbose
    fprintf("Parsed %d tracks and extracted %d edges.\n", height(tracks), height(edges));
end

%% -------------------------
%  Parse Spots (optionally only those used by parsed tracks)
% --------------------------
spotNodes = doc.getElementsByTagName('Spot');
nSpotsAll = spotNodes.getLength;

requiredSpotIds = nan(0,1);
if ~isempty(edges) && height(edges) > 0
    requiredSpotIds = unique([edges.SPOT_SOURCE_ID; edges.SPOT_TARGET_ID]);
end

keepTrackedOnly = opts.parseTrackedSpotsOnly && ~isempty(requiredSpotIds);

if opts.verbose
    if keepTrackedOnly
        fprintf("Found %d Spot nodes; parsing tracked subset (~%d IDs).\n", nSpotsAll, numel(requiredSpotIds));
    else
        fprintf("Found %d Spot nodes; parsing all.\n", nSpotsAll);
    end
end

nAlloc = nSpotsAll;
if keepTrackedOnly
    nAlloc = numel(requiredSpotIds);
end

ID        = nan(nAlloc,1);
FRAME     = nan(nAlloc,1);
POS_T     = nan(nAlloc,1);
X         = nan(nAlloc,1);
Y         = nan(nAlloc,1);
Z         = nan(nAlloc,1);
QUALITY   = nan(nAlloc,1);
RADIUS    = nan(nAlloc,1);
AREA      = nan(nAlloc,1);
PERIM     = nan(nAlloc,1);
CIRC      = nan(nAlloc,1);
SOLIDITY  = nan(nAlloc,1);
MEAN_I1   = nan(nAlloc,1);
TOTAL_I1  = nan(nAlloc,1);
ROI_NPTS  = nan(nAlloc,1);
ROI_XY    = cell(nAlloc,1);

if keepTrackedOnly
    requiredMap = containers.Map('KeyType','double','ValueType','logical');
    for rr = 1:numel(requiredSpotIds)
        requiredMap(requiredSpotIds(rr)) = true;
    end
else
    requiredMap = [];
end

nKeep = 0;

for i = 1:nSpotsAll
    node = spotNodes.item(i-1);
    sid = getAttrNum(node, 'ID');

    if keepTrackedOnly
        if ~isKey(requiredMap, sid)
            continue;
        end
    end

    nKeep = nKeep + 1;

    ID(nKeep)       = sid;
    FRAME(nKeep)    = getAttrNum(node, 'FRAME');
    POS_T(nKeep)    = getAttrNum(node, 'POSITION_T');
    X(nKeep)        = getAttrNum(node, 'POSITION_X');
    Y(nKeep)        = getAttrNum(node, 'POSITION_Y');
    Z(nKeep)        = getAttrNum(node, 'POSITION_Z');
    QUALITY(nKeep)  = getAttrNum(node, 'QUALITY');
    RADIUS(nKeep)   = getAttrNum(node, 'RADIUS');
    AREA(nKeep)     = getAttrNum(node, 'AREA');
    PERIM(nKeep)    = getAttrNum(node, 'PERIMETER');
    CIRC(nKeep)     = getAttrNum(node, 'CIRCULARITY');
    SOLIDITY(nKeep) = getAttrNum(node, 'SOLIDITY');
    MEAN_I1(nKeep)  = getAttrNum(node, 'MEAN_INTENSITY_CH1');
    TOTAL_I1(nKeep) = getAttrNum(node, 'TOTAL_INTENSITY_CH1');

    ROI_NPTS(nKeep) = getAttrNum(node, 'ROI_N_POINTS');

    if opts.parseROI
        txt = char(node.getTextContent);
        nums = str2double(strsplit(strtrim(txt)));
        if numel(nums) >= 2 && mod(numel(nums),2)==0
            ROI_XY{nKeep} = reshape(nums, 2, []).';
        else
            ROI_XY{nKeep} = [];
        end
    end

    if keepTrackedOnly
        remove(requiredMap, sid);
        if requiredMap.Count == 0
            break;
        end
    end
end

ID        = ID(1:nKeep);
FRAME     = FRAME(1:nKeep);
POS_T     = POS_T(1:nKeep);
X         = X(1:nKeep);
Y         = Y(1:nKeep);
Z         = Z(1:nKeep);
QUALITY   = QUALITY(1:nKeep);
RADIUS    = RADIUS(1:nKeep);
AREA      = AREA(1:nKeep);
PERIM     = PERIM(1:nKeep);
CIRC      = CIRC(1:nKeep);
SOLIDITY  = SOLIDITY(1:nKeep);
MEAN_I1   = MEAN_I1(1:nKeep);
TOTAL_I1  = TOTAL_I1(1:nKeep);
ROI_NPTS  = ROI_NPTS(1:nKeep);
ROI_XY    = ROI_XY(1:nKeep);

spots = table(ID, FRAME, POS_T, X, Y, Z, QUALITY, RADIUS, AREA, PERIM, CIRC, SOLIDITY, MEAN_I1, TOTAL_I1, ROI_NPTS);
if opts.parseROI
    spots.ROI_XY = ROI_XY;
end

if isempty(edges) || height(edges) == 0
    out = struct();
    out.spots = spots;
    out.tracks = tracks;
    out.edges = edges;
    out.trajectories = struct('TRACK_ID',{},'spotIds',{},'frame',{},'t',{},'x',{},'y',{}, ...
        'x_phys',{},'y_phys',{},'speed_phys',{},'ds_phys',{},'dt',{});
    out.meta = build_parser_meta(opts, useFilteredTrackSet, hasFilteredTracks, filteredTrackIds, nTracksAll, tracks, edges, spots);
    return;
end

%% -------------------------
%  Reconstruct ordered trajectories (per track)
% --------------------------
trajTemplate = struct('TRACK_ID', NaN, 'spotIds', [], 'frame', [], 't', [], 'x', [], 'y', [], ...
    'x_phys', [], 'y_phys', [], 'speed_phys', [], 'ds_phys', [], 'dt', []);

uniqueTracks = unique(edges.TRACK_ID);
traj = repmat(trajTemplate, numel(uniqueTracks), 1);

% Pre-sort spot IDs once for fast lookup
spotIDs = spots.ID;
[spotIDsSorted, idxSorted] = sort(spotIDs);

for k = 1:numel(uniqueTracks)
    tid = uniqueTracks(k);
    sub = edges(edges.TRACK_ID == tid, :);

    src = sub.SPOT_SOURCE_ID;
    tgt = sub.SPOT_TARGET_ID;

    startCandidates = setdiff(src, tgt);
    if isempty(startCandidates)
        start = src(1);
    else
        start = startCandidates(1);
    end

    % Build mapping source->target (first occurrence)
    % (TrackMate tracks typically linear; if branches exist we take first)
    [uSrc, ia] = unique(src, 'stable');
    uTgt = tgt(ia);

    % Walk the chain
    ordered = start;
    while true
        j = find(uSrc == ordered(end), 1, 'first');
        if isempty(j), break; end
        nxt = uTgt(j);
        if any(ordered == nxt), break; end
        ordered(end+1) = nxt; %#ok<AGROW>
    end

    n = numel(ordered);

    % Map ordered spot IDs to spots table rows
    % ismember on sorted list, then convert to original row via idxSorted
    [tf, loc] = ismember(ordered, spotIDsSorted);
    rowIdx = nan(n,1);
    rowIdx(tf) = idxSorted(loc(tf));

    fr = nan(n,1); tt = nan(n,1); xx = nan(n,1); yy = nan(n,1);
    good = isfinite(rowIdx);

    fr(good) = spots.FRAME(rowIdx(good));
    tt(good) = spots.POS_T(rowIdx(good));
    xx(good) = spots.X(rowIdx(good));
    yy(good) = spots.Y(rowIdx(good));

    % Sort by frame to be safe
    [fr, ordIdx] = sort(fr);
    ordered = ordered(ordIdx);
    tt = tt(ordIdx); xx = xx(ordIdx); yy = yy(ordIdx);

    xPhys = xx * opts.pixelSize;
    yPhys = yy * opts.pixelSize;

    % Time axis (fully toolbox/version safe)
    ttv = tt(isfinite(tt));
    if isempty(ttv) || (max(ttv) - min(ttv)) == 0
        tt = fr * opts.dt;
    else
        dtt = diff(tt);
        dtt = dtt(isfinite(dtt));
        if ~isempty(dtt) && abs(median(dtt) - 1) < 1e-6
            tt = tt * opts.dt;
        end
    end

    dX = diff(xPhys); dY = diff(yPhys);
    ds = hypot(dX, dY);
    dT = diff(tt);
    sp = ds ./ dT;

    traj(k).TRACK_ID = tid;
    traj(k).spotIds  = ordered(:);
    traj(k).frame    = fr(:);
    traj(k).t        = tt(:);
    traj(k).x        = xx(:);
    traj(k).y        = yy(:);
    traj(k).x_phys   = xPhys(:);
    traj(k).y_phys   = yPhys(:);
    traj(k).ds_phys  = ds(:);
    traj(k).dt       = dT(:);
    traj(k).speed_phys = sp(:);
end

%% Optional plots (OFF by default)
if opts.makePlots
    trackLengths = arrayfun(@(s) numel(s.spotIds), traj);
    figure('Name','Track length distribution');
    histogram(trackLengths);
    xlabel('Number of spots per track'); ylabel('Count'); grid on;

    figure('Name','Trajectories overlay'); hold on; axis ij; axis equal;
    title('Example trajectories (first 50 tracks)'); xlabel('x (px)'); ylabel('y (px)');
    nShow = min(50, numel(traj));
    for i = 1:nShow
        plot(traj(i).x, traj(i).y, '-');
    end
    hold off; grid on;

    allSpeeds = vertcat(traj.speed_phys);
    allSpeeds = allSpeeds(isfinite(allSpeeds));
    figure('Name','Speed histogram');
    histogram(allSpeeds, 50);
    xlabel('Speed (unit/s)'); ylabel('Count'); grid on;
end

%% Pack outputs
out = struct();
out.spots = spots;
out.tracks = tracks;
out.edges = edges;
out.trajectories = traj;
out.meta = build_parser_meta(opts, useFilteredTrackSet, hasFilteredTracks, filteredTrackIds, nTracksAll, tracks, edges, spots);

end

%% -------------------------
% Helpers
% --------------------------
function meta = build_parser_meta(opts, useFilteredTrackSet, hasFilteredTracks, filteredTrackIds, nTracksAll, tracks, edges, spots)
meta = struct();
meta.parserVersion = 2;
meta.parseTrackedSpotsOnly = logical(opts.parseTrackedSpotsOnly);
meta.parseFilteredTracksOnly = logical(useFilteredTrackSet);
meta.hasFilteredTracks = logical(hasFilteredTracks);
meta.nFilteredTrackIds = numel(filteredTrackIds);
meta.nTracksAll = nTracksAll;
meta.nTracksParsed = height(tracks);
meta.nEdgesParsed = height(edges);
meta.nSpotsParsed = height(spots);
end

function [trackIds, hasFilteredTracks] = parse_filtered_track_ids(doc)
trackIds = nan(0,1);
hasFilteredTracks = false;

filteredNodes = doc.getElementsByTagName('FilteredTracks');
if filteredNodes.getLength < 1
    return;
end

hasFilteredTracks = true;
filteredNode = filteredNodes.item(0);
trackIdNodes = filteredNode.getElementsByTagName('TrackID');

for i = 1:trackIdNodes.getLength
    n = trackIdNodes.item(i-1);
    tid = getAttrNum(n, 'TRACK_ID');
    if isfinite(tid)
        trackIds(end+1,1) = tid; %#ok<AGROW>
    end
end

if ~isempty(trackIds)
    trackIds = unique(trackIds, 'stable');
end
end

function xmlPath = sanitize_xml_path(xmlPath)
xmlPath = char(string(xmlPath));

% Normalize hidden control chars often introduced in batch/cluster path assembly.
xmlPath = strrep(xmlPath, sprintf('\r'), ' ');
xmlPath = strrep(xmlPath, sprintf('\n'), ' ');
xmlPath = strtrim(xmlPath);
end

function doc = read_trackmate_xml_document(xmlFile)
% Preflight readability before touching Java XML parsers.
fInfo = dir(xmlFile);
if isempty(fInfo)
    error('XML file metadata unavailable: %s', xmlFile);
end

[fid, fopenMsg] = fopen(xmlFile, 'r');
if fid < 0
    error('Cannot open XML for reading: %s\nfopen: %s', xmlFile, fopenMsg);
end
fclose(fid);

% Primary parser (MATLAB helper).
ME1 = [];
try
    doc = xmlread(xmlFile);
    return;
catch ME1
    % Fall through to Java parser fallback.
end

% Fallback parser (direct Java DOM builder).
try
    dbf = javax.xml.parsers.DocumentBuilderFactory.newInstance();
    dbf.setNamespaceAware(true);
    try
        dbf.setFeature('http://apache.org/xml/features/nonvalidating/load-external-dtd', false);
    catch
        % Feature may be unavailable in some JVM/parser builds.
    end
    db = dbf.newDocumentBuilder();
    doc = db.parse(java.io.File(xmlFile));
    return;
catch ME2
    headPreview = safe_text_head(xmlFile, 6);
    msg1 = exception_chain_to_text(ME1);
    msg2 = exception_chain_to_text(ME2);
    error(['xmlread failed for %s\n', ...
        'Primary parser:\n%s\n', ...
        'Fallback parser:\n%s\n', ...
        'File size: %.3f MB\n', ...
        'File head preview:\n%s'], ...
        xmlFile, msg1, msg2, fInfo.bytes / 1024 / 1024, headPreview);
end
end

function txt = safe_text_head(filePath, maxLines)
txt = "<unavailable>";
if nargin < 2 || ~isfinite(maxLines) || maxLines < 1
    maxLines = 6;
end

[fid, msg] = fopen(filePath, 'r');
if fid < 0
    txt = sprintf('<fopen failed: %s>', msg);
    return;
end

c = onCleanup(@() fclose(fid));
lines = strings(0,1);
for i = 1:maxLines
    ln = fgetl(fid);
    if ~ischar(ln)
        break;
    end
    lines(end+1,1) = string(ln); %#ok<AGROW>
end

if isempty(lines)
    txt = "<empty file or non-text prefix>";
else
    txt = char(strjoin(lines, newline));
end
end

function msg = exception_chain_to_text(ME)
if isempty(ME)
    msg = "<no exception>";
    return;
end

parts = strings(0,1);
parts(end+1,1) = string(ME.message); %#ok<AGROW>

for ci = 1:numel(ME.cause)
    c = ME.cause{ci};
    parts(end+1,1) = sprintf("Cause %d: %s", ci, string(c.message)); %#ok<AGROW>
end

stackTop = "";
if ~isempty(ME.stack)
    s = ME.stack(1);
    stackTop = sprintf(" (at %s:%d)", string(s.name), s.line);
end

msg = strjoin(parts, newline) + stackTop;
msg = char(msg);
end

function v = getAttrNum(node, attrName)
v = NaN;
if node.hasAttributes
    attrs = node.getAttributes;
    a = attrs.getNamedItem(attrName);
    if ~isempty(a)
        s = char(a.getValue);
        if strcmpi(s,'Infinity'); v = Inf; return; end
        if strcmpi(s,'NaN');      v = NaN; return; end
        tmp = str2double(s);
        if ~isnan(tmp) || strcmpi(strtrim(s),'0')
            v = tmp;
        end
    end
end
end
