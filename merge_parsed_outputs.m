function outMerged = merge_parsed_outputs(outs, xmlFiles)
%MERGE_PARSED_OUTPUTS  Merge parsed TrackMate XML outputs into one case output.
%
%   outMerged = merge_parsed_outputs(outs, xmlFiles)
%
%   Merges the parsed outputs returned by analyze_trackmate_xml_arc so the
%   combined case behaves like one longer sequential record. Spot and track
%   identifiers are offset to avoid collisions across chunks, while frame
%   and time axes are offset so later chunks follow earlier ones.

if nargin < 1 || isempty(outs)
    error('merge_parsed_outputs:MissingInput', ...
        'At least one parsed output struct is required.');
end

nChunks = numel(outs);
if nargin < 2 || isempty(xmlFiles)
    xmlFiles = repmat({''}, nChunks, 1);
elseif numel(xmlFiles) ~= nChunks
    error('merge_parsed_outputs:InputSizeMismatch', ...
        'xmlFiles must have the same length as outs.');
end
xmlFiles = cellfun(@char, cellstr(string(xmlFiles(:))), 'UniformOutput', false);

baseOut = outs{1};
allSpots = baseOut.spots([],:);
allTracks = baseOut.tracks([],:);
allEdges = baseOut.edges([],:);
allTraj = baseOut.trajectories([]);

chunkTemplate = struct( ...
    'chunkIndex', NaN, ...
    'xmlFile', "", ...
    'nSpots', 0, ...
    'nTracks', 0, ...
    'frameMinOriginal', NaN, ...
    'frameMaxOriginal', NaN, ...
    'frameOffset', 0, ...
    'timeMinOriginal', NaN, ...
    'timeMaxOriginal', NaN, ...
    'timeOffset', 0, ...
    'timeStepUsed', NaN, ...
    'spotIdOffset', 0, ...
    'trackIdOffset', 0);
chunkInfo = repmat(chunkTemplate, nChunks, 1);

nextSpotIdOffset = 0;
nextTrackIdOffset = 0;
combinedFrameMax = NaN;
combinedTimeMax = NaN;
lastPositiveTimeStep = NaN;

for c = 1:nChunks
    outChunk = outs{c};
    if ~isstruct(outChunk)
        error('merge_parsed_outputs:InvalidChunk', ...
            'Chunk %d is not a parsed output struct.', c);
    end

    spots = outChunk.spots;
    tracks = outChunk.tracks;
    edges = outChunk.edges;
    traj = outChunk.trajectories;

    frameVals = collect_numeric_field(spots, 'FRAME');
    frameVals = [frameVals; collect_traj_field(traj, 'frame')];
    frameVals = frameVals(isfinite(frameVals));
    frameMin = NaN;
    frameMax = NaN;
    if ~isempty(frameVals)
        frameMin = min(frameVals);
        frameMax = max(frameVals);
    end

    timeVals = collect_numeric_field(spots, 'POS_T');
    timeVals = [timeVals; collect_traj_field(traj, 't')];
    timeVals = timeVals(isfinite(timeVals));
    timeMin = NaN;
    timeMax = NaN;
    if ~isempty(timeVals)
        timeMin = min(timeVals);
        timeMax = max(timeVals);
    end

    frameOffset = 0;
    if c > 1 && isfinite(combinedFrameMax)
        if isfinite(frameMin)
            frameOffset = combinedFrameMax + 1 - frameMin;
        else
            frameOffset = combinedFrameMax + 1;
        end
    end

    timeStep = estimate_time_step(outChunk);
    if ~isfinite(timeStep) || timeStep <= 0
        timeStep = lastPositiveTimeStep;
    end
    if ~isfinite(timeStep) || timeStep <= 0
        timeStep = 0;
    end

    timeOffset = 0;
    if c > 1 && isfinite(combinedTimeMax)
        if isfinite(timeMin)
            timeOffset = combinedTimeMax + timeStep - timeMin;
        else
            timeOffset = combinedTimeMax + timeStep;
        end
    end

    spotIdOffset = nextSpotIdOffset;
    trackIdOffset = nextTrackIdOffset;

    if ~isempty(spots) && height(spots) > 0
        if ismember('ID', spots.Properties.VariableNames)
            spots.ID = spots.ID + spotIdOffset;
        end
        if ismember('FRAME', spots.Properties.VariableNames)
            spots.FRAME = spots.FRAME + frameOffset;
        end
        if ismember('POS_T', spots.Properties.VariableNames)
            spots.POS_T = spots.POS_T + timeOffset;
        end
    end

    for k = 1:numel(traj)
        if isfield(traj, 'TRACK_ID') && isfinite(traj(k).TRACK_ID)
            traj(k).TRACK_ID = traj(k).TRACK_ID + trackIdOffset;
        end
        if isfield(traj, 'spotIds') && ~isempty(traj(k).spotIds)
            traj(k).spotIds = traj(k).spotIds + spotIdOffset;
        end
        if isfield(traj, 'frame') && ~isempty(traj(k).frame)
            traj(k).frame = traj(k).frame + frameOffset;
        end
        if isfield(traj, 't') && ~isempty(traj(k).t)
            traj(k).t = traj(k).t + timeOffset;
        end
    end

    if ~isempty(tracks) && height(tracks) > 0
        if ismember('TRACK_ID', tracks.Properties.VariableNames)
            tracks.TRACK_ID = tracks.TRACK_ID + trackIdOffset;
        end
        if ismember('TRACK_START', tracks.Properties.VariableNames)
            tracks.TRACK_START = tracks.TRACK_START + frameOffset;
        end
        if ismember('TRACK_STOP', tracks.Properties.VariableNames)
            tracks.TRACK_STOP = tracks.TRACK_STOP + frameOffset;
        end
    end

    if ~isempty(edges) && height(edges) > 0
        if ismember('TRACK_ID', edges.Properties.VariableNames)
            edges.TRACK_ID = edges.TRACK_ID + trackIdOffset;
        end
        if ismember('SPOT_SOURCE_ID', edges.Properties.VariableNames)
            edges.SPOT_SOURCE_ID = edges.SPOT_SOURCE_ID + spotIdOffset;
        end
        if ismember('SPOT_TARGET_ID', edges.Properties.VariableNames)
            edges.SPOT_TARGET_ID = edges.SPOT_TARGET_ID + spotIdOffset;
        end
        if ismember('EDGE_TIME', edges.Properties.VariableNames)
            edges.EDGE_TIME = edges.EDGE_TIME + timeOffset;
        end
    end

    allSpots = [allSpots; spots]; %#ok<AGROW>
    allTracks = [allTracks; tracks]; %#ok<AGROW>
    allEdges = [allEdges; edges]; %#ok<AGROW>
    allTraj = [allTraj; traj(:)]; %#ok<AGROW>

    if ~isempty(spots) && height(spots) > 0 && ismember('ID', spots.Properties.VariableNames)
        nextSpotIdOffset = max(spots.ID) + 1;
    end

    trackIds = nan(0,1);
    if ~isempty(tracks) && height(tracks) > 0 && ismember('TRACK_ID', tracks.Properties.VariableNames)
        trackIds = [trackIds; tracks.TRACK_ID];
    end
    trajTrackIds = collect_traj_scalar_field(traj, 'TRACK_ID');
    if ~isempty(trajTrackIds)
        trackIds = [trackIds; trajTrackIds];
    end
    trackIds = trackIds(isfinite(trackIds));
    if ~isempty(trackIds)
        nextTrackIdOffset = max(trackIds) + 1;
    end

    if isfinite(frameMax)
        combinedFrameMax = frameMax + frameOffset;
    end
    if isfinite(timeMax)
        combinedTimeMax = timeMax + timeOffset;
    end
    if isfinite(timeStep) && timeStep > 0
        lastPositiveTimeStep = timeStep;
    end

    chunkInfo(c).chunkIndex = c;
    chunkInfo(c).xmlFile = string(xmlFiles{c});
    chunkInfo(c).nSpots = height(spots);
    chunkInfo(c).nTracks = numel(traj);
    chunkInfo(c).frameMinOriginal = frameMin;
    chunkInfo(c).frameMaxOriginal = frameMax;
    chunkInfo(c).frameOffset = frameOffset;
    chunkInfo(c).timeMinOriginal = timeMin;
    chunkInfo(c).timeMaxOriginal = timeMax;
    chunkInfo(c).timeOffset = timeOffset;
    chunkInfo(c).timeStepUsed = timeStep;
    chunkInfo(c).spotIdOffset = spotIdOffset;
    chunkInfo(c).trackIdOffset = trackIdOffset;
end

outMerged = baseOut;
outMerged.spots = allSpots;
outMerged.tracks = allTracks;
outMerged.edges = allEdges;
outMerged.trajectories = allTraj;
if ~isfield(outMerged, 'meta') || ~isstruct(outMerged.meta)
    outMerged.meta = struct();
end
outMerged.meta.nChunksMerged = nChunks;
outMerged.meta.sourceXmlFiles = string(xmlFiles);
outMerged.meta.chunkInfo = chunkInfo;
outMerged.meta.isMergedSequentialChunks = true;
end

function vals = collect_numeric_field(tbl, fieldName)
vals = nan(0,1);
if isempty(tbl) || ~istable(tbl) || height(tbl) == 0
    return;
end
if ~ismember(fieldName, tbl.Properties.VariableNames)
    return;
end
vals = tbl.(fieldName);
vals = vals(:);
end

function vals = collect_traj_field(traj, fieldName)
vals = nan(0,1);
if isempty(traj) || ~isstruct(traj) || ~isfield(traj, fieldName)
    return;
end
for k = 1:numel(traj)
    raw = traj(k).(fieldName);
    if isempty(raw)
        continue;
    end
    vals = [vals; raw(:)]; %#ok<AGROW>
end
end

function vals = collect_traj_scalar_field(traj, fieldName)
vals = nan(0,1);
if isempty(traj) || ~isstruct(traj) || ~isfield(traj, fieldName)
    return;
end
vals = nan(numel(traj), 1);
for k = 1:numel(traj)
    raw = traj(k).(fieldName);
    if isempty(raw)
        continue;
    end
    vals(k) = raw(1);
end
end

function dtVal = estimate_time_step(outChunk)
dtVal = NaN;

dtCandidates = nan(0,1);

if isfield(outChunk, 'trajectories') && ~isempty(outChunk.trajectories)
    for k = 1:numel(outChunk.trajectories)
        tr = outChunk.trajectories(k);
        if isfield(tr, 'dt') && ~isempty(tr.dt)
            dtCandidates = [dtCandidates; tr.dt(:)]; %#ok<AGROW>
        elseif isfield(tr, 't') && numel(tr.t) >= 2
            dtCandidates = [dtCandidates; diff(tr.t(:))]; %#ok<AGROW>
        end
    end
end

if isempty(dtCandidates) && isfield(outChunk, 'spots') && istable(outChunk.spots) && ...
        ismember('POS_T', outChunk.spots.Properties.VariableNames)
    tVals = unique(outChunk.spots.POS_T(:));
    tVals = tVals(isfinite(tVals));
    if numel(tVals) >= 2
        dtCandidates = diff(sort(tVals));
    end
end

dtCandidates = dtCandidates(isfinite(dtCandidates) & dtCandidates > 0);
if isempty(dtCandidates)
    return;
end

dtVal = min(dtCandidates);
end
