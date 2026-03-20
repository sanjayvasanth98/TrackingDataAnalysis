function [selectedIdx, selectedTrackIds] = resolve_diagnostic_track_indices(trackCatalog, trackRequests, caseDef, callerTag)
selectedIdx = nan(0,1);
selectedTrackIds = nan(0,1);

if nargin < 1 || isempty(trackCatalog)
    return;
end
if nargin < 2 || isempty(trackRequests)
    trackRequests = [];
end
if nargin < 3 || isempty(caseDef)
    caseDef = struct('name', "case");
end
if nargin < 4 || isempty(callerTag)
    callerTag = 'diagnostic';
end

trackIds = [trackCatalog.TRACK_ID];
caseName = "case";
if isstruct(caseDef) && isfield(caseDef, 'name')
    caseName = string(caseDef.name);
end

if isempty(trackRequests)
    selectedIdx = (1:numel(trackCatalog)).';
    selectedTrackIds = trackIds(:);
    return;
end

trackRequests = trackRequests(:).';
for req = trackRequests
    idx = find(trackIds == req, 1, 'first');
    if isempty(idx) && isfinite(req) && req >= 1 && req <= numel(trackCatalog) && abs(req - round(req)) < 1e-9
        idx = round(req);
    end

    if isempty(idx)
        warning('%s: requested track %g not found for case %s.', callerTag, req, char(caseName));
        continue;
    end

    selectedIdx(end+1,1) = idx; %#ok<AGROW>
end

selectedIdx = unique(selectedIdx, 'stable');
if isempty(selectedIdx)
    return;
end
selectedTrackIds = trackIds(selectedIdx).';
end
