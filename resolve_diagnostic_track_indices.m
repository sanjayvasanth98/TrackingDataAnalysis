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
excludedMask = false(numel(trackCatalog), 1);
for i = 1:numel(trackCatalog)
    % Exclude origin-box tracks
    if isfield(trackCatalog, 'isOriginExcludedBox')
        v = trackCatalog(i).isOriginExcludedBox;
        if (islogical(v) && any(v(:))) || (isnumeric(v) && any(v(:) ~= 0))
            excludedMask(i) = true;
            continue;
        end
    end
    % Exclude tracks that are not basic-valid (short, non-finite, etc.)
    if isfield(trackCatalog, 'isBasicValid')
        v = trackCatalog(i).isBasicValid;
        if (islogical(v) && ~any(v(:))) || (isnumeric(v) && ~any(v(:) ~= 0))
            excludedMask(i) = true;
        end
    end
end
caseName = "case";
if isstruct(caseDef) && isfield(caseDef, 'name')
    caseName = string(caseDef.name);
end

if isempty(trackRequests)
    selectedIdx = find(~excludedMask);
    selectedTrackIds = trackIds(:);
    selectedTrackIds = selectedTrackIds(selectedIdx);
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
if any(excludedMask(selectedIdx))
    nDrop = sum(excludedMask(selectedIdx));
    warning('%s: dropped %d requested track(s) excluded by origin box for case %s.', callerTag, nDrop, char(caseName));
    selectedIdx = selectedIdx(~excludedMask(selectedIdx));
end
if isempty(selectedIdx)
    return;
end
selectedTrackIds = trackIds(selectedIdx).';
end
