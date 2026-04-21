function caseKey = build_case_analysis_cache_key(caseDef, xmlFiles, maxTracks, policyTag)
%BUILD_CASE_ANALYSIS_CACHE_KEY Cache key for a full derived case analysis.

if nargin < 4 || isempty(policyTag)
    policyTag = "default";
end

if nargin < 2 || isempty(xmlFiles)
    xmlFiles = {caseDef.xmlFile};
end
xmlFiles = cellstr(string(xmlFiles(:)));

maxTracksStr = "all";
if isfinite(maxTracks)
    maxTracksStr = string(max(0, floor(maxTracks)));
end

xmlParts = strings(numel(xmlFiles), 1);
for i = 1:numel(xmlFiles)
    xmlParts(i) = sprintf("%s@%s", string(xmlFiles{i}), file_signature(xmlFiles{i}));
end

caseKey = sprintf("case=%s|Re=%.12g|kD=%.12g|px=%.12g|dt=%.12g|max=%s|xmls=%s|policy=%s", ...
    string(caseDef.name), caseDef.Re, caseDef.kD, caseDef.pixelSize, caseDef.dt, ...
    maxTracksStr, strjoin(xmlParts, ";;"), string(policyTag));
end


function sig = file_signature(filePath)
if isfile(filePath)
    d = dir(filePath);
    sig = sprintf("%d_%0.0f", d.bytes, d.datenum * 1e6);
else
    sig = "missing";
end
end
