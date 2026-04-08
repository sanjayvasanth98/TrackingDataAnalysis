function xmlFiles = detect_chunk_xml_files(xmlPathOrBase)
%DETECT_CHUNK_XML_FILES  Find numbered chunk XMLs from a folder or base XML.
%
%   xmlFiles = detect_chunk_xml_files('path/to/P10S100_48lit.xml')
%   xmlFiles = detect_chunk_xml_files('path/to/P10S100 2000')
%
%   Given a base XML path or a folder path, auto-detects numbered chunk
%   files in the target directory using the pattern <baseName>_N.xml
%   (N = 1, 2, 3, ...).
%
%   Preferred behavior:
%     - If numbered chunk files exist, return only those:
%         { '.../P10S100 2000_1.xml',
%           '.../P10S100 2000_2.xml', ... }
%     - If no numbered chunk files exist, fall back to the unnumbered base:
%         { '.../P10S100 2000.xml' }
%
%   If the input is a folder, the folder name is used as the XML stem and
%   the search is performed inside that folder. If no exact base XML is
%   found but there is exactly one XML in the folder, that XML is used as a
%   single-file fallback.

xmlPathOrBase = char(string(xmlPathOrBase));
if isempty(strtrim(xmlPathOrBase))
    error('detect_chunk_xml_files:MissingPath', 'A base XML path is required.');
end
[xmlDir, baseName] = resolve_xml_search_root(xmlPathOrBase);
if isempty(xmlDir) || isempty(baseName)
    error('detect_chunk_xml_files:InvalidPath', ...
        'Expected an XML path or folder path: %s', xmlPathOrBase);
end

listing = [ ...
    dir(fullfile(xmlDir, [baseName '_*.xml'])); ...
    dir(fullfile(xmlDir, [baseName '_*.XML']))];

chunkFiles = {};
chunkNums  = [];
escapedBase = regexptranslate('escape', baseName);
pattern = ['^' escapedBase '_(\d+)\.(xml|XML)$'];

for k = 1:numel(listing)
    tok = regexp(listing(k).name, pattern, 'tokens');
    if ~isempty(tok)
        chunkFiles{end+1} = fullfile(xmlDir, listing(k).name); %#ok<AGROW>
        chunkNums(end+1)   = str2double(tok{1}{1});             %#ok<AGROW>
    end
end

% Sort by chunk number
[~, sortIdx] = sort(chunkNums);
chunkFiles = chunkFiles(sortIdx);

if ~isempty(chunkFiles)
    xmlFiles = chunkFiles;
    return;
end

baseCandidates = { ...
    fullfile(xmlDir, [baseName '.xml']), ...
    fullfile(xmlDir, [baseName '.XML'])};
for i = 1:numel(baseCandidates)
    if isfile(baseCandidates{i})
        xmlFiles = baseCandidates(i);
        return;
    end
end

allXml = [dir(fullfile(xmlDir, '*.xml')); dir(fullfile(xmlDir, '*.XML'))];
if numel(allXml) == 1
    xmlFiles = {fullfile(allXml(1).folder, allXml(1).name)};
    return;
end

error('detect_chunk_xml_files:NoChunkFilesFound', ...
    'No numbered chunk XMLs and no usable base XML found for: %s', xmlPathOrBase);
end

function [xmlDir, baseName] = resolve_xml_search_root(xmlPathOrBase)
xmlDir = '';
baseName = '';

if isfolder(xmlPathOrBase)
    folderPath = char(string(xmlPathOrBase));
    folderPath = regexprep(folderPath, '[\\/]+$', '');
    [~, folderName] = fileparts(folderPath);
    xmlDir = folderPath;
    baseName = folderName;
    return;
end

[xmlDir, baseName, ext] = fileparts(xmlPathOrBase);
if isempty(xmlDir)
    xmlDir = '.';
end
if ~strcmpi(ext, '.xml')
    xmlDir = '';
    baseName = '';
end
end
