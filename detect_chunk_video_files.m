function videoFiles = detect_chunk_video_files(videoPathOrBase)
%DETECT_CHUNK_VIDEO_FILES  Find numbered chunk videos from a folder or base video.
%
%   videoFiles = detect_chunk_video_files('path/to/case.avi')
%   videoFiles = detect_chunk_video_files('path/to/case folder')
%
%   Preferred behavior:
%     - If numbered files exist, return only:
%         stem_1.avi, stem_2.avi, ...
%     - Otherwise, if the unnumbered base exists, return that single file.
%     - If nothing exists or the input is empty, return {}.
%
%   If the input is a folder, the folder name is used as the base stem and
%   the search is performed inside that folder. For example:
%     input folder:  E:\...\P10S100 2000\
%     chunk videos:  E:\...\P10S100 2000\P10S100 2000_1.avi
%                    E:\...\P10S100 2000\P10S100 2000_2.avi

videoFiles = {};

if nargin < 1
    return;
end

videoPathOrBase = char(string(videoPathOrBase));
if isempty(strtrim(videoPathOrBase))
    return;
end

[videoDir, baseName, ext] = resolve_video_search_root(videoPathOrBase);
if isempty(videoDir) || isempty(baseName)
    return;
end

listing = [ ...
    dir(fullfile(videoDir, [baseName '_*.avi'])); ...
    dir(fullfile(videoDir, [baseName '_*.AVI']))];
chunkFiles = {};
chunkNums = [];
escapedBase = regexptranslate('escape', baseName);
pattern = ['^' escapedBase '_(\d+)\.(avi|AVI)$'];

for k = 1:numel(listing)
    tok = regexp(listing(k).name, pattern, 'tokens');
    if ~isempty(tok)
        chunkFiles{end+1} = fullfile(videoDir, listing(k).name); %#ok<AGROW>
        chunkNums(end+1) = str2double(tok{1}{1}); %#ok<AGROW>
    end
end

if ~isempty(chunkFiles)
    [~, sortIdx] = sort(chunkNums);
    videoFiles = chunkFiles(sortIdx);
    return;
end

baseCandidates = { ...
    fullfile(videoDir, [baseName '.avi']), ...
    fullfile(videoDir, [baseName '.AVI'])};
for i = 1:numel(baseCandidates)
    if isfile(baseCandidates{i})
        videoFiles = baseCandidates(i);
        return;
    end
end
end

function [videoDir, baseName, ext] = resolve_video_search_root(videoPathOrBase)
videoDir = '';
baseName = '';
ext = '';

if isfolder(videoPathOrBase)
    folderPath = char(string(videoPathOrBase));
    folderPath = regexprep(folderPath, '[\\/]+$', '');
    [~, folderName] = fileparts(folderPath);
    videoDir = folderPath;
    baseName = folderName;
    return;
end

[videoDir, baseName, ext] = fileparts(videoPathOrBase);
if isempty(videoDir)
    videoDir = '.';
end
if isempty(ext)
    return;
end
end
