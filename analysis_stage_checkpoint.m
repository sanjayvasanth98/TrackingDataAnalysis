function varargout = analysis_stage_checkpoint(action, cacheFile, stageName, cacheKey, varargin)
%ANALYSIS_STAGE_CHECKPOINT Save/load derived-analysis checkpoints.
%
%   This is intentionally separate from the raw TrackMate XML parser cache.
%   It stores expensive downstream analysis results, such as full per-case
%   metrics or one XML's breakup-event extraction.

action = lower(char(string(action)));
switch action
    case 'file'
        varargout{1} = checkpoint_file(cacheFile, stageName, cacheKey);

    case 'load'
        checkpointFile = checkpoint_file(cacheFile, stageName, cacheKey);
        [ok, cachedData] = load_checkpoint(checkpointFile, cacheKey);
        varargout{1} = ok;
        varargout{2} = cachedData;
        varargout{3} = checkpointFile;

    case 'save'
        if isempty(varargin)
            error('analysis_stage_checkpoint:MissingData', ...
                'The save action requires cached data.');
        end
        cachedData = varargin{1};
        checkpointFile = checkpoint_file(cacheFile, stageName, cacheKey);
        ok = save_checkpoint(checkpointFile, cacheKey, cachedData);
        varargout{1} = ok;
        varargout{2} = checkpointFile;

    otherwise
        error('analysis_stage_checkpoint:UnknownAction', ...
            'Unknown checkpoint action: %s', action);
end
end


function checkpointFile = checkpoint_file(cacheFile, stageName, cacheKey)
cacheFile = char(string(cacheFile));
[cacheDir, ~, ~] = fileparts(cacheFile);
if isempty(cacheDir)
    cacheDir = pwd;
end

stageToken = sanitize_file_token(stageName);
if isempty(stageToken)
    stageToken = 'analysis';
end
checkpointDir = fullfile(cacheDir, [stageToken '_checkpoints']);

keyText = char(string(cacheKey));
caseToken = 'case';
tok = regexp(keyText, 'case=([^|]+)', 'tokens', 'once');
if ~isempty(tok)
    caseToken = tok{1};
else
    tok = regexp(keyText, 'xml=([^|]+)', 'tokens', 'once');
    if ~isempty(tok)
        [~, caseToken, ~] = fileparts(tok{1});
    end
end
caseToken = sanitize_file_token(caseToken);
if isempty(caseToken)
    caseToken = 'case';
end

checkpointFile = fullfile(checkpointDir, sprintf('%s_%s.mat', ...
    caseToken, cache_key_hash(keyText)));
end


function [ok, cachedData] = load_checkpoint(checkpointFile, expectedKey)
ok = false;
cachedData = [];
if ~isfile(checkpointFile)
    return;
end

try
    S = load(checkpointFile, 'cacheKey', 'cachedData');
catch ME
    warning('analysis_stage_checkpoint:LoadFailed', ...
        'Could not load analysis checkpoint %s: %s', checkpointFile, ME.message);
    return;
end

if ~isfield(S, 'cacheKey') || ~strcmp(string(S.cacheKey), string(expectedKey))
    return;
end
if ~isfield(S, 'cachedData')
    return;
end

cachedData = S.cachedData;
ok = true;
end


function ok = save_checkpoint(checkpointFile, cacheKey, cachedData)
ok = false;
[checkpointDir, ~, ~] = fileparts(checkpointFile);
if ~isfolder(checkpointDir)
    mkdir(checkpointDir);
end

cacheMeta = struct();
cacheMeta.savedAt = datestr(now, 31);
cacheMeta.checkpointVersion = 1;

tmpFile = [tempname(checkpointDir), '.mat'];
try
    save(tmpFile, 'cacheKey', 'cachedData', 'cacheMeta', '-v7.3');
    movefile(tmpFile, checkpointFile, 'f');
    ok = true;
catch ME
    warning('analysis_stage_checkpoint:SaveFailed', ...
        'Could not save analysis checkpoint %s: %s', checkpointFile, ME.message);
    if isfile(tmpFile)
        delete(tmpFile);
    end
end
end


function token = sanitize_file_token(token)
token = char(string(token));
token = regexprep(token, '[^A-Za-z0-9_.-]+', '_');
token = regexprep(token, '^_+', '');
token = regexprep(token, '_+$', '');
if numel(token) > 80
    token = token(1:80);
end
end


function hText = cache_key_hash(keyText)
h = 0;
for ii = 1:numel(keyText)
    h = mod(h * 131 + double(uint16(keyText(ii))), 2^32);
end
hText = sprintf('%08X', round(h));
end
