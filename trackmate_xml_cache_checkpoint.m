function varargout = trackmate_xml_cache_checkpoint(action, cacheFile, caseKey, varargin)
%TRACKMATE_XML_CACHE_CHECKPOINT  Save/load one parsed TrackMate XML checkpoint.
%
%   This is intentionally one checkpoint per XML cache key.  If a long ARC
%   run dies after parsing sample 1 of N, the next run can resume from the
%   sample-level checkpoint instead of waiting for the end-of-script cache.

action = lower(char(string(action)));
switch action
    case 'file'
        varargout{1} = checkpoint_file(cacheFile, caseKey);

    case 'load'
        parserOpts = struct();
        if ~isempty(varargin)
            parserOpts = varargin{1};
        end
        checkpointFile = checkpoint_file(cacheFile, caseKey);
        [ok, cachedOut] = load_checkpoint(checkpointFile, caseKey, parserOpts);
        varargout{1} = ok;
        varargout{2} = cachedOut;
        varargout{3} = checkpointFile;

    case 'save'
        if isempty(varargin)
            error('trackmate_xml_cache_checkpoint:MissingOutput', ...
                'The save action requires the parsed output struct.');
        end
        parsedOut = varargin{1};
        checkpointFile = checkpoint_file(cacheFile, caseKey);
        ok = save_checkpoint(checkpointFile, caseKey, parsedOut);
        varargout{1} = ok;
        varargout{2} = checkpointFile;

    otherwise
        error('trackmate_xml_cache_checkpoint:UnknownAction', ...
            'Unknown checkpoint action: %s', action);
end
end


function checkpointFile = checkpoint_file(cacheFile, caseKey)
cacheFile = char(string(cacheFile));
[cacheDir, ~, ~] = fileparts(cacheFile);
if isempty(cacheDir)
    cacheDir = pwd;
end
checkpointDir = fullfile(cacheDir, 'xml_parse_checkpoints');

keyText = char(string(caseKey));
xmlStem = 'xml';
tok = regexp(keyText, 'xml=([^|]+)', 'tokens', 'once');
if ~isempty(tok)
    [~, xmlStem, ~] = fileparts(tok{1});
end
xmlStem = sanitize_file_token(xmlStem);
if isempty(xmlStem)
    xmlStem = 'xml';
end

checkpointFile = fullfile(checkpointDir, sprintf('%s_%s.mat', xmlStem, cache_key_hash(keyText)));
end


function [ok, cachedOut] = load_checkpoint(checkpointFile, expectedKey, parserOpts)
ok = false;
cachedOut = [];
if ~isfile(checkpointFile)
    return;
end

try
    S = load(checkpointFile, 'caseKey', 'cachedOut');
catch ME
    warning('trackmate_xml_cache_checkpoint:LoadFailed', ...
        'Could not load XML checkpoint %s: %s', checkpointFile, ME.message);
    return;
end

if ~isfield(S, 'caseKey') || ~strcmp(string(S.caseKey), string(expectedKey))
    return;
end
if ~isfield(S, 'cachedOut')
    return;
end
if ~is_cache_entry_compatible_checkpoint(S.cachedOut, parserOpts)
    return;
end

cachedOut = S.cachedOut;
ok = true;
end


function ok = save_checkpoint(checkpointFile, caseKey, cachedOut)
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
    save(tmpFile, 'caseKey', 'cachedOut', 'cacheMeta', '-v7.3');
    movefile(tmpFile, checkpointFile, 'f');
    ok = true;
catch ME
    warning('trackmate_xml_cache_checkpoint:SaveFailed', ...
        'Could not save XML checkpoint %s: %s', checkpointFile, ME.message);
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


function tf = is_cache_entry_compatible_checkpoint(out, parserOpts)
tf = false;
if ~isstruct(out) || ~isfield(out, 'meta') || ~isstruct(out.meta)
    return;
end

meta = out.meta;
if isfield(parserOpts, 'parserVersion')
    if ~isfield(meta, 'parserVersion') || ~isfinite(meta.parserVersion)
        return;
    end
    if meta.parserVersion < parserOpts.parserVersion
        return;
    end
end

if isfield(parserOpts, 'parseTrackedSpotsOnly') && parserOpts.parseTrackedSpotsOnly
    if ~isfield(meta, 'parseTrackedSpotsOnly') || ~is_true_flag_checkpoint(meta.parseTrackedSpotsOnly)
        return;
    end
end

if isfield(parserOpts, 'parseFilteredTracksOnly') && parserOpts.parseFilteredTracksOnly
    if ~isfield(meta, 'parseFilteredTracksOnly') || ~is_true_flag_checkpoint(meta.parseFilteredTracksOnly)
        return;
    end
end

tf = true;
end


function tf = is_true_flag_checkpoint(v)
tf = false;
if islogical(v)
    tf = any(v(:));
elseif isnumeric(v)
    tf = any(v(:) ~= 0);
elseif isstring(v) || ischar(v)
    tf = any(strcmpi(string(v), ["true", "1", "yes"]));
end
end
