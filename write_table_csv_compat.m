function write_table_csv_compat(tbl, outCsv)
%WRITE_TABLE_CSV_COMPAT  Write a MATLAB table to CSV with old-version fallback.
%
% Uses writetable when available. If unavailable, falls back to manual CSV
% serialization for scalar table entries.

if exist('writetable', 'file') == 2
    try
        if ~requires_custom_csv_headers(tbl)
            writetable(tbl, outCsv);
            return;
        end
    catch
        % Fall back to manual writer below.
    end
end

fid = fopen(outCsv, 'w');
if fid < 0
    error('Cannot open CSV for writing: %s', outCsv);
end

varNames = tbl.Properties.VariableNames;
nCols = numel(varNames);
nRows = height(tbl);
headerNames = display_csv_header_names(varNames);

% Header
header = cell(1, nCols);
for c = 1:nCols
    header{c} = csv_escape(headerNames{c});
end
fprintf(fid, '%s\n', strjoin(header, ','));

% Rows
for r = 1:nRows
    rowVals = cell(1, nCols);
    for c = 1:nCols
        v = tbl{r, c};
        if iscell(v)
            if isempty(v)
                v = [];
            else
                v = v{1};
            end
        end
        rowVals{c} = csv_escape(to_scalar_text(v));
    end
    fprintf(fid, '%s\n', strjoin(rowVals, ','));
end

fclose(fid);
end

function tf = requires_custom_csv_headers(tbl)
tf = false;
if ~istable(tbl)
    return;
end
varNames = tbl.Properties.VariableNames;
customNames = {'leftMovingActivated_pct', 'AE_leftMoving_pct'};
tf = any(ismember(varNames, customNames));
end

function headerNames = display_csv_header_names(varNames)
headerNames = varNames;
for c = 1:numel(headerNames)
    switch headerNames{c}
        case 'leftMovingActivated_pct'
            headerNames{c} = 'leftMovingActivated %';
        case 'AE_leftMoving_pct'
            headerNames{c} = 'AE_leftmoving %';
    end
end
end

function s = to_scalar_text(v)
if isempty(v)
    s = '';
    return;
end

if isstring(v)
    if numel(v) == 1
        s = char(v);
    else
        s = char(strjoin(v(:).', ';'));
    end
    return;
end

if ischar(v)
    s = v;
    return;
end

if isnumeric(v) || islogical(v)
    if isscalar(v)
        s = sprintf('%.15g', double(v));
    else
        parts = arrayfun(@(x) sprintf('%.15g', double(x)), v(:), 'UniformOutput', false);
        s = strjoin(parts, ';');
    end
    return;
end

if has_function('isdatetime') && isdatetime(v)
    if isscalar(v)
        s = char(v);
    else
        s = char(strjoin(string(v(:).'), ';'));
    end
    return;
end

if has_function('isduration') && isduration(v)
    s = char(string(v));
    return;
end

if has_function('iscategorical') && iscategorical(v)
    s = char(string(v));
    return;
end

% Generic fallback
try
    s = char(string(v));
catch
    s = '<unsupported>';
end
end

function out = csv_escape(in)
if isstring(in)
    in = char(in);
end
if isempty(in)
    out = '';
    return;
end

out = strrep(in, '"', '""');
if any(out == ',') || any(out == '"') || any(out == sprintf('\n')) || any(out == sprintf('\r'))
    out = ['"' out '"'];
end
end

function tf = has_function(fnName)
tf = (exist(fnName, 'builtin') == 5) || (exist(fnName, 'file') == 2);
end
