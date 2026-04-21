function sig = cache_value_signature(v)
%CACHE_VALUE_SIGNATURE Compact deterministic signature for cache keys.
%
%   Large numeric arrays, such as ROI masks, are summarized rather than fully
%   expanded so cache keys stay small while still changing when key settings
%   or masks change.

sig = signature_value(v, 0);
end


function sig = signature_value(v, depth)
if nargin < 2
    depth = 0;
end
if depth > 12
    sig = sprintf('<max-depth:%s>', class(v));
    return;
end

if isstruct(v)
    sig = signature_struct(v, depth);
elseif istable_compat(v)
    sig = sprintf('table[%s]h%dw%d', class(v), size(v, 1), size(v, 2));
elseif iscell(v)
    sig = signature_cell(v, depth);
elseif isnumeric(v) || islogical(v)
    sig = signature_numeric(v);
elseif isstring(v)
    sig = ['string:', strjoin(cellstr(v(:)), '|')];
elseif ischar(v)
    sig = ['char:', v];
elseif isa(v, 'function_handle')
    sig = ['function:', func2str(v)];
else
    sig = ['class:', class(v)];
end
end


function sig = signature_struct(s, depth)
if isempty(s)
    sig = sprintf('struct:%s:empty', class(s));
    return;
end
names = sort(fieldnames(s));
parts = cell(numel(names) + 1, 1);
parts{1} = sprintf('struct[%s]', size_text(size(s)));
for i = 1:numel(names)
    name = names{i};
    vals = cell(numel(s), 1);
    for k = 1:numel(s)
        vals{k} = signature_value(s(k).(name), depth + 1);
    end
    parts{i + 1} = [name '=' strjoin(vals, ',')];
end
sig = strjoin(parts, ';');
end


function sig = signature_cell(c, depth)
if isempty(c)
    sig = sprintf('cell[%s]:empty', size_text(size(c)));
    return;
end
n = numel(c);
nKeep = min(n, 32);
parts = cell(nKeep + 1, 1);
parts{1} = sprintf('cell[%s]n%d', size_text(size(c)), n);
for i = 1:nKeep
    parts{i + 1} = signature_value(c{i}, depth + 1);
end
if nKeep < n
    parts{end} = sprintf('%s;truncated=%d', parts{end}, n - nKeep);
end
sig = strjoin(parts, ';');
end


function sig = signature_numeric(v)
sz = size(v);
cls = class(v);
if isempty(v)
    sig = sprintf('%s[%s]:empty', cls, size_text(sz));
    return;
end

if numel(v) <= 64
    sig = sprintf('%s[%s]:%s', cls, size_text(sz), mat2str(v));
    return;
end

if islogical(v)
    numericVals = double(v(:));
else
    numericVals = double(v(:));
end
finiteVals = numericVals(isfinite(numericVals));
if isempty(finiteVals)
    sig = sprintf('%s[%s]:n%d;finite=0;nnz=%d', ...
        cls, size_text(sz), numel(v), nnz(numericVals));
    return;
end

sig = sprintf('%s[%s]:n%d;finite=%d;nnz=%d;sum=%.17g;mean=%.17g;min=%.17g;max=%.17g', ...
    cls, size_text(sz), numel(v), numel(finiteVals), nnz(numericVals), ...
    sum(finiteVals), mean(finiteVals), min(finiteVals), max(finiteVals));
end


function txt = size_text(sz)
txt = sprintf('%dx', sz);
if ~isempty(txt)
    txt = txt(1:end-1);
end
end


function tf = istable_compat(v)
tf = false;
try
    tf = istable(v);
catch
    tf = false;
end
end
