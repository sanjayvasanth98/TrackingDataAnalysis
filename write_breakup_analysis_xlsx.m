function write_breakup_analysis_xlsx(allBreakup, xlsxFile)
%WRITE_BREAKUP_ANALYSIS_XLSX  Save breakup event data to an Excel file.
%   One sheet per case, with every child-parent pair as one row.
%
%   allBreakup: struct array with fields .caseName, optional .Re, .kD, .events

if isempty(allBreakup)
    warning('write_breakup_analysis_xlsx: no cases, nothing written.');
    return;
end

nCases = numel(allBreakup);

headers = {'Case','Re','kD','TrackID','ParentSpotID','ParentFrame', ...
    'ParentArea_px2','ParentAR','nChildren', ...
    'ChildSpotID','ChildFrame','ChildArea_px2','gamma','dRatio'};

% Delete existing file so old sheets don't persist.
if isfile(xlsxFile), delete(xlsxFile); end

usedSheetNames = strings(0,1);
for ci = 1:nCases
    ev = allBreakup(ci).events;
    reVal = get_case_re(allBreakup, ci);
    [sheetName, usedSheetNames] = make_unique_sheet_name(allBreakup(ci).caseName, reVal, ci, usedSheetNames);

    if isempty(ev)
        writecell(headers, xlsxFile, 'Sheet', sheetName, 'Range', 'A1');
        continue;
    end

    nEv  = numel(ev);
    rows = cell(nEv, numel(headers));
    for ei = 1:nEv
        rows(ei,:) = { ...
            char(allBreakup(ci).caseName), reVal, allBreakup(ci).kD, ...
            ev(ei).trackID,    ev(ei).parentID,   ev(ei).parentFrame, ...
            ev(ei).parentArea, ev(ei).parentAR,   ev(ei).nChildren, ...
            ev(ei).childID,    ev(ei).childFrame,  ev(ei).childArea, ...
            ev(ei).gamma,      ev(ei).dRatio};
    end

    T = cell2table(rows, 'VariableNames', headers);
    writetable(T, xlsxFile, 'Sheet', sheetName);
end

fprintf('Saved breakup analysis XLSX: %s\n', xlsxFile);
end


function reVal = get_case_re(allBreakup, ci)
reVal = NaN;
if isfield(allBreakup, 'Re') && numel(allBreakup) >= ci && ...
        ~isempty(allBreakup(ci).Re) && isfinite(allBreakup(ci).Re)
    reVal = double(allBreakup(ci).Re);
end
end


function [sheetName, usedSheetNames] = make_unique_sheet_name(caseName, reVal, ci, usedSheetNames)
if isfinite(reVal)
    base = sprintf('Re%g_%s', reVal, char(string(caseName)));
else
    base = char(string(caseName));
end
base = regexprep(base, '[:\\/\?\*\[\]]', '_');
base = regexprep(base, '\s+', '_');
base = strtrim(base);
if isempty(base)
    base = sprintf('Case%d', ci);
end

base = truncate_sheet_name(base, 31);
candidate = base;
suffixIdx = 1;
while any(strcmpi(usedSheetNames, string(candidate)))
    suffix = sprintf('_%d', suffixIdx);
    candidate = [truncate_sheet_name(base, 31 - numel(suffix)), suffix];
    suffixIdx = suffixIdx + 1;
end

sheetName = candidate;
usedSheetNames(end+1,1) = string(sheetName); %#ok<AGROW>
end


function out = truncate_sheet_name(in, maxLen)
if maxLen < 1
    out = '';
    return;
end
if numel(in) > maxLen
    out = in(1:maxLen);
else
    out = in;
end
end
