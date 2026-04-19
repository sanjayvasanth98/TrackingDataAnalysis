function csvFiles = write_breakup_analysis_csv(allBreakup, outPath)
%WRITE_BREAKUP_ANALYSIS_CSV  Save breakup event data as CSV files.
%
%   Writes one combined CSV plus one CSV per case.  This replaces the old
%   Excel-sheet export because some ARC MATLAB environments do not provide
%   writetable/writecell.
%
%   If outPath is ".../breakup_events.csv":
%     .../breakup_events.csv
%     .../breakup_events_by_case/<case>.csv
%
%   If outPath is an old ".../breakup_events.xlsx" path, the extension is
%   changed to CSV and the same CSV outputs are written.

csvFiles = strings(0,1);
if isempty(allBreakup)
    warning('write_breakup_analysis_csv: no cases, nothing written.');
    return;
end

[outDir, stem, ext] = fileparts(char(string(outPath)));
if isempty(outDir)
    outDir = pwd;
end
if isempty(stem)
    stem = 'breakup_events';
end
if isempty(ext)
    combinedCsv = fullfile(outDir, [stem, '.csv']);
else
    combinedCsv = fullfile(outDir, [stem, '.csv']);
end
caseDir = fullfile(outDir, [stem, '_by_case']);
if ~isfolder(outDir), mkdir(outDir); end
if ~isfolder(caseDir), mkdir(caseDir); end

allRows = table();
usedCaseTokens = strings(0,1);
for ci = 1:numel(allBreakup)
    caseTbl = breakup_case_to_table(allBreakup, ci);
    allRows = append_table_compat_local(allRows, caseTbl);

    [caseToken, usedCaseTokens] = make_unique_case_token(allBreakup(ci).caseName, ...
        get_case_re(allBreakup, ci), ci, usedCaseTokens);
    caseCsv = fullfile(caseDir, [caseToken, '.csv']);
    write_table_csv_compat(caseTbl, caseCsv);
    csvFiles(end+1,1) = string(caseCsv); %#ok<AGROW>
end

write_table_csv_compat(allRows, combinedCsv);
csvFiles = [string(combinedCsv); csvFiles(:)];
fprintf('Saved breakup analysis CSVs: %s and %d per-case file(s) in %s\n', ...
    combinedCsv, numel(csvFiles) - 1, caseDir);
end


function T = breakup_case_to_table(allBreakup, ci)
headers = {'Case','Re','kD','TrackID','ParentSpotID','ParentFrame', ...
    'ParentArea_px2','ParentAR','nChildren', ...
    'ChildSpotID','ChildFrame','ChildArea_px2','gamma','dRatio'};

ev = allBreakup(ci).events;
reVal = get_case_re(allBreakup, ci);
if isempty(ev)
    T = cell2table(cell(0, numel(headers)), 'VariableNames', headers);
    return;
end

nEv  = numel(ev);
rows = cell(nEv, numel(headers));
for ei = 1:nEv
    rows(ei,:) = { ...
        char(allBreakup(ci).caseName), reVal, allBreakup(ci).kD, ...
        get_event_field(ev(ei), 'trackID'), ...
        get_event_field(ev(ei), 'parentID'), ...
        get_event_field(ev(ei), 'parentFrame'), ...
        get_event_field(ev(ei), 'parentArea'), ...
        get_event_field(ev(ei), 'parentAR'), ...
        get_event_field(ev(ei), 'nChildren'), ...
        get_event_field(ev(ei), 'childID'), ...
        get_event_field(ev(ei), 'childFrame'), ...
        get_event_field(ev(ei), 'childArea'), ...
        get_event_field(ev(ei), 'gamma'), ...
        get_event_field(ev(ei), 'dRatio')};
end

T = cell2table(rows, 'VariableNames', headers);
end


function v = get_event_field(ev, fieldName)
v = NaN;
if isstruct(ev) && isfield(ev, fieldName) && ~isempty(ev.(fieldName))
    raw = ev.(fieldName);
    if isnumeric(raw) || islogical(raw)
        v = raw(1);
    else
        v = raw;
    end
end
end


function reVal = get_case_re(allBreakup, ci)
reVal = NaN;
if isfield(allBreakup, 'Re') && numel(allBreakup) >= ci && ...
        ~isempty(allBreakup(ci).Re) && isfinite(allBreakup(ci).Re)
    reVal = double(allBreakup(ci).Re);
end
end


function [token, usedTokens] = make_unique_case_token(caseName, reVal, ci, usedTokens)
if isfinite(reVal)
    base = sprintf('Re%g_%s', reVal, char(string(caseName)));
else
    base = char(string(caseName));
end
base = regexprep(base, '[^A-Za-z0-9_.-]+', '_');
base = regexprep(base, '^_+', '');
base = regexprep(base, '_+$', '');
if isempty(base)
    base = sprintf('Case%d', ci);
end

token = base;
suffixIdx = 1;
while any(strcmpi(usedTokens, string(token)))
    token = sprintf('%s_%d', base, suffixIdx);
    suffixIdx = suffixIdx + 1;
end
usedTokens(end+1,1) = string(token); %#ok<AGROW>
end


function outTbl = append_table_compat_local(outTbl, newRows)
if isempty(outTbl) || width(outTbl) == 0
    outTbl = newRows;
    return;
end

if isempty(newRows) || height(newRows) == 0
    return;
end

outTbl = [outTbl; newRows]; %#ok<AGROW>
end
