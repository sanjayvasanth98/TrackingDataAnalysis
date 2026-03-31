function write_breakup_analysis_xlsx(allBreakup, xlsxFile)
%WRITE_BREAKUP_ANALYSIS_XLSX  Save breakup event data to an Excel file.
%   One sheet per case, with every child-parent pair as one row.
%
%   allBreakup: struct array with fields .caseName, .kD, .events

if isempty(allBreakup)
    warning('write_breakup_analysis_xlsx: no cases, nothing written.');
    return;
end

nCases = numel(allBreakup);

headers = {'Case','kD','TrackID','ParentSpotID','ParentFrame', ...
    'ParentArea_px2','ParentAR','nChildren', ...
    'ChildSpotID','ChildFrame','ChildArea_px2','gamma','dRatio'};

% Delete existing file so old sheets don't persist.
if isfile(xlsxFile), delete(xlsxFile); end

for ci = 1:nCases
    ev = allBreakup(ci).events;
    sheetName = char(allBreakup(ci).caseName);

    if isempty(ev)
        writecell(headers, xlsxFile, 'Sheet', sheetName, 'Range', 'A1');
        continue;
    end

    nEv  = numel(ev);
    rows = cell(nEv, numel(headers));
    for ei = 1:nEv
        rows(ei,:) = { ...
            char(allBreakup(ci).caseName), allBreakup(ci).kD, ...
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
