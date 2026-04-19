function csvFiles = write_breakup_analysis_xlsx(allBreakup, xlsxFile)
%WRITE_BREAKUP_ANALYSIS_XLSX  Backward-compatible wrapper for CSV export.
%
%   Older code called this function to create one Excel sheet per case.
%   ARC MATLAB may not have writetable/writecell available, so this wrapper
%   now writes CSV files instead:
%     <requested_stem>.csv
%     <requested_stem>_by_case/*.csv

warning('write_breakup_analysis_xlsx:CSVOnly', ...
    ['Excel export has been replaced with CSV export for ARC compatibility. ', ...
    'Writing CSV files using the requested path stem.']);
csvFiles = write_breakup_analysis_csv(allBreakup, xlsxFile);
end
