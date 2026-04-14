function T = write_collapse_analysis_csv(allCollapse, csvFile)
%WRITE_COLLAPSE_ANALYSIS_CSV  Per-case collapse summary to CSV.
%
%   One row per case with event counts, rates, and peak-area statistics.

T = collapse_analysis_to_table(allCollapse);
if isempty(T) || height(T) == 0
    warning('write_collapse_analysis_csv: no cases, nothing written.');
    return;
end

write_table_csv_compat(T, csvFile);
fprintf('Saved collapse analysis CSV: %s\n', csvFile);
end
