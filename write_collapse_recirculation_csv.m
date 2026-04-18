function T = write_collapse_recirculation_csv(allCollapseRecirculation, csvFile)
%WRITE_COLLAPSE_RECIRCULATION_CSV  Write per-case collapse recirculation summary.

T = collapse_recirculation_to_table(allCollapseRecirculation);
if isempty(T) || height(T) == 0
    warning('write_collapse_recirculation_csv: no cases, nothing written.');
    return;
end

outDir = fileparts(csvFile);
if ~isempty(outDir) && ~isfolder(outDir)
    mkdir(outDir);
end

write_table_csv_compat(T, csvFile);
fprintf('Saved collapse recirculation CSV: %s\n', csvFile);
end
