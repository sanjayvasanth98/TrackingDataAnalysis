function write_summary_csv_legacy(summaryRows, outCsv)
fid = fopen(outCsv, 'w');
if fid < 0
    error('Cannot open for writing: %s', outCsv);
end

fprintf(fid, 'Case,Sa,nTracksTotal,nInjected,nActivated,A_over_I,A/I_absolute,Tau_median\n');

n = numel(summaryRows);
for i = 1:n
    aOverIAbsolute = legacy_numeric_field(summaryRows(i), 'A_over_I_absolute', summaryRows(i).A_over_I);
    tauMedian = legacy_numeric_field(summaryRows(i), 'Tau_median', ...
        legacy_numeric_field(summaryRows(i), 'tau_median', NaN));
    fprintf(fid, '%s,%.15g,%d,%d,%d,%.15g,%.15g,%.15g\n', ...
        summaryRows(i).Case, ...
        summaryRows(i).Sa, ...
        summaryRows(i).nTracksTotal, ...
        summaryRows(i).nInjected, ...
        summaryRows(i).nActivated, ...
        summaryRows(i).A_over_I, ...
        aOverIAbsolute, ...
        tauMedian);
end

fclose(fid);
end

function val = legacy_numeric_field(s, fieldName, defaultVal)
if isfield(s, fieldName) && ~isempty(s.(fieldName)) && isscalar(s.(fieldName)) && isfinite(s.(fieldName))
    val = s.(fieldName);
else
    val = defaultVal;
end
end
