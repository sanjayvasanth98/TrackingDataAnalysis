function write_summary_csv_legacy(summaryRows, outCsv)
fid = fopen(outCsv, 'w');
if fid < 0
    error('Cannot open for writing: %s', outCsv);
end

fprintf(fid, 'Case,Sa,nTracksTotal,nInjected,nActivated,A_over_I\n');

n = numel(summaryRows);
for i = 1:n
    fprintf(fid, '%s,%.15g,%d,%d,%d,%.15g\n', ...
        summaryRows(i).Case, ...
        summaryRows(i).Sa, ...
        summaryRows(i).nTracksTotal, ...
        summaryRows(i).nInjected, ...
        summaryRows(i).nActivated, ...
        summaryRows(i).A_over_I);
end

fclose(fid);
end