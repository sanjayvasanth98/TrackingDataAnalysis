function write_collapse_analysis_csv(allCollapse, csvFile)
%WRITE_COLLAPSE_ANALYSIS_CSV  Per-case collapse summary to CSV.
%
%   One row per case with event counts, rates, and peak-area statistics.

nCases = numel(allCollapse.caseName);
if nCases == 0
    warning('write_collapse_analysis_csv: no cases, nothing written.');
    return;
end

hasPixelSize = isfield(allCollapse, 'pixelSize') && numel(allCollapse.pixelSize) >= nCases;

headers = {'Case','Re','kD','dt_s','frameRate_Hz', ...
    'nTotal','nQualified','nTruncated','nROIRejected', ...
    'ratePerFrame','ratePerSec', ...
    'peakArea_mean_px2','peakArea_median_px2','peakArea_std_px2'};
if hasPixelSize
    headers = [headers, {'peakDiam_mean_um','peakDiam_median_um'}];
end

nCols = numel(headers);
rows  = cell(nCases, nCols);

for ci = 1:nCases
    cd = allCollapse.data{ci};

    if isempty(cd)
        rows(ci,:) = [{allCollapse.caseName{ci}}, ...
            num2cell([allCollapse.Re(ci), allCollapse.kD(ci)]), ...
            repmat({NaN}, 1, nCols - 3)];
        continue;
    end

    dt_s = allCollapse.dt(ci);
    frameRate = 1 / max(dt_s, eps);

    % Peak-area statistics
    pa = cd.peakArea_px2(:);
    pa = pa(isfinite(pa) & pa > 0);
    if isempty(pa)
        paMean = NaN; paMedian = NaN; paStd = NaN;
        dMean  = NaN; dMedian  = NaN;
    else
        paMean   = mean(pa);
        paMedian = median(pa);
        paStd    = std(pa);
        if hasPixelSize
            d_eq    = 2 * sqrt(pa / pi) * allCollapse.pixelSize(ci) * 1000; % um
            dMean   = mean(d_eq);
            dMedian = median(d_eq);
        else
            dMean  = NaN;
            dMedian = NaN;
        end
    end

    baseRow = { ...
        allCollapse.caseName{ci}, ...
        allCollapse.Re(ci), ...
        allCollapse.kD(ci), ...
        dt_s, ...
        frameRate, ...
        cd.nTotal, ...
        cd.nQualified, ...
        cd.nTruncated, ...
        cd.nROIRejected, ...
        cd.ratePerFrame, ...
        cd.ratePerSec, ...
        paMean, paMedian, paStd};

    if hasPixelSize
        baseRow = [baseRow, {dMean, dMedian}]; %#ok<AGROW>
    end

    rows(ci,:) = baseRow;
end

T = cell2table(rows, 'VariableNames', headers);
write_table_csv_compat(T, csvFile);
fprintf('Saved collapse analysis CSV: %s\n', csvFile);
end
