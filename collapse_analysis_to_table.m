function T = collapse_analysis_to_table(allCollapse)
%COLLAPSE_ANALYSIS_TO_TABLE  Convert per-case collapse summaries to a table.
%
%   T = collapse_analysis_to_table(allCollapse)
%
%   Shared table builder used by write_collapse_analysis_csv and by the
%   chunk-level aggregate CSVs.

nCases = numel(allCollapse.caseName);
if nCases == 0
    T = table();
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
end
