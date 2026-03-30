function write_collapse_analysis_csv(allCollapse, csvFile)
%WRITE_COLLAPSE_ANALYSIS_CSV  Save per-case collapse rate summary and dominant
%   FFT frequencies to a CSV file.  One row per case.
%
%   Columns: Case, Re, kD, dt_s, frameRate_Hz,
%            nTotal, nQualified, nTruncated, nROIRejected,
%            ratePerFrame, ratePerSec,
%            domFreq1_Hz, domFreq1_power, ..., domFreqN_Hz, domFreqN_power

nCases = numel(allCollapse.caseName);
if nCases == 0
    warning('write_collapse_analysis_csv: no cases, nothing written.');
    return;
end

% Maximum number of dominant frequencies reported across all cases
nDomMax = 0;
for ci = 1:nCases
    cd = allCollapse.data{ci};
    if ~isempty(cd) && isfield(cd,'domFreqs_Hz')
        nDomMax = max(nDomMax, numel(cd.domFreqs_Hz));
    end
end

% Build header
header = {'Case','Re','kD','dt_s','frameRate_Hz', ...
    'nTotal','nQualified','nTruncated','nROIRejected', ...
    'ratePerFrame','ratePerSec'};
for j = 1:nDomMax
    header{end+1} = sprintf('domFreq%d_Hz',  j); %#ok<AGROW>
    header{end+1} = sprintf('domFreq%d_power', j); %#ok<AGROW>
end

% Build rows
rows = cell(nCases, numel(header));
for ci = 1:nCases
    cd  = allCollapse.data{ci};
    dt  = allCollapse.dt(ci);
    fsHz = 1 / max(dt, eps);

    if isempty(cd)
        cd = struct('nTotal',0,'nQualified',0,'nTruncated',0,'nROIRejected',0, ...
            'ratePerFrame',NaN,'ratePerSec',NaN, ...
            'domFreqs_Hz',nan(0,1),'domFreqPowers',nan(0,1));
    end

    col = 1;
    rows{ci,col} = allCollapse.caseName{ci}; col=col+1;
    rows{ci,col} = allCollapse.Re(ci);       col=col+1;
    rows{ci,col} = allCollapse.kD(ci);       col=col+1;
    rows{ci,col} = dt;                        col=col+1;
    rows{ci,col} = fsHz;                      col=col+1;
    rows{ci,col} = cd.nTotal;                 col=col+1;
    rows{ci,col} = cd.nQualified;             col=col+1;
    rows{ci,col} = cd.nTruncated;             col=col+1;
    rows{ci,col} = cd.nROIRejected;           col=col+1;
    rows{ci,col} = cd.ratePerFrame;           col=col+1;
    rows{ci,col} = cd.ratePerSec;             col=col+1;

    domF = cd.domFreqs_Hz(:);
    domP = cd.domFreqPowers(:);
    for j = 1:nDomMax
        if j <= numel(domF)
            rows{ci,col} = domF(j); col=col+1;
            rows{ci,col} = domP(j); col=col+1;
        else
            rows{ci,col} = NaN; col=col+1;
            rows{ci,col} = NaN; col=col+1;
        end
    end
end

T = cell2table(rows, 'VariableNames', header);
write_table_csv_compat(T, csvFile);
fprintf('Saved collapse analysis CSV: %s\n', csvFile);
end
