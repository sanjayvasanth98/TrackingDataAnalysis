%% covergence_checkcode.m
% Inspect the convergence summary CSV file, check convergence by case, and
% save diagnostic plots/tables for the main metrics.
%
% Output folder:
%   convergence_check_outputs
%
% Main settings to tune:
%   relativeTolerancePercent  - maximum allowed step-to-step percent change
%   consecutiveStepsRequired  - number of consecutive increments that must
%                               satisfy the tolerance to call a case converged

clear; close all; clc;

%% User settings
relativeTolerancePercent = 5;      % Percent change tolerance for convergence
consecutiveStepsRequired = 3;      % Number of final increments that must be stable
makeFiguresVisible = false;        % false saves figures without opening windows
makeFullRangePlots = false;        % true also saves all-chunk plots in the main output folder

inputPath = 'E:\March Re 90,000 inception data\Processed images\results\results 35\convergence_summary_all_cases.csv';
outputFolderName = 'convergence_check_outputs6';

% Optional extra plot window for only a subset of chunks.
% Examples:
%   []       -> do not make extra chunk-window plots
%   [6 Inf]  -> chunks 6 through the last chunk in each case
%   [3 8]    -> chunks 3 through 8
% Note: these extra plots zoom/filter the CSV's cumulative values. They
% do not recalculate metrics as if earlier chunks were removed.
chunkRangeForExtraPlots = [1 Inf];

% Metrics emphasized in the convergence decision. The code only uses columns
% that are present in the CSV.
primaryMetricCandidates = {
    'A_over_I'
    'leftMovingActivated %'
    'AE_leftMoving %'
    'leftMovingFrac_total'
    'activationFrac_valid'
    };

% Metrics plotted for inspection. The code skips any missing columns.
plotMetricCandidates = {
    'A_over_I'
    'leftMovingActivated %'
    'AE_leftMoving %'
    'leftMovingFrac_total'
    'activationFrac_valid'
    'tau_mean'
    'activatedVelocity_mean_m_s'
    'nValidTracks'
    'nLeftMovingTracks'
    'nLeftMovingActivated'
    'leftMovingFrameExposure'
    };

%% Locate inputs and create output folder
scriptPath = mfilename('fullpath');
if isempty(scriptPath)
    rootDir = pwd;
else
    rootDir = fileparts(scriptPath);
end

outDir = fullfile(rootDir, outputFolderName);

if exist(inputPath, 'file') ~= 2
    error('Input CSV file not found: %s', inputPath);
end
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

fprintf('\nReading convergence CSV:\n  %s\n', inputPath);

%% Read the CSV while preserving original headers when possible
try
    opts = detectImportOptions(inputPath, 'VariableNamingRule', 'preserve');
catch
    opts = detectImportOptions(inputPath);
end
T = readtable(inputPath, opts);

if isempty(T)
    error('The CSV was read, but the table is empty.');
end

varNames = T.Properties.VariableNames;
fprintf('Loaded %d rows and %d columns.\n', height(T), width(T));

%% Identify core columns
caseVar = pickVariable(varNames, {'Case'});
frameVar = pickVariable(varNames, {'CumulativeFramesUsed', 'CumulativeFrameTarget', 'FrameBin'});
frameBinVar = pickVariable(varNames, {'FrameBin'});
reVar = pickVariable(varNames, {'Re'});
kdVar = pickVariable(varNames, {'kD'});

if isempty(caseVar)
    error('Could not find a Case column in the CSV.');
end
if isempty(frameVar)
    error('Could not find a frame column. Expected CumulativeFramesUsed, CumulativeFrameTarget, or FrameBin.');
end

caseValues = asStringColumn(T.(caseVar));
frames = asNumericColumn(T.(frameVar));
if isempty(frameBinVar)
    chunkBins = nan(height(T), 1);
else
    chunkBins = asNumericColumn(T.(frameBinVar));
end

if isempty(reVar)
    reValues = nan(height(T), 1);
else
    reValues = asNumericColumn(T.(reVar));
end
if isempty(kdVar)
    kdValues = nan(height(T), 1);
else
    kdValues = asNumericColumn(T.(kdVar));
end

% Stable grouping label used throughout the script.
groupLabelsAllRows = caseValues + " | Re=" + compactNumber(reValues) + " | kD=" + compactNumber(kdValues);
[originalGroupLabels, originalGroupIndex] = stableStringGroups(groupLabelsAllRows);

% Sort rows so each case keeps the CSV order, with frame order inside it.
[~, sortIdx] = sortrows(table(originalGroupIndex, frames));
T = T(sortIdx, :);
caseValues = caseValues(sortIdx);
frames = frames(sortIdx);
chunkBins = chunkBins(sortIdx);
reValues = reValues(sortIdx);
kdValues = kdValues(sortIdx);
groupLabelsAllRows = groupLabelsAllRows(sortIdx);

groupLabels = originalGroupLabels;
groupIndex = originalGroupIndex(sortIdx);
nGroups = numel(groupLabels);

fprintf('Detected %d case group(s).\n', nGroups);
for g = 1:nGroups
    rows = find(groupIndex == g);
    fprintf('  %-32s  rows=%d, frames=%g to %g\n', char(groupLabels(g)), numel(rows), min(frames(rows)), max(frames(rows)));
end

primaryMetricVars = keepExistingVariables(varNames, primaryMetricCandidates);
plotMetricVars = keepExistingVariables(varNames, plotMetricCandidates);

if isempty(primaryMetricVars)
    error('None of the requested primary convergence metrics were found in the CSV.');
end

fprintf('\nPrimary convergence metrics:\n');
fprintf('  %s\n', primaryMetricVars{:});
fprintf('\nConvergence rule: abs(step percent change) <= %.3g%% for the last %d increment(s).\n', ...
    relativeTolerancePercent, consecutiveStepsRequired);

%% Check convergence by group and metric
summaryRows = table();
metricRows = table();

for g = 1:nGroups
    rows = find(groupIndex == g);
    thisFrames = frames(rows);
    thisCase = caseValues(rows(1));
    thisRe = reValues(rows(1));
    thisKd = kdValues(rows(1));

    stableMatrix = false(numel(rows), numel(primaryMetricVars));
    primaryLastPct = nan(1, numel(primaryMetricVars));
    primaryLastWindowMax = nan(1, numel(primaryMetricVars));

    for m = 1:numel(primaryMetricVars)
        metricName = primaryMetricVars{m};
        y = asNumericColumn(T.(metricName));
        y = y(rows);
        pctChange = getStepPercentChange(T, rows, metricName);
        stableMatrix(:, m) = pctChange <= relativeTolerancePercent;
        primaryLastPct(m) = pctChange(end);
        if numel(rows) >= consecutiveStepsRequired
            primaryLastWindowMax(m) = max(pctChange(end-consecutiveStepsRequired+1:end), [], 'omitnan');
        end
    end

    allPrimaryStable = all(stableMatrix, 2);
    [isConvergedAtEnd, firstConvergedLocalIndex] = findConvergencePoint(allPrimaryStable, consecutiveStepsRequired);

    if isnan(firstConvergedLocalIndex)
        firstConvergedFrame = nan;
    else
        firstConvergedFrame = thisFrames(firstConvergedLocalIndex);
    end

    finalFrame = thisFrames(end);

    newSummary = table(thisCase, thisRe, thisKd, string(groupLabels(g)), numel(rows), finalFrame, ...
        isConvergedAtEnd, firstConvergedFrame, relativeTolerancePercent, consecutiveStepsRequired, ...
        'VariableNames', {'Case', 'Re', 'kD', 'GroupLabel', 'Rows', 'FinalFramesUsed', ...
        'ConvergedAtEnd', 'FirstConvergedFrame', 'TolerancePercent', 'ConsecutiveSteps'});
    summaryRows = [summaryRows; newSummary]; %#ok<AGROW>

    metricsForDiagnostics = unique([primaryMetricVars(:); plotMetricVars(:)], 'stable');
    for m = 1:numel(metricsForDiagnostics)
        metricName = metricsForDiagnostics{m};
        yAll = asNumericColumn(T.(metricName));
        y = yAll(rows);
        pctChange = getStepPercentChange(T, rows, metricName);

        if numel(rows) >= consecutiveStepsRequired
            lastWindowMax = max(pctChange(end-consecutiveStepsRequired+1:end), [], 'omitnan');
        else
            lastWindowMax = nan;
        end

        isPrimary = any(strcmp(metricName, primaryMetricVars));
        metricConvergedAtEnd = lastWindowMax <= relativeTolerancePercent;

        newMetricRow = table(thisCase, thisRe, thisKd, string(groupLabels(g)), string(metricName), isPrimary, ...
            y(end), pctChange(end), lastWindowMax, metricConvergedAtEnd, ...
            'VariableNames', {'Case', 'Re', 'kD', 'GroupLabel', 'Metric', 'UsedForDecision', ...
            'FinalValue', 'FinalStepPctChange', 'MaxPctChangeLastWindow', 'MetricConvergedAtEnd'});
        metricRows = [metricRows; newMetricRow]; %#ok<AGROW>
    end
end

%% Save convergence tables
summaryCsv = fullfile(outDir, 'convergence_status_by_case.csv');
metricCsv = fullfile(outDir, 'convergence_metric_diagnostics.csv');
writetable(summaryRows, summaryCsv);
writetable(metricRows, metricCsv);

summaryXlsx = fullfile(outDir, 'convergence_check_results.xlsx');
try
    writetable(summaryRows, summaryXlsx, 'Sheet', 'case_status');
    writetable(metricRows, summaryXlsx, 'Sheet', 'metric_diagnostics');
catch ME
    warning('Could not write XLSX summary (%s). CSV files were still written.', ME.message);
end

fprintf('\nSaved convergence tables:\n  %s\n  %s\n', summaryCsv, metricCsv);

%% Make plots
figVisibility = 'off';
if makeFiguresVisible
    figVisibility = 'on';
end

if makeFullRangePlots || isempty(chunkRangeForExtraPlots)
    plotImportantMetrics(T, groupIndex, groupLabels, frames, plotMetricVars, outDir, figVisibility);
    plotAOverIWithCI(T, groupIndex, groupLabels, frames, outDir, figVisibility);
    plotPrimaryPercentChanges(T, groupIndex, groupLabels, frames, primaryMetricVars, ...
        relativeTolerancePercent, outDir, figVisibility);
    plotFinalComparison(summaryRows, metricRows, outDir, figVisibility);
else
    fprintf('\nSkipped full-range plots because makeFullRangePlots = false. Existing full-range PNGs in the main output folder may be from an earlier run.\n');
end

if ~isempty(chunkRangeForExtraPlots)
    if isempty(frameBinVar)
        warning('chunkRangeForExtraPlots was requested, but no FrameBin column was found. Skipping chunk-window plots.');
    else
        [chunkMask, chunkRangeLabel] = chunkRangeMask(chunkBins, chunkRangeForExtraPlots);
        if ~any(chunkMask)
            warning('No rows matched chunkRangeForExtraPlots = [%g %g]. Skipping chunk-window plots.', ...
                chunkRangeForExtraPlots(1), chunkRangeForExtraPlots(2));
        else
            chunkOutDir = fullfile(outDir, chunkRangeLabel);
            if ~exist(chunkOutDir, 'dir')
                mkdir(chunkOutDir);
            end
            fprintf('\nChunk-window selection %s matched %d row(s). Saving filtered plots here:\n  %s\n', ...
                chunkRangeLabel, nnz(chunkMask), chunkOutDir);

            [Tchunk, groupIndexChunk, groupLabelsChunk, framesChunk, caseValuesChunk, reValuesChunk, kdValuesChunk] = ...
                subsetGroupedRows(T, groupIndex, groupLabels, frames, caseValues, reValues, kdValues, chunkMask);

            [summaryRowsChunk, metricRowsChunk] = buildConvergenceTables(Tchunk, groupIndexChunk, groupLabelsChunk, ...
                framesChunk, caseValuesChunk, reValuesChunk, kdValuesChunk, primaryMetricVars, plotMetricVars, ...
                relativeTolerancePercent, consecutiveStepsRequired);

            writetable(summaryRowsChunk, fullfile(chunkOutDir, 'convergence_status_by_case.csv'));
            writetable(metricRowsChunk, fullfile(chunkOutDir, 'convergence_metric_diagnostics.csv'));

            plotImportantMetrics(Tchunk, groupIndexChunk, groupLabelsChunk, framesChunk, plotMetricVars, chunkOutDir, figVisibility);
            plotAOverIWithCI(Tchunk, groupIndexChunk, groupLabelsChunk, framesChunk, chunkOutDir, figVisibility);
            plotPrimaryPercentChanges(Tchunk, groupIndexChunk, groupLabelsChunk, framesChunk, primaryMetricVars, ...
                relativeTolerancePercent, chunkOutDir, figVisibility);
            plotFinalComparison(summaryRowsChunk, metricRowsChunk, chunkOutDir, figVisibility);

            fprintf('\nSaved extra chunk-window plots/tables for %s to:\n  %s\n', ...
                chunkRangeLabel, chunkOutDir);
        end
    end
end

%% Print concise command-window summary
fprintf('\nConvergence summary by case:\n');
for i = 1:height(summaryRows)
    if summaryRows.ConvergedAtEnd(i)
        statusText = 'CONVERGED';
    else
        statusText = 'NOT converged';
    end
    fprintf('  %-32s  %-13s  final frames=%g', char(summaryRows.GroupLabel(i)), statusText, summaryRows.FinalFramesUsed(i));
    if ~isnan(summaryRows.FirstConvergedFrame(i))
        fprintf(', first converged at %g frames', summaryRows.FirstConvergedFrame(i));
    end
    fprintf('\n');
end

fprintf('\nMain output folder:\n  %s\n', outDir);
if ~isempty(chunkRangeForExtraPlots)
    fprintf('Chunk-window plots are saved in the chunks_* subfolder listed above.\n');
end
fprintf('Finished convergence check.\n\n');

%% Local functions
function [summaryRows, metricRows] = buildConvergenceTables(T, groupIndex, groupLabels, frames, ...
    caseValues, reValues, kdValues, primaryMetricVars, plotMetricVars, relativeTolerancePercent, consecutiveStepsRequired)
summaryRows = table();
metricRows = table();
metricsForDiagnostics = unique([primaryMetricVars(:); plotMetricVars(:)], 'stable');

for g = 1:numel(groupLabels)
    rows = find(groupIndex == g);
    if isempty(rows)
        continue;
    end

    thisFrames = frames(rows);
    thisCase = caseValues(rows(1));
    thisRe = reValues(rows(1));
    thisKd = kdValues(rows(1));

    stableMatrix = false(numel(rows), numel(primaryMetricVars));
    for m = 1:numel(primaryMetricVars)
        metricName = primaryMetricVars{m};
        pctChange = getStepPercentChange(T, rows, metricName);
        stableMatrix(:, m) = pctChange <= relativeTolerancePercent;
    end

    allPrimaryStable = all(stableMatrix, 2);
    [isConvergedAtEnd, firstConvergedLocalIndex] = findConvergencePoint(allPrimaryStable, consecutiveStepsRequired);

    if isnan(firstConvergedLocalIndex)
        firstConvergedFrame = nan;
    else
        firstConvergedFrame = thisFrames(firstConvergedLocalIndex);
    end

    newSummary = table(thisCase, thisRe, thisKd, string(groupLabels(g)), numel(rows), thisFrames(end), ...
        isConvergedAtEnd, firstConvergedFrame, relativeTolerancePercent, consecutiveStepsRequired, ...
        'VariableNames', {'Case', 'Re', 'kD', 'GroupLabel', 'Rows', 'FinalFramesUsed', ...
        'ConvergedAtEnd', 'FirstConvergedFrame', 'TolerancePercent', 'ConsecutiveSteps'});
    summaryRows = [summaryRows; newSummary]; %#ok<AGROW>

    for m = 1:numel(metricsForDiagnostics)
        metricName = metricsForDiagnostics{m};
        yAll = asNumericColumn(T.(metricName));
        y = yAll(rows);
        pctChange = getStepPercentChange(T, rows, metricName);

        if numel(rows) >= consecutiveStepsRequired
            lastWindowMax = max(pctChange(end-consecutiveStepsRequired+1:end), [], 'omitnan');
        else
            lastWindowMax = nan;
        end

        isPrimary = any(strcmp(metricName, primaryMetricVars));
        metricConvergedAtEnd = lastWindowMax <= relativeTolerancePercent;

        newMetricRow = table(thisCase, thisRe, thisKd, string(groupLabels(g)), string(metricName), isPrimary, ...
            y(end), pctChange(end), lastWindowMax, metricConvergedAtEnd, ...
            'VariableNames', {'Case', 'Re', 'kD', 'GroupLabel', 'Metric', 'UsedForDecision', ...
            'FinalValue', 'FinalStepPctChange', 'MaxPctChangeLastWindow', 'MetricConvergedAtEnd'});
        metricRows = [metricRows; newMetricRow]; %#ok<AGROW>
    end
end
end

function [mask, label] = chunkRangeMask(chunkBins, chunkRange)
if numel(chunkRange) ~= 2
    error('chunkRangeForExtraPlots must be empty or a two-value vector such as [6 Inf] or [3 8].');
end
chunkStart = chunkRange(1);
chunkEnd = chunkRange(2);
if chunkEnd < chunkStart
    error('chunkRangeForExtraPlots end value must be greater than or equal to the start value.');
end

mask = chunkBins >= chunkStart & chunkBins <= chunkEnd;
if isinf(chunkEnd)
    label = sprintf('chunks_%g_to_end', chunkStart);
else
    label = sprintf('chunks_%g_to_%g', chunkStart, chunkEnd);
end
label = regexprep(label, '[^A-Za-z0-9_\\-]+', '_');
end

function [Tsub, groupIndexSub, groupLabelsSub, framesSub, caseValuesSub, reValuesSub, kdValuesSub] = ...
    subsetGroupedRows(T, groupIndex, groupLabels, frames, caseValues, reValues, kdValues, mask)
Tsub = T(mask, :);
framesSub = frames(mask);
caseValuesSub = caseValues(mask);
reValuesSub = reValues(mask);
kdValuesSub = kdValues(mask);
oldGroupIndexSub = groupIndex(mask);

oldGroupIds = [];
for i = 1:numel(oldGroupIndexSub)
    if ~any(oldGroupIds == oldGroupIndexSub(i))
        oldGroupIds(end+1, 1) = oldGroupIndexSub(i); %#ok<AGROW>
    end
end

groupLabelsSub = strings(numel(oldGroupIds), 1);
groupIndexSub = zeros(numel(oldGroupIndexSub), 1);
for g = 1:numel(oldGroupIds)
    groupLabelsSub(g) = groupLabels(oldGroupIds(g));
    groupIndexSub(oldGroupIndexSub == oldGroupIds(g)) = g;
end
end

function selected = keepExistingVariables(varNames, candidates)
selected = {};
for i = 1:numel(candidates)
    v = pickVariable(varNames, candidates(i));
    if ~isempty(v)
        selected{end+1, 1} = v; %#ok<AGROW>
    end
end
selected = unique(selected, 'stable');
end

function [groupLabels, groupIndex] = stableStringGroups(values)
% Keep first-seen group order without relying on unique(..., 'stable') with
% three outputs, which can vary across MATLAB releases.
values = string(values(:));
groupLabels = strings(0, 1);
groupIndex = zeros(numel(values), 1);
for i = 1:numel(values)
    matchIdx = find(groupLabels == values(i), 1, 'first');
    if isempty(matchIdx)
        groupLabels(end+1, 1) = values(i); %#ok<AGROW>
        matchIdx = numel(groupLabels);
    end
    groupIndex(i) = matchIdx;
end
end

function varName = pickVariable(varNames, candidates)
varName = '';
if ischar(candidates) || isstring(candidates)
    candidates = cellstr(candidates);
end
normVars = normalizeName(varNames);
for i = 1:numel(candidates)
    target = normalizeName(candidates{i});
    idx = find(strcmp(normVars, target), 1, 'first');
    if ~isempty(idx)
        varName = varNames{idx};
        return;
    end
end
for i = 1:numel(candidates)
    target = normalizeName(candidates{i});
    idx = find(contains(normVars, target), 1, 'first');
    if ~isempty(idx)
        varName = varNames{idx};
        return;
    end
end
end

function out = normalizeName(in)
if iscell(in)
    out = cell(size(in));
    for k = 1:numel(in)
        out{k} = normalizeName(in{k});
    end
    return;
end
out = lower(char(in));
out = strrep(out, '%', 'pct');
out = regexprep(out, '[^a-z0-9]+', '');
out = strrep(out, 'percent', 'pct');
end

function x = asNumericColumn(v)
if isnumeric(v)
    x = double(v);
elseif islogical(v)
    x = double(v);
elseif iscell(v)
    x = str2double(string(v));
elseif iscategorical(v)
    x = str2double(string(v));
elseif isstring(v) || ischar(v)
    x = str2double(string(v));
else
    try
        x = double(v);
    catch
        x = str2double(string(v));
    end
end
x = x(:);
end

function s = asStringColumn(v)
if iscell(v)
    s = string(v);
elseif iscategorical(v)
    s = string(v);
elseif isstring(v)
    s = v;
elseif ischar(v)
    s = string(cellstr(v));
else
    s = string(v);
end
s = s(:);
end

function s = compactNumber(x)
s = strings(size(x));
for i = 1:numel(x)
    if isnan(x(i))
        s(i) = "NA";
    else
        s(i) = string(sprintf('%.6g', x(i)));
    end
end
end

function pctChange = getStepPercentChange(T, rows, metricName)
varNames = T.Properties.VariableNames;
pctVar = pctChangeVariableForMetric(varNames, metricName);
if ~isempty(pctVar)
    pctChange = abs(asNumericColumn(T.(pctVar)));
    pctChange = pctChange(rows);
else
    yAll = asNumericColumn(T.(metricName));
    y = yAll(rows);
    pctChange = nan(size(y));
    for i = 2:numel(y)
        denom = max(abs(y(i-1)), eps);
        pctChange(i) = abs((y(i) - y(i-1)) ./ denom) .* 100;
    end
end
pctChange(~isfinite(pctChange)) = nan;
end

function pctVar = pctChangeVariableForMetric(varNames, metricName)
pctVar = '';
metricNorm = normalizeName(metricName);

% Workbook convention: pctChange_A_over_I, pctChange_leftMovingActivated_pct, etc.
candidates = {
    ['pctChange_' metricName]
    ['pctChange_' strrep(metricName, ' %', '_pct')]
    ['pctChange_' strrep(metricName, '%', 'pct')]
    };
for i = 1:numel(candidates)
    v = pickVariable(varNames, candidates{i});
    if ~isempty(v)
        pctVar = v;
        return;
    end
end

normVars = normalizeName(varNames);
patterns = {
    ['pctchange' metricNorm]
    ['pctchange' strrep(metricNorm, 'pct', '') 'pct']
    };
for p = 1:numel(patterns)
    idx = find(strcmp(normVars, patterns{p}), 1, 'first');
    if ~isempty(idx)
        pctVar = varNames{idx};
        return;
    end
end
end

function [isConvergedAtEnd, firstConvergedIndex] = findConvergencePoint(allStable, consecutiveStepsRequired)
firstConvergedIndex = nan;
if numel(allStable) < consecutiveStepsRequired + 1
    isConvergedAtEnd = false;
    return;
end
for i = consecutiveStepsRequired + 1:numel(allStable)
    window = allStable(i-consecutiveStepsRequired+1:i);
    if all(window)
        firstConvergedIndex = i;
        break;
    end
end
lastWindow = allStable(end-consecutiveStepsRequired+1:end);
isConvergedAtEnd = all(lastWindow);
end

function plotImportantMetrics(T, groupIndex, groupLabels, frames, metricVars, outDir, figVisibility)
if isempty(metricVars)
    return;
end
nMetrics = numel(metricVars);
nCols = min(3, nMetrics);
nRows = ceil(nMetrics / nCols);
fig = figure('Color', 'w', 'Visible', figVisibility, 'Position', [80 80 1600 430*nRows]);
colors = lines(numel(groupLabels));
legendHandles = gobjects(numel(groupLabels), 1);

for m = 1:nMetrics
    ax = subplot(nRows, nCols, m); %#ok<LAXES>
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    yAll = asNumericColumn(T.(metricVars{m}));
    for g = 1:numel(groupLabels)
        rows = find(groupIndex == g);
        h = plot(ax, frames(rows), yAll(rows), '-o', 'LineWidth', 1.3, 'MarkerSize', 4, ...
            'Color', colors(g, :), 'DisplayName', char(groupLabels(g)));
        if m == 1
            legendHandles(g) = h;
        end
    end
    xlabel(ax, 'Cumulative frames used');
    ylabel(ax, prettifyMetricName(metricVars{m}));
    title(ax, prettifyMetricName(metricVars{m}), 'Interpreter', 'none');
    makeAxisReadable(ax);
end

% Put the legend in an unused tile when possible. Attaching a bestoutside
% legend to the final data axis can squeeze that last plot badly.
if nMetrics < nRows * nCols
    legendAx = subplot(nRows, nCols, nMetrics + 1); %#ok<LAXES>
    hold(legendAx, 'on');
    axis(legendAx, 'off');
    dummyLegendHandles = gobjects(numel(groupLabels), 1);
    for g = 1:numel(groupLabels)
        dummyLegendHandles(g) = plot(legendAx, nan, nan, '-o', 'LineWidth', 1.3, 'MarkerSize', 4, ...
            'Color', colors(g, :));
    end
    axis(legendAx, 'off');
    legend(legendAx, dummyLegendHandles, cellstr(groupLabels), 'Location', 'northwest', 'Interpreter', 'none');
else
    legend(legendHandles, cellstr(groupLabels), 'Location', 'bestoutside', 'Interpreter', 'none');
end
saveFigure(fig, fullfile(outDir, 'important_metrics_overview.png'));
close(fig);

% Also save one larger plot per metric for easier reporting.
for m = 1:nMetrics
    fig = figure('Color', 'w', 'Visible', figVisibility, 'Position', [120 120 1050 720]);
    ax = axes(fig); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    yAll = asNumericColumn(T.(metricVars{m}));
    for g = 1:numel(groupLabels)
        rows = find(groupIndex == g);
        plot(ax, frames(rows), yAll(rows), '-o', 'LineWidth', 1.6, 'MarkerSize', 5, ...
            'Color', colors(g, :), 'DisplayName', char(groupLabels(g)));
    end
    xlabel(ax, 'Cumulative frames used');
    ylabel(ax, prettifyMetricName(metricVars{m}));
    title(ax, ['Convergence trend: ' prettifyMetricName(metricVars{m})], 'Interpreter', 'none');
    legend(ax, 'Location', 'bestoutside', 'Interpreter', 'none');
    saveFigure(fig, fullfile(outDir, ['trend_' safeFileName(metricVars{m}) '.png']));
    close(fig);
end
end

function plotAOverIWithCI(T, groupIndex, groupLabels, frames, outDir, figVisibility)
varNames = T.Properties.VariableNames;
aoiVar = pickVariable(varNames, {'A_over_I'});
if isempty(aoiVar)
    return;
end
lowVar = pickVariable(varNames, {'A_over_I_err_low'});
highVar = pickVariable(varNames, {'A_over_I_err_high'});
ciLowVar = pickVariable(varNames, {'A_over_I_ci_low'});
ciHighVar = pickVariable(varNames, {'A_over_I_ci_high'});

fig = figure('Color', 'w', 'Visible', figVisibility, 'Position', [120 120 1100 760]);
ax = axes(fig); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
colors = lines(numel(groupLabels));
yAll = asNumericColumn(T.(aoiVar));

for g = 1:numel(groupLabels)
    rows = find(groupIndex == g);
    x = frames(rows);
    y = yAll(rows);
    if ~isempty(lowVar) && ~isempty(highVar)
        eLow = asNumericColumn(T.(lowVar));
        eHigh = asNumericColumn(T.(highVar));
        errorbar(ax, x, y, eLow(rows), eHigh(rows), '-o', 'LineWidth', 1.4, ...
            'MarkerSize', 4, 'Color', colors(g, :), 'DisplayName', char(groupLabels(g)));
    elseif ~isempty(ciLowVar) && ~isempty(ciHighVar)
        lo = asNumericColumn(T.(ciLowVar));
        hi = asNumericColumn(T.(ciHighVar));
        errorbar(ax, x, y, y-lo(rows), hi(rows)-y, '-o', 'LineWidth', 1.4, ...
            'MarkerSize', 4, 'Color', colors(g, :), 'DisplayName', char(groupLabels(g)));
    else
        plot(ax, x, y, '-o', 'LineWidth', 1.5, 'MarkerSize', 5, ...
            'Color', colors(g, :), 'DisplayName', char(groupLabels(g)));
    end
end
xlabel(ax, 'Cumulative frames used');
ylabel(ax, 'A over I');
title(ax, 'A over I convergence with confidence/error bounds', 'Interpreter', 'none');
legend(ax, 'Location', 'bestoutside', 'Interpreter', 'none');
saveFigure(fig, fullfile(outDir, 'A_over_I_with_CI.png'));
close(fig);
end

function plotPrimaryPercentChanges(T, groupIndex, groupLabels, frames, primaryMetricVars, tolerance, outDir, figVisibility)
nGroups = numel(groupLabels);
nCols = min(3, nGroups);
nRows = ceil(nGroups / nCols);
colors = lines(numel(primaryMetricVars));
fig = figure('Color', 'w', 'Visible', figVisibility, 'Position', [60 60 1500 400*nRows]);

for g = 1:nGroups
    ax = subplot(nRows, nCols, g); %#ok<LAXES>
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    rows = find(groupIndex == g);
    for m = 1:numel(primaryMetricVars)
        pctChange = getStepPercentChange(T, rows, primaryMetricVars{m});
        semilogy(ax, frames(rows), pctChange, '-o', 'LineWidth', 1.2, 'MarkerSize', 4, ...
            'Color', colors(m, :), 'DisplayName', prettifyMetricName(primaryMetricVars{m}));
    end
    xLimits = xlim(ax);
    semilogy(ax, xLimits, [tolerance tolerance], 'k--', 'LineWidth', 1.1, ...
        'DisplayName', sprintf('%.3g%% tolerance', tolerance));
    xlim(ax, xLimits);
    xlabel(ax, 'Cumulative frames used');
    ylabel(ax, 'Absolute step change (%)');
    title(ax, char(groupLabels(g)), 'Interpreter', 'none');
end
legend(gca, 'Location', 'bestoutside', 'Interpreter', 'none');
saveFigure(fig, fullfile(outDir, 'primary_metric_percent_changes.png'));
close(fig);
end

function plotFinalComparison(summaryRows, metricRows, outDir, figVisibility)
keyMetrics = {'A_over_I', 'leftMovingActivated %', 'AE_leftMoving %', 'tau_mean', 'activatedVelocity_mean_m_s'};
fig = figure('Color', 'w', 'Visible', figVisibility, 'Position', [120 120 1500 850]);

availableMetrics = strings(0, 1);
for i = 1:numel(keyMetrics)
    if any(metricRows.Metric == string(keyMetrics{i}))
        availableMetrics(end+1, 1) = string(keyMetrics{i}); %#ok<AGROW>
    end
end
if isempty(availableMetrics)
    close(fig);
    return;
end

nMetrics = numel(availableMetrics);
nCols = min(3, nMetrics);
nRows = ceil(nMetrics / nCols);

for m = 1:nMetrics
    ax = subplot(nRows, nCols, m); %#ok<LAXES>
    thisMetricRows = metricRows(metricRows.Metric == availableMetrics(m), :);
    [~, order] = ismember(summaryRows.GroupLabel, thisMetricRows.GroupLabel);
    order = order(order > 0);
    vals = thisMetricRows.FinalValue(order);
    bar(ax, vals);
    grid(ax, 'on'); box(ax, 'on');
    set(ax, 'XTick', 1:numel(order), 'XTickLabel', cellstr(thisMetricRows.Case(order)));
    xtickangle(ax, 35);
    ylabel(ax, ['Final ' prettifyMetricName(char(availableMetrics(m)))]);
    title(ax, ['Final ' prettifyMetricName(char(availableMetrics(m)))], 'Interpreter', 'none');
end
saveFigure(fig, fullfile(outDir, 'final_metric_comparison.png'));
close(fig);
end

function label = prettifyMetricName(metricName)
label = char(metricName);
label = strrep(label, '_', ' ');
label = strrep(label, 'pct', '%');
label = strtrim(label);
end

function makeAxisReadable(ax)
% Avoid title collisions with automatic x10^n text on dense overview plots.
try
    ax.YAxis.Exponent = 0;
catch
end
try
    ax.Title.FontSize = 11;
    ax.XLabel.FontSize = 10;
    ax.YLabel.FontSize = 10;
catch
end
end

function name = safeFileName(metricName)
name = char(metricName);
name = strrep(name, '%', 'pct');
name = regexprep(name, '[^A-Za-z0-9_\-]+', '_');
name = regexprep(name, '_+', '_');
name = regexprep(name, '^_|_$', '');
end

function saveFigure(fig, filePath)
try
    exportgraphics(fig, filePath, 'Resolution', 220);
catch
    saveas(fig, filePath);
end
end
