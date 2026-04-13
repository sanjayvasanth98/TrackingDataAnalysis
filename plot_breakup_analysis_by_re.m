function plot_breakup_analysis_by_re(allBreakup, figDir, plotOpts, matDir, arThresholds)
%PLOT_BREAKUP_ANALYSIS_BY_RE  Save breakup plots split by Reynolds number.
%
%   The main comparison plots are handled elsewhere. Breakup plots are
%   intentionally written as individual Re-specific plot sets when multiple
%   Reynolds numbers are present, so cases with identical k/d do not pool
%   across Re.

if nargin < 3
    plotOpts = struct();
end
if nargin < 4
    matDir = "";
end
if nargin < 5 || isempty(arThresholds)
    arThresholds = [];
end
hasMatDir = strlength(string(matDir)) > 0;

if isempty(allBreakup)
    warning('plot_breakup_analysis_by_re: no breakup data. Skipping.');
    return;
end

if ~isfolder(figDir), mkdir(figDir); end

[groups, groupLabels] = breakup_groups_by_re(allBreakup);
nGroups = numel(groups);

for gi = 1:nGroups
    breakupNow = allBreakup(groups{gi});
    if isempty(breakupNow)
        continue;
    end

    if nGroups > 1
        outFigDir = fullfile(figDir, groupLabels(gi));
        if hasMatDir
            outMatDir = fullfile(matDir, "BreakupAnalysis", groupLabels(gi));
        else
            outMatDir = "";
        end
    else
        outFigDir = figDir;
        outMatDir = matDir;
    end

    if ~isfolder(outFigDir), mkdir(outFigDir); end
    if strlength(string(outMatDir)) > 0 && ~isfolder(outMatDir), mkdir(outMatDir); end

    plot_breakup_gamma_vs_dratio(breakupNow, outFigDir, plotOpts);
    plot_breakup_gamma_beeswarm_vs_kd(breakupNow, outFigDir, plotOpts, "", outMatDir);

    for arThr = reshape(arThresholds, 1, [])
        arTag = sprintf('AR%s', strrep(sprintf('%.1f', arThr), '.', 'p'));
        filteredBreakup = filter_breakup_by_ar(breakupNow, arThr); %#ok<NASGU>
        plot_breakup_gamma_beeswarm_vs_kd(filteredBreakup, outFigDir, plotOpts, arTag, outMatDir);
        if strlength(string(outMatDir)) > 0
            save(fullfile(outMatDir, sprintf("breakup_analysis_by_case_%s.mat", arTag)), 'filteredBreakup');
        end
    end

    plot_breakup_gamma_scatter_vs_ar(breakupNow, outFigDir, plotOpts, "", outMatDir);
    if nGroups > 1
        fprintf('  Breakup plots for %s saved to: %s\n', groupLabels(gi), outFigDir);
    end
end
end


function [groups, groupLabels] = breakup_groups_by_re(allBreakup)
nCases = numel(allBreakup);
groups = {1:nCases};
groupLabels = "all";

if ~isfield(allBreakup, 'Re')
    return;
end

reVals = nan(nCases, 1);
for ci = 1:nCases
    if ~isempty(allBreakup(ci).Re) && isfinite(allBreakup(ci).Re)
        reVals(ci) = double(allBreakup(ci).Re);
    end
end

uniqueRe = unique(reVals(isfinite(reVals)));
if numel(uniqueRe) <= 1
    return;
end

groups = cell(numel(uniqueRe), 1);
groupLabels = strings(numel(uniqueRe), 1);
for ri = 1:numel(uniqueRe)
    groups{ri} = find(reVals == uniqueRe(ri));
    groupLabels(ri) = sanitize_re_label(uniqueRe(ri));
end
end


function label = sanitize_re_label(ReVal)
raw = sprintf('Re_%g', ReVal);
label = string(regexprep(raw, '[^\w.-]', '_'));
end
