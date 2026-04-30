%% test_event_level_activation_logistic.m
% Event-level logistic activation model.
%
% Model target:
%   Y = 1 for activated strict-primary/left-moving tracks
%   Y = 0 for non-activated strict-primary/left-moving tracks
%
% Predictors are built from lagrangian_acceleration_by_case.mat. Older .mat
% files may not contain track-level d0/tau/velocity fields; the script
% automatically drops unavailable or constant predictors.
%
% Important interpretation note:
%   aWindowStar is peak trigger-window a* for activated tracks, and peak
%   random-window a* for non-activated tracks. This makes the predictor
%   comparable across Y=1 and Y=0 rows.

clear; clc;

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot);
addpath(fileparts(mfilename('fullpath')));

%% Paths
% Option 1: point this to your plot_data_mat folder from a completed run.
matDir = test_plotmat_location("lagrangian_acceleration_by_case.mat");

% Option 2: or point directly to lagrangian_acceleration_by_case.mat.
% If non-empty, this takes priority over matDir.
matFile = "";

outDir = fullfile(fileparts(mfilename('fullpath')), ...
    'test_outputs', 'EventActivationLogistic');
if ~isfolder(outDir), mkdir(outDir); end

%% Plot/model options
set(0, 'DefaultFigureVisible', 'on');
plotOpts = struct();
plotOpts.enableNormalTheme = true;
plotOpts.enablePosterTheme = false;
plotOpts.savePNG = true;
plotOpts.saveSVG = false;
plotOpts.themes = "normal";
plotOpts.keepFiguresOpen = true;

modelOpts = struct();
modelOpts.minFinitePerPredictor = 10;
modelOpts.minRowsPerClass = 5;
modelOpts.nTriggerBins = 8;
modelOpts.candidatePredictors = ["kD", "Re", "d0_um", "tau_star", ...
    "aWindowStar", "upstreamVelocity_m_s"];

%% Load data and build event table
[allLagAccel, dataLabel] = load_lagrangian_mat(matDir, matFile);
eventTable = build_event_activation_table(allLagAccel);
fprintf('Using %s Lagrangian acceleration data.\n', dataLabel);
fprintf('Built %d event-level rows.\n', height(eventTable));

if height(eventTable) == 0
    error('No event-level rows were built. Check that allLagAccel contains track-level data.');
end

write_table_csv_compat(eventTable, fullfile(outDir, 'event_activation_logistic_rows.csv'));
save(fullfile(outDir, 'event_activation_logistic_rows.mat'), 'eventTable', 'modelOpts');

%% Fit logistic model
[fitResult, fitTable] = fit_event_logistic(eventTable, modelOpts);
fprintf('\nLogistic model method: %s\n', fitResult.method);
fprintf('Rows used in model: %d activated, %d non-activated\n', ...
    sum(fitTable.Y == 1), sum(fitTable.Y == 0));
fprintf('Predictors used:\n');
disp(fitResult.predictorNames(:));
disp(fitResult.coefTable);

write_table_csv_compat(fitResult.coefTable, fullfile(outDir, 'event_activation_logistic_coefficients.csv'));
save(fullfile(outDir, 'event_activation_logistic_model.mat'), 'fitResult', 'fitTable', 'modelOpts');

%% Plots
plot_logistic_odds_ratios(fitResult.coefTable, outDir, plotOpts);
plot_activation_probability_by_kd(fitTable, fitResult.predictedProbability, outDir, plotOpts);
plot_activation_probability_by_trigger(fitTable, fitResult.predictedProbability, modelOpts, outDir, plotOpts);

fprintf('Done. Output in: %s\n', outDir);


% =========================================================================
function [allLagAccel, dataLabel] = load_lagrangian_mat(matDir, matFile)
matDir = string(matDir);
matFile = string(matFile);

if strlength(strtrim(matFile)) == 0
    matFile = fullfile(matDir, "lagrangian_acceleration_by_case.mat");
end

if ~isfile(matFile)
    error(['Could not find Lagrangian acceleration MAT file:\n  %s\n' ...
        'Set matDir to a folder containing lagrangian_acceleration_by_case.mat, ' ...
        'or set matFile directly.'], matFile);
end

S = load(matFile);
if ~isfield(S, 'allLagAccel')
    error('MAT file does not contain allLagAccel: %s', matFile);
end
allLagAccel = S.allLagAccel;
dataLabel = "saved .mat";
fprintf('Loaded: %s\n', matFile);
end


% =========================================================================
function T = build_event_activation_table(allLagAccel)
caseName = strings(0, 1);
Re = nan(0, 1);
kD = nan(0, 1);
trackId = nan(0, 1);
Y = nan(0, 1);
aWindowStar = nan(0, 1);
aWindowType = strings(0, 1);
nWindowSamples = nan(0, 1);
d0_um = nan(0, 1);
d0Source = strings(0, 1);
tau_s = nan(0, 1);
tau_star = nan(0, 1);
upstreamVelocity_m_s = nan(0, 1);
upstreamVelocitySource = strings(0, 1);
growthRatio = nan(0, 1);

for ci = 1:numel(allLagAccel)
    d = allLagAccel(ci);
    if ~isfield(d, 'tracks') || isempty(d.tracks)
        continue;
    end

    caseD0_m = scalar_field(d, 'dMean_m', NaN);
    caseU_m_s = scalar_field(d, 'U_ref_m_s', NaN);
    tracks = d.tracks(:);

    for ti = 1:numel(tracks)
        tr = tracks(ti);
        isAct = logical(scalar_field(tr, 'isActivated', false));
        [aWin, winType, nWin] = activation_window_astar(tr, isAct);

        if isAct
            yVal = 1;
        else
            yVal = 0;
        end

        [d0Val_m, d0Src] = track_or_case_scalar(tr, 'd0_m', caseD0_m, "case_dMean");
        [uVal, uSrc] = track_or_case_scalar(tr, 'upstreamVelocity_m_s', caseU_m_s, "case_Uref");

        caseName(end+1, 1) = string(scalar_string_field(d, 'caseName', sprintf('case_%d', ci))); %#ok<AGROW>
        Re(end+1, 1) = scalar_field(d, 'Re', NaN); %#ok<AGROW>
        kD(end+1, 1) = scalar_field(d, 'kD', NaN); %#ok<AGROW>
        trackId(end+1, 1) = scalar_field(tr, 'TRACK_ID', NaN); %#ok<AGROW>
        Y(end+1, 1) = yVal; %#ok<AGROW>
        aWindowStar(end+1, 1) = aWin; %#ok<AGROW>
        aWindowType(end+1, 1) = winType; %#ok<AGROW>
        nWindowSamples(end+1, 1) = nWin; %#ok<AGROW>
        d0_um(end+1, 1) = d0Val_m * 1e6; %#ok<AGROW>
        d0Source(end+1, 1) = d0Src; %#ok<AGROW>
        tau_s(end+1, 1) = scalar_field(tr, 'tau_s', NaN); %#ok<AGROW>
        tau_star(end+1, 1) = scalar_field(tr, 'tau_star', NaN); %#ok<AGROW>
        upstreamVelocity_m_s(end+1, 1) = uVal; %#ok<AGROW>
        upstreamVelocitySource(end+1, 1) = uSrc; %#ok<AGROW>
        growthRatio(end+1, 1) = scalar_field(tr, 'growthRatio', NaN); %#ok<AGROW>
    end
end

T = table(caseName, Re, kD, trackId, Y, aWindowStar, aWindowType, ...
    nWindowSamples, d0_um, d0Source, tau_s, tau_star, ...
    upstreamVelocity_m_s, upstreamVelocitySource, growthRatio);
end


% =========================================================================
function [aWin, winType, nWin] = activation_window_astar(tr, isAct)
aWin = NaN;
nWin = 0;
if isAct
    winType = "trigger";
    aWin = scalar_field(tr, 'peakTriggerAstar', NaN);
    if ~isfinite(aWin)
        aWin = max_finite(vector_field(tr, 'triggerAstar'));
    end
    nWin = scalar_field(tr, 'nTriggerSamples', numel(vector_field(tr, 'triggerAstar')));
else
    winType = "random_control";
    aWin = max_finite(vector_field(tr, 'randomWindowAstar'));
    nWin = scalar_field(tr, 'nRandomWindowSamples', numel(vector_field(tr, 'randomWindowAstar')));
end
end


% =========================================================================
function [fitResult, Tfit] = fit_event_logistic(eventTable, opts)
predictorNames = select_available_predictors(eventTable, opts.candidatePredictors, opts);
if isempty(predictorNames)
    error('No usable predictors were available for the logistic model.');
end

valid = isfinite(eventTable.Y) & (eventTable.Y == 0 | eventTable.Y == 1);
for i = 1:numel(predictorNames)
    valid = valid & isfinite(eventTable.(char(predictorNames(i))));
end

Tfit = eventTable(valid, :);
if sum(Tfit.Y == 1) < opts.minRowsPerClass || sum(Tfit.Y == 0) < opts.minRowsPerClass
    error('Not enough rows per class after predictor filtering. Activated=%d, non-activated=%d.', ...
        sum(Tfit.Y == 1), sum(Tfit.Y == 0));
end

[predictorNames, Tfit, standardization] = standardize_predictors(Tfit, predictorNames);
if isempty(predictorNames)
    error('All predictors became constant after row filtering.');
end

zNames = standardization.zName;
formula = "Y ~ " + strjoin(zNames, " + ");

usedFitglm = exist('fitglm', 'file') == 2;
if usedFitglm
    try
        mdl = fitglm(Tfit, char(formula), 'Distribution', 'binomial', 'Link', 'logit');
        pred = predict(mdl, Tfit);
        coefStats = mdl.Coefficients;
        coefNames = string(coefStats.Properties.RowNames);
        beta = coefStats.Estimate;
        se = coefStats.SE;
        zStat = beta ./ se;
        pVal = coefStats.pValue;
        coefTable = make_coef_table(coefNames, beta, se, zStat, pVal, standardization);

        fitResult = struct();
        fitResult.method = "fitglm";
        fitResult.formula = formula;
        fitResult.model = mdl;
        fitResult.predictorNames = predictorNames;
        fitResult.standardization = standardization;
        fitResult.coefTable = coefTable;
        fitResult.predictedProbability = pred;
        return;
    catch ME
        warning('fitglm failed (%s). Falling back to local ridge-IRLS logistic fit.', ME.message);
    end
end

[beta, se, zStat, pVal, pred] = fit_logistic_irls(Tfit, zNames);
coefNames = ["Intercept"; zNames(:)];
coefTable = make_coef_table(coefNames, beta, se, zStat, pVal, standardization);

fitResult = struct();
fitResult.method = "local ridge-IRLS";
fitResult.formula = formula;
fitResult.model = [];
fitResult.predictorNames = predictorNames;
fitResult.standardization = standardization;
fitResult.coefTable = coefTable;
fitResult.predictedProbability = pred;
end


% =========================================================================
function predictorNames = select_available_predictors(T, candidatePredictors, opts)
predictorNames = strings(0, 1);
for i = 1:numel(candidatePredictors)
    name = string(candidatePredictors(i));
    if ~ismember(char(name), T.Properties.VariableNames)
        continue;
    end
    x = T.(char(name));
    x = x(isfinite(x));
    if numel(x) < opts.minFinitePerPredictor
        fprintf('Skipping predictor %s: only %d finite rows.\n', name, numel(x));
        continue;
    end
    if finite_std(x) <= 0
        fprintf('Skipping predictor %s: no variation.\n', name);
        continue;
    end
    predictorNames(end+1, 1) = name; %#ok<AGROW>
end
end


% =========================================================================
function [predictorNames, Tfit, standardization] = standardize_predictors(Tfit, predictorNames)
keep = true(numel(predictorNames), 1);
mu = nan(numel(predictorNames), 1);
sigma = nan(numel(predictorNames), 1);
zName = strings(numel(predictorNames), 1);

for i = 1:numel(predictorNames)
    x = Tfit.(char(predictorNames(i)));
    mu(i) = mean(x);
    sigma(i) = std(x, 0);
    if ~(isfinite(sigma(i)) && sigma(i) > 0)
        keep(i) = false;
        continue;
    end
    zName(i) = "z_" + predictorNames(i);
    Tfit.(char(zName(i))) = (x - mu(i)) ./ sigma(i);
end

predictorNames = predictorNames(keep);
mu = mu(keep);
sigma = sigma(keep);
zName = zName(keep);
standardization = table(predictorNames(:), zName(:), mu(:), sigma(:), ...
    'VariableNames', {'predictorName', 'zName', 'mean', 'std'});
end


% =========================================================================
function [beta, se, zStat, pVal, pred] = fit_logistic_irls(Tfit, zNames)
y = double(Tfit.Y(:));
n = numel(y);
p = numel(zNames) + 1;
X = ones(n, p);
for j = 1:numel(zNames)
    X(:, j + 1) = Tfit.(char(zNames(j)));
end

beta = zeros(p, 1);
ridge = 1e-6;
ridgeMat = ridge * eye(p);
ridgeMat(1, 1) = 0;

for iter = 1:100
    eta = X * beta;
    pred = logistic_link(eta);
    w = max(pred .* (1 - pred), eps);
    z = eta + (y - pred) ./ w;
    Xw = X .* w;
    betaNew = (X' * Xw + ridgeMat) \ (X' * (w .* z));
    if max(abs(betaNew - beta)) < 1e-8
        beta = betaNew;
        break;
    end
    beta = betaNew;
end

pred = logistic_link(X * beta);
w = max(pred .* (1 - pred), eps);
covBeta = pinv(X' * (X .* w) + ridgeMat);
se = sqrt(max(diag(covBeta), 0));
zStat = beta ./ se;
pVal = erfc(abs(zStat) ./ sqrt(2));
end


% =========================================================================
function coefTable = make_coef_table(coefNames, beta, se, zStat, pVal, standardization)
coefNames = string(coefNames(:));
predictorName = strings(numel(coefNames), 1);
predictorName(1) = "Intercept";
for i = 2:numel(coefNames)
    idx = find(standardization.zName == coefNames(i), 1, 'first');
    if isempty(idx)
        predictorName(i) = coefNames(i);
    else
        predictorName(i) = standardization.predictorName(idx);
    end
end

oddsRatio = exp(beta);
ciLow = exp(beta - 1.96 .* se);
ciHigh = exp(beta + 1.96 .* se);
displayName = arrayfun(@display_predictor_name, predictorName, 'UniformOutput', false);
displayName = string(displayName(:));

coefTable = table(coefNames, predictorName, displayName, beta(:), se(:), ...
    zStat(:), pVal(:), oddsRatio(:), ciLow(:), ciHigh(:), ...
    'VariableNames', {'Coefficient', 'Predictor', 'DisplayName', 'Beta', ...
    'SE', 'zStat', 'pValue', 'OddsRatio_per1SD', ...
    'OddsRatio_CI_low', 'OddsRatio_CI_high'});
end


% =========================================================================
function plot_logistic_odds_ratios(coefTable, outDir, plotOpts)
mask = coefTable.Predictor ~= "Intercept" & ...
    isfinite(coefTable.OddsRatio_per1SD) & ...
    isfinite(coefTable.OddsRatio_CI_low) & ...
    isfinite(coefTable.OddsRatio_CI_high);
C = coefTable(mask, :);
if isempty(C)
    warning('No finite odds ratios to plot.');
    return;
end

for theme = reshape(plotOpts.themes, 1, [])
    fontName = resolve_plot_font_name();
    f = figure('Color', 'w', 'Position', [120 120 920 620]);
    ax = axes(f);
    hold(ax, 'on');

    y = (1:height(C)).';
    for i = 1:height(C)
        plot(ax, [C.OddsRatio_CI_low(i), C.OddsRatio_CI_high(i)], [y(i), y(i)], ...
            '-', 'Color', [0.20 0.20 0.20], 'LineWidth', 1.4);
    end
    scatter(ax, C.OddsRatio_per1SD, y, 70, [0.05 0.38 0.78], ...
        'filled', 'MarkerEdgeColor', [0 0 0], 'LineWidth', 0.8);

    yl = [0.5, height(C) + 0.5];
    plot(ax, [1 1], yl, '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.0);
    ylim(ax, yl);
    set(ax, 'YTick', y, 'YTickLabel', cellstr(C.DisplayName), ...
        'YDir', 'reverse', 'XScale', 'log', 'FontName', fontName);
    xlabel(ax, 'Odds ratio per 1 SD increase', 'Interpreter', 'latex');
    ylabel(ax, '');
    title(ax, 'Event-level activation logistic model', 'FontName', fontName);
    grid(ax, 'off');
    box(ax, 'on');

    apply_plot_theme(ax, char(theme));
    save_fig_dual_safe(f, fullfile(outDir, "EventActivationLogistic_odds_ratios_" + theme), plotOpts);
    close_if_needed(f, plotOpts);
end
end


% =========================================================================
function plot_activation_probability_by_kd(Tfit, pred, outDir, plotOpts)
ReVals = unique(Tfit.Re(:));
ReVals = ReVals(isfinite(ReVals));
if isempty(ReVals)
    return;
end

for theme = reshape(plotOpts.themes, 1, [])
    fontName = resolve_plot_font_name();
    f = figure('Color', 'w', 'Position', [120 120 980 680]);
    ax = axes(f);
    hold(ax, 'on');

    cmap = lines(max(numel(ReVals), 1));
    markerStyles = {'o', 's', '^', 'd', 'v', '>', '<', 'p', 'h'};
    lgd = gobjects(0, 1);
    lgdTxt = strings(0, 1);

    for r = 1:numel(ReVals)
        Rei = ReVals(r);
        maskRe = Tfit.Re == Rei;
        kVals = unique(Tfit.kD(maskRe));
        kVals = sort(kVals(isfinite(kVals)));
        if isempty(kVals)
            continue;
        end

        obs = nan(numel(kVals), 1);
        pr = nan(numel(kVals), 1);
        lo = nan(numel(kVals), 1);
        hi = nan(numel(kVals), 1);
        n = nan(numel(kVals), 1);
        for i = 1:numel(kVals)
            m = maskRe & Tfit.kD == kVals(i);
            n(i) = sum(m);
            obs(i) = mean(Tfit.Y(m));
            pr(i) = mean(pred(m));
            [lo(i), hi(i)] = wilson_ci(sum(Tfit.Y(m) == 1), n(i));
        end

        col = cmap(r, :);
        markerStyle = markerStyles{mod(r - 1, numel(markerStyles)) + 1};
        hObs = errorbar(ax, kVals, obs, obs - lo, hi - obs, markerStyle, ...
            'LineStyle', '-', 'Color', col, 'LineWidth', 1.2, ...
            'CapSize', 8, 'MarkerSize', 9, 'MarkerFaceColor', col, ...
            'MarkerEdgeColor', [0 0 0]);
        plot(ax, kVals, pr, '--', 'Color', col, 'LineWidth', 2.0, ...
            'HandleVisibility', 'off');
        lgd(end+1,1) = hObs; %#ok<AGROW>
        lgdTxt(end+1,1) = sprintf('Observed + model, Re=%g', Rei); %#ok<AGROW>
    end

    xlabel(ax, '$k/d$', 'Interpreter', 'latex');
    ylabel(ax, 'Activation probability', 'Interpreter', 'latex');
    ylim(ax, [0 1]);
    title(ax, 'Observed and logistic-model activation probability', 'FontName', fontName);
    grid(ax, 'off');
    box(ax, 'on');
    if ~isempty(lgd)
        leg = legend(ax, lgd, cellstr(lgdTxt), 'Location', 'northeast', 'Box', 'off');
    else
        leg = [];
    end
    apply_plot_theme(ax, char(theme));
    style_legend_for_theme(leg, char(theme));
    save_fig_dual_safe(f, fullfile(outDir, "EventActivationLogistic_probability_vs_kD_" + theme), plotOpts);
    close_if_needed(f, plotOpts);
end
end


% =========================================================================
function plot_activation_probability_by_trigger(Tfit, pred, opts, outDir, plotOpts)
if ~ismember('aWindowStar', Tfit.Properties.VariableNames)
    return;
end
valid = isfinite(Tfit.aWindowStar) & Tfit.aWindowStar > 0 & isfinite(pred);
if sum(valid) < 20 || finite_std(Tfit.aWindowStar(valid)) <= 0
    return;
end

x = Tfit.aWindowStar(valid);
y = Tfit.Y(valid);
p = pred(valid);
[xBin, yObs, yPred, lo, hi, n] = equal_count_bins(x, y, p, opts.nTriggerBins); %#ok<ASGLU>
if isempty(xBin)
    return;
end

for theme = reshape(plotOpts.themes, 1, [])
    fontName = resolve_plot_font_name();
    f = figure('Color', 'w', 'Position', [120 120 920 650]);
    ax = axes(f);
    hold(ax, 'on');

    errorbar(ax, xBin, yObs, yObs - lo, hi - yObs, 'o-', ...
        'Color', [0.05 0.34 0.72], 'MarkerFaceColor', [0.05 0.34 0.72], ...
        'MarkerEdgeColor', [0 0 0], 'LineWidth', 1.4, 'CapSize', 8, ...
        'MarkerSize', 8);
    plot(ax, xBin, yPred, 's--', ...
        'Color', [0.83 0.22 0.08], 'MarkerFaceColor', [0.83 0.22 0.08], ...
        'MarkerEdgeColor', [0 0 0], 'LineWidth', 1.8, 'MarkerSize', 7);

    set(ax, 'XScale', 'log', 'FontName', fontName);
    xlabel(ax, 'Peak trigger/control-window $|a^*|$', 'Interpreter', 'latex');
    ylabel(ax, 'Activation probability', 'Interpreter', 'latex');
    ylim(ax, [0 1]);
    title(ax, 'Activation probability vs window acceleration', 'FontName', fontName);
    grid(ax, 'off');
    box(ax, 'on');
    leg = legend(ax, {'Observed bins', 'Logistic model'}, 'Location', 'northeast', 'Box', 'off');
    apply_plot_theme(ax, char(theme));
    style_legend_for_theme(leg, char(theme));
    save_fig_dual_safe(f, fullfile(outDir, "EventActivationLogistic_probability_vs_aWindowStar_" + theme), plotOpts);
    close_if_needed(f, plotOpts);
end
end


% =========================================================================
function [xBin, yObs, yPred, lo, hi, n] = equal_count_bins(x, y, pred, nBins)
[xSorted, ord] = sort(x(:));
ySorted = y(ord);
pSorted = pred(ord);
nTotal = numel(xSorted);
nBins = min(max(1, round(nBins)), nTotal);
edges = unique(round(linspace(1, nTotal + 1, nBins + 1)));
if numel(edges) < 3
    xBin = []; yObs = []; yPred = []; lo = []; hi = []; n = [];
    return;
end

xBin = nan(numel(edges) - 1, 1);
yObs = nan(size(xBin));
yPred = nan(size(xBin));
lo = nan(size(xBin));
hi = nan(size(xBin));
n = nan(size(xBin));
for b = 1:(numel(edges) - 1)
    idx = edges(b):(edges(b + 1) - 1);
    idx = idx(idx >= 1 & idx <= nTotal);
    if isempty(idx)
        continue;
    end
    n(b) = numel(idx);
    xBin(b) = median(xSorted(idx));
    yObs(b) = mean(ySorted(idx));
    yPred(b) = mean(pSorted(idx));
    [lo(b), hi(b)] = wilson_ci(sum(ySorted(idx) == 1), n(b));
end
valid = isfinite(xBin) & isfinite(yObs) & isfinite(yPred);
xBin = xBin(valid);
yObs = yObs(valid);
yPred = yPred(valid);
lo = lo(valid);
hi = hi(valid);
n = n(valid);
end


% =========================================================================
function [lo, hi] = wilson_ci(k, n)
if ~(isfinite(k) && isfinite(n) && n > 0)
    lo = NaN;
    hi = NaN;
    return;
end
z = 1.96;
p = k / n;
den = 1 + z^2 / n;
center = (p + z^2 / (2 * n)) / den;
half = z * sqrt((p * (1 - p) + z^2 / (4 * n)) / n) / den;
lo = max(0, center - half);
hi = min(1, center + half);
end


% =========================================================================
function y = logistic_link(eta)
eta = max(min(eta, 40), -40);
y = 1 ./ (1 + exp(-eta));
end


% =========================================================================
function [val, source] = track_or_case_scalar(s, fieldName, fallback, fallbackSource)
val = scalar_field(s, fieldName, NaN);
if isfinite(val)
    source = "track_" + string(fieldName);
else
    val = fallback;
    source = fallbackSource;
end
end


% =========================================================================
function val = scalar_field(s, fieldName, defaultVal)
val = defaultVal;
if isstruct(s) && isfield(s, fieldName) && ~isempty(s.(fieldName))
    raw = s.(fieldName);
    if isnumeric(raw) || islogical(raw)
        raw = double(raw(:));
        raw = raw(isfinite(raw));
        if ~isempty(raw)
            val = raw(1);
        end
    end
end
end


% =========================================================================
function val = scalar_string_field(s, fieldName, defaultVal)
val = string(defaultVal);
if isstruct(s) && isfield(s, fieldName) && ~isempty(s.(fieldName))
    raw = string(s.(fieldName));
    if ~isempty(raw)
        val = raw(1);
    end
end
end


% =========================================================================
function v = vector_field(s, fieldName)
v = nan(0, 1);
if isstruct(s) && isfield(s, fieldName) && ~isempty(s.(fieldName)) && isnumeric(s.(fieldName))
    v = double(s.(fieldName)(:));
end
end


% =========================================================================
function m = max_finite(x)
x = x(isfinite(x));
if isempty(x)
    m = NaN;
else
    m = max(x);
end
end


% =========================================================================
function s = finite_std(x)
x = x(isfinite(x));
if numel(x) < 2
    s = 0;
else
    s = std(x, 0);
end
end


% =========================================================================
function label = display_predictor_name(name)
name = string(name);
switch char(name)
    case 'kD'
        label = 'k/d';
    case 'Re'
        label = 'Re';
    case 'd0_um'
        label = 'd0 (um)';
    case 'tau_star'
        label = 'tau*';
    case 'aWindowStar'
        label = 'window |a*|';
    case 'upstreamVelocity_m_s'
        label = 'upstream velocity';
    case 'Intercept'
        label = 'Intercept';
    otherwise
        label = char(name);
end
end


% =========================================================================
function style_legend_for_theme(leg, theme)
if isempty(leg) || ~isgraphics(leg)
    return;
end
if strcmp(theme, 'poster')
    leg.TextColor = [1 1 1];
    leg.Color = [0 0 0];
else
    leg.TextColor = [0 0 0];
    leg.Color = 'none';
end
end


% =========================================================================
function close_if_needed(f, plotOpts)
if ~isfield(plotOpts, 'keepFiguresOpen') || ~plotOpts.keepFiguresOpen
    close(f);
end
end
