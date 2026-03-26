function plot_ai_vs_kdh_re(summaryTable, figDir, fitTxtFile, plotOpts)

if nargin < 4 || ~isfield(plotOpts, 'themes') || isempty(plotOpts.themes)
    plotOpts.themes = "normal";
end

ReVals = unique(summaryTable.Re(:));
if isempty(ReVals)
    warning('No rows found in summary table. Skipping A/I vs k/d plot.');
    return;
end

% Write fit summary by Reynolds set (once, outside theme loop)
fid = fopen(fitTxtFile, 'w');
if fid < 0
    warning('Could not open fit text file: %s', fitTxtFile);
else
    fprintf(fid, 'Fit model by Re: log10(A/I) = a + b*(k/d)\n\n');
    for r = 1:numel(ReVals)
        Rei = ReVals(r);
        sub = summaryTable(summaryTable.Re == Rei, :);
        sub = sortrows(sub, 'kD');
        x = sub.kD(:);
        y = sub.A_over_I(:);
        valid = isfinite(x) & isfinite(y) & (y > 0);
        x = x(valid);
        y = y(valid);
        if numel(x) >= 2
            p = polyfit(x, log10(y), 1);
            b = p(1);
            a = p(2);
            Y = log10(y);
            Yhat = polyval(p, x);
            SSres = sum((Y - Yhat).^2);
            SStot = sum((Y - mean(Y)).^2);
            R2 = 1 - SSres / max(SStot, eps);
            fprintf(fid, 'Re = %g\n', Rei);
            fprintf(fid, 'a = %.8g\n', a);
            fprintf(fid, 'b = %.8g\n', b);
            fprintf(fid, 'R2 = %.6f\n', R2);
            fprintf(fid, 'A/I = %.8g * 10^(%.8g * (k/d))\n\n', 10^a, b);
        else
            fprintf(fid, 'Re = %g\n', Rei);
            fprintf(fid, 'Not enough points for fit (need >=2).\n\n');
        end
    end
    fclose(fid);
end

for theme = reshape(plotOpts.themes, 1, [])
    f = figure('Color', 'w', 'Position', [100 100 1000 700]);
    ax = axes(f);
    hold(ax, 'on');

    markerFaceColor = [0 0 1];
    markerEdgeColor = [0 0 0];
    errorBarColor = [0 0 0];
    fitColor = [0.35 0.35 0.35];
    lgd = gobjects(0,1);
    lgdTxt = strings(0,1);

    for r = 1:numel(ReVals)
        Rei = ReVals(r);
        sub = summaryTable(summaryTable.Re == Rei, :);
        sub = sortrows(sub, 'kD');

        x = sub.kD(:);
        y = sub.A_over_I(:);
        ciLow = sub.A_over_I_ci_low(:);
        ciHigh = sub.A_over_I_ci_high(:);

        valid = isfinite(x) & isfinite(y) & (y > 0);
        x = x(valid);
        y = y(valid);
        ciLow = ciLow(valid);
        ciHigh = ciHigh(valid);

        if isempty(x)
            continue;
        end

        if ~isempty(ciLow)
            ciLow(~isfinite(ciLow)) = y(~isfinite(ciLow));
            ciHigh(~isfinite(ciHigh)) = y(~isfinite(ciHigh));
        end

        tinyBound = max(realmin, y .* 1e-6);
        lowerBound = max(ciLow, tinyBound);
        upperBound = max(ciHigh, y);
        errLow = max(0, y - lowerBound);
        errHigh = max(0, upperBound - y);

        hErr = errorbar(ax, x, y, errLow, errHigh, 'o', ...
            'LineStyle', 'none', ...
            'LineWidth', 0.75, ...
            'MarkerSize', 12, ...
            'CapSize', 8, ...
            'Color', errorBarColor, ...
            'MarkerFaceColor', markerFaceColor, ...
            'MarkerEdgeColor', markerEdgeColor);
        lgd(end+1,1) = hErr; %#ok<AGROW>
        lgdTxt(end+1,1) = sprintf('Data, Re=%g', Rei); %#ok<AGROW>

        if numel(x) >= 2
            p = polyfit(x, log10(y), 1); % log10(y)=b*x+a
            b = p(1);
            a = p(2);

            xFit = linspace(min(x), max(x), 200).';
            yFit = 10.^(a + b*xFit);
            hFit = semilogy(ax, xFit, yFit, '--', 'LineWidth', 1.8, 'Color', fitColor);
            lgd(end+1,1) = hFit; %#ok<AGROW>
            lgdTxt(end+1,1) = sprintf('Fit, Re=%g', Rei); %#ok<AGROW>

        end
    end

    xlabel(ax, '$k/d$', 'Interpreter', 'latex');
    ylabel(ax, '$A/I$', 'Interpreter', 'latex');
    title(ax, 'Activation/Injection vs k/d');
    set(ax, 'YScale', 'log');
    grid(ax, 'off');
    box(ax, 'on');

    if ~isempty(lgd)
        leg = legend(ax, lgd, cellstr(lgdTxt), 'Location', 'southoutside', 'NumColumns', 2, 'Box', 'off');
    else
        leg = [];
    end

    apply_plot_theme(ax, char(theme));
    style_legend_for_theme(leg, char(theme));

    outBase = fullfile(figDir, "AI_vs_kD_by_Re_" + theme);
    save_fig_dual_safe(f, outBase, plotOpts);
    close(f);
end

fprintf('Saved fit info: %s\n', fitTxtFile);
end

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
