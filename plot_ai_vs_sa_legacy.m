function plot_ai_vs_sa_legacy(summaryRows, figDir, fitTxtFile)
% summaryRows: struct array with fields Sa, A_over_I

Sa = [summaryRows.Sa].';
y  = [summaryRows.A_over_I].';

valid = isfinite(Sa) & isfinite(y) & (y > 0);
Sa = Sa(valid);
y  = y(valid);

Y = log10(y);
p = polyfit(Sa, Y, 1);   % Y = b*Sa + a
b = p(1);
a = p(2);

Yhat = polyval(p, Sa);
SSres = sum((Y - Yhat).^2);
SStot = sum((Y - mean(Y)).^2);
R2 = 1 - SSres/max(SStot, eps);

SaFit = linspace(min(Sa), max(Sa), 200).';
yFit  = 10.^(a + b*SaFit);

% Save fit info text
fid = fopen(fitTxtFile, 'w');
if fid < 0
    warning('Could not write fit file: %s', fitTxtFile);
else
    fprintf(fid, 'Fit (semilogy): log10(A/I) = a + b*Sa\n');
    fprintf(fid, 'a = %.6g\n', a);
    fprintf(fid, 'b = %.6g\n', b);
    fprintf(fid, 'R2 = %.6f\n\n', R2);
    fprintf(fid, 'Equivalent: A/I = 10^a * 10^(b*Sa)\n');
    fprintf(fid, '         = %.6g * 10^(%.6g*Sa)\n', 10^a, b);
    fclose(fid);
end

themes = {'normal','poster'};
for ti = 1:numel(themes)
    theme = themes{ti};

    f = figure('Color','w','Position',[100 100 900 650]);
    ax = axes('Parent',f); hold(ax,'on');

    semilogy(ax, Sa, y, 'o', 'LineWidth', 1.8, 'MarkerSize', 8);
    semilogy(ax, SaFit, yFit, '-', 'LineWidth', 2.2);

    grid(ax,'on'); box(ax,'on');
    xlabel(ax, 'Sa (roughness)');
    ylabel(ax, 'A/I (Activation / Injection)');
    title(ax, 'Activation/Injection vs Roughness');

    % Legend (keep simple for compatibility)
    legend(ax, {'Cases','Fit'}, 'Location','southoutside');

    apply_plot_theme(ax, theme);

    outBase = fullfile(figDir, ['AI_vs_Sa_' theme]);
    save_fig_dual_safe(f, outBase);

    close(f);
end
end