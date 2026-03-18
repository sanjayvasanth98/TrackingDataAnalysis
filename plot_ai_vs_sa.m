function plot_ai_vs_sa(summaryTable, figDir, fitTxtFile)

Sa = summaryTable.Sa(:);
y  = summaryTable.A_over_I(:);

% Remove zeros or negatives for log plotting/fit
valid = isfinite(Sa) & isfinite(y) & (y > 0);
Sa = Sa(valid);
y  = y(valid);

% Fit: log10(y) = a + b*Sa
Y = log10(y);
p = polyfit(Sa, Y, 1);
b = p(1);
a = p(2);

Yhat = polyval(p, Sa);
SSres = sum((Y - Yhat).^2);
SStot = sum((Y - mean(Y)).^2);
R2 = 1 - SSres/max(SStot, eps);

% For plotting fitted line
SaFit = linspace(min(Sa), max(Sa), 200).';
yFit  = 10.^(a + b*SaFit);

% Save fit info
fid = fopen(fitTxtFile, 'w');
fprintf(fid, "Fit (semilogy): log10(A/I) = a + b*Sa\n");
fprintf(fid, "a = %.6g\n", a);
fprintf(fid, "b = %.6g\n", b);
fprintf(fid, "R2 = %.6f\n\n", R2);
fprintf(fid, "Equivalent: A/I = 10^a * 10^(b*Sa)\n");
fprintf(fid, "         = %.6g * 10^(%.6g*Sa)\n", 10^a, b);
fclose(fid);

% Make plots in two themes
for theme = ["normal","poster"]
    f = figure('Color','w','Position',[100 100 900 650]);
    ax = axes(f); hold(ax,'on');

    semilogy(ax, Sa, y, 'o', 'LineWidth', 1.8, 'MarkerSize', 8);
    semilogy(ax, SaFit, yFit, '-', 'LineWidth', 2.2);

    grid(ax,'on'); box(ax,'on');
    xlabel(ax, 'Sa (roughness)');
    ylabel(ax, 'A/I (Activation / Injection)');
    title(ax, 'Activation/Injection vs Roughness');
    legend(ax, {'Cases','Fit'}, 'Location','southoutside', 'Box','off');

    apply_plot_theme(ax, theme);

    outBase = fullfile(figDir, "AI_vs_Sa_" + theme);
    save_fig_dual(f, outBase);

    close(f);
end

fprintf("Saved fit info: %s\n", fitTxtFile);
end