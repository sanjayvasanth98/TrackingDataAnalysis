function apply_plot_theme(ax, theme)
set(ax, 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman');

if strcmp(theme,'normal')
    set(gcf,'Color','w');
    set(ax,'Color','w', 'XColor','k','YColor','k');
    ax.GridColor = [0.7 0.7 0.7];
    ax.MinorGridColor = [0.85 0.85 0.85];
else
    set(gcf,'Color','k');
    set(ax,'Color','k', 'XColor','w','YColor','w');
    ax.GridColor = [1 1 1]*0.35;
    ax.MinorGridColor = [1 1 1]*0.22;
end

set(ax,'GridAlpha',0.25,'MinorGridAlpha',0.12);
grid(ax,'on');
box(ax,'on');
end

function save_fig_dual_safe(figHandle, outBase)
% PNG always; SVG try print then saveas; never crash job on save errors.

[outDir,~,~] = fileparts(outBase);
if ~exist(outDir,'dir'), mkdir(outDir); end

% Force software rendering on headless nodes (helps SVG/PNG reliability)
try
    opengl('save','software');
catch
end

% PNG 600 dpi
try
    print(figHandle, [outBase '.png'], '-dpng', '-r600');
catch ME
    warning('PNG save failed for %s: %s', outBase, ME.message);
end

% SVG
try
    print(figHandle, [outBase '.svg'], '-dsvg');
catch
    try
        saveas(figHandle, [outBase '.svg']);
    catch ME
        warning('SVG save failed for %s: %s', outBase, ME.message);
    end
end
end