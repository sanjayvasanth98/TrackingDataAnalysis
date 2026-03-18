function apply_plot_theme(ax, theme)
%APPLY_PLOT_THEME  'normal' or 'poster' theme (poster = black bg + white axes)

set(ax, 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Times New Roman');

if strcmp(theme,'normal')
    set(gcf,'Color','w');
    set(ax,'Color','w', 'XColor','k','YColor','k');
    ax.GridColor = [0.7 0.7 0.7];
    ax.MinorGridColor = [0.85 0.85 0.85];
    textColor = [0 0 0];
else
    set(gcf,'Color','k');
    set(ax,'Color','k', 'XColor','w','YColor','w');
    ax.GridColor = [1 1 1]*0.35;
    ax.MinorGridColor = [1 1 1]*0.22;
    textColor = [1 1 1];
end

ax.Title.Color = textColor;
ax.XLabel.Color = textColor;
ax.YLabel.Color = textColor;

set(ax,'GridAlpha',0.25,'MinorGridAlpha',0.12);
grid(ax,'on');
box(ax,'on');
end
