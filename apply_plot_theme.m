function apply_plot_theme(ax, theme)
%APPLY_PLOT_THEME  'normal' or 'poster' theme (poster = black bg + white axes)

if nargin < 2 || isempty(theme)
    theme = 'normal';
end
theme = lower(char(string(theme)));

fontName = resolve_plot_font_name();
set(ax, 'LineWidth', 1.2, 'FontSize', 14, 'FontName', fontName);
fig = ancestor(ax, 'figure');

if strcmp(theme,'normal')
    set(fig,'Color','w');
    set(ax,'Color','w', 'XColor','k','YColor','k');
    ax.GridColor = [0.7 0.7 0.7];
    ax.MinorGridColor = [0.85 0.85 0.85];
    textColor = [0 0 0];
else
    set(fig,'Color','k');
    set(ax,'Color','k', 'XColor','w','YColor','w');
    ax.GridColor = [1 1 1]*0.35;
    ax.MinorGridColor = [1 1 1]*0.22;
    textColor = [1 1 1];
end

if isgraphics(ax.Title)
    ax.Title.Color = textColor;
    ax.Title.FontName = fontName;
end
if isgraphics(ax.XLabel)
    ax.XLabel.Color = textColor;
    ax.XLabel.FontName = fontName;
end
if isgraphics(ax.YLabel)
    ax.YLabel.Color = textColor;
    ax.YLabel.FontName = fontName;
end

set(ax,'GridAlpha',0.25,'MinorGridAlpha',0.12);
box(ax,'on');
end
