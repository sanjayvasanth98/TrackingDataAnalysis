function prepare_figure_for_export(figHandle)
%PREPARE_FIGURE_FOR_EXPORT  Normalize fonts and force layout before capture/export.

if nargin < 1 || isempty(figHandle) || ~isgraphics(figHandle)
    return;
end

if ~isgraphics(figHandle, 'figure')
    figHandle = ancestor(figHandle, 'figure');
end
if isempty(figHandle) || ~isgraphics(figHandle, 'figure')
    return;
end

fontName = resolve_plot_font_name();
set(figHandle, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');

axHandles = findall(figHandle, 'Type', 'axes');
for i = 1:numel(axHandles)
    ax = axHandles(i);
    if isprop(ax, 'FontName')
        ax.FontName = fontName;
    end
    local_set_font(ax.Title, fontName);
    local_set_font(ax.XLabel, fontName);
    local_set_font(ax.YLabel, fontName);
    if isprop(ax, 'ZLabel')
        local_set_font(ax.ZLabel, fontName);
    end
end

legendHandles = findall(figHandle, 'Type', 'Legend');
for i = 1:numel(legendHandles)
    leg = legendHandles(i);
    if isprop(leg, 'FontName')
        leg.FontName = fontName;
    end
    if isprop(leg, 'AutoUpdate')
        leg.AutoUpdate = 'off';
    end
end

colorbarHandles = findall(figHandle, 'Type', 'ColorBar');
for i = 1:numel(colorbarHandles)
    cb = colorbarHandles(i);
    if isprop(cb, 'FontName')
        cb.FontName = fontName;
    end
    if isprop(cb, 'Label')
        local_set_font(cb.Label, fontName);
    end
end

textHandles = findall(figHandle, 'Type', 'Text');
for i = 1:numel(textHandles)
    local_set_font(textHandles(i), fontName);
end

drawnow;
drawnow;
end

function local_set_font(h, fontName)
if isempty(h) || ~isgraphics(h)
    return;
end
if isprop(h, 'FontName')
    h.FontName = fontName;
end
end
