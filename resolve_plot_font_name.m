function fontName = resolve_plot_font_name()
%RESOLVE_PLOT_FONT_NAME  Pick a stable cross-platform plot font.

persistent cachedFontName

if ~isempty(cachedFontName)
    fontName = cachedFontName;
    return;
end

preferredFonts = [ ...
    "Times New Roman"
    "Liberation Serif"
    "Nimbus Roman"
    "TeX Gyre Termes"
    "DejaVu Serif"
    "Times"
    "Arial"
    "Helvetica"];

fontName = 'Helvetica';
availableFonts = strings(0,1);
try
    availableFonts = string(listfonts);
catch
    availableFonts = strings(0,1);
end

availableFonts = strip(availableFonts);
for i = 1:numel(preferredFonts)
    if any(strcmpi(preferredFonts(i), availableFonts))
        fontName = char(preferredFonts(i));
        break;
    end
end

cachedFontName = fontName;
end
