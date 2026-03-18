function save_fig_dual_safe(figHandle, outBase, opts)
%SAVE_FIG_DUAL_SAFE  Save PNG (600 dpi) + SVG with fallback, never crash job.
%
% Usage:
%   save_fig_dual_safe(gcf, '/path/figure_name_no_ext')
%   save_fig_dual_safe(gcf, '/path/figure_name_no_ext', struct('savePNG', true, 'saveSVG', false))

if isstring(outBase)
    outBase = char(outBase);
end

if nargin < 3 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, 'savePNG'), opts.savePNG = true; end
if ~isfield(opts, 'saveSVG'), opts.saveSVG = true; end

pngPath = [outBase '.png'];
svgPath = [outBase '.svg'];

[outDir,~,~] = fileparts(outBase);
if ~exist(outDir,'dir'), mkdir(outDir); end

if opts.savePNG
    % PNG 600 dpi
    try
        print(figHandle, pngPath, '-dpng', '-r600');
    catch ME
        try
            exportgraphics(figHandle, pngPath, 'Resolution', 600);
        catch ME2
            warning('PNG save failed for %s: %s | fallback: %s', outBase, ME.message, ME2.message);
        end
    end
end

if opts.saveSVG
    % SVG: try print, else saveas
    try
        print(figHandle, svgPath, '-dsvg');
    catch
        try
            saveas(figHandle, svgPath);
        catch ME
            warning('SVG save failed for %s: %s', outBase, ME.message);
        end
    end
end
end
