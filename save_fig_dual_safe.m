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
if ~isempty(outDir) && ~exist(outDir,'dir')
    mkdir(outDir);
end

prepare_figure_for_export(figHandle);

if opts.savePNG
    % PNG 600 dpi
    try
        exportgraphics(figHandle, pngPath, 'Resolution', 600);
    catch ME
        try
            print(figHandle, pngPath, '-dpng', '-r600', '-opengl');
        catch ME2
            warning('PNG save failed for %s: %s | fallback: %s', outBase, ME.message, ME2.message);
        end
    end
end

if opts.saveSVG
    % SVG: try exportgraphics, then print, else saveas
    try
        exportgraphics(figHandle, svgPath, 'ContentType', 'vector');
    catch ME
        try
            print(figHandle, svgPath, '-dsvg');
        catch ME2
            try
                saveas(figHandle, svgPath);
            catch ME3
                warning('SVG save failed for %s: %s | fallback: %s | saveas: %s', ...
                    outBase, ME.message, ME2.message, ME3.message);
            end
        end
    end
end
end
