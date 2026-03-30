function plot_collapse_power_spectrum(allCollapse, figDir, plotOpts)
%PLOT_COLLAPSE_POWER_SPECTRUM  FFT power spectrum of per-frame collapse-count
%   signal for each case. All cases in one figure.  Y-axis: log scale.
%
%   The signal is the per-frame collapse-count vector computed by
%   analyze_collapse_events.  It is Hann-windowed and zero-meaned before
%   FFT so the plot shows the AC (periodic) content only.
%
%   Dominant peaks (top-N by power, pre-identified in analyze_collapse_events)
%   are marked with filled triangles on each case's spectrum line.
%
%   X-axis in kHz (frame rate ~102 kfps → Nyquist ~51 kHz).

if nargin < 3 || ~isfield(plotOpts,'themes') || isempty(plotOpts.themes)
    plotOpts.themes = "normal";
end

nCases = numel(allCollapse.caseName);
if nCases == 0
    warning('No collapse data for power spectrum plot.');
    return;
end

cmap = lines(max(nCases,1));

for theme = reshape(plotOpts.themes, 1, [])
    fontName = resolve_plot_font_name();
    f  = figure('Color','w','Position',[100 100 1100 650]);
    ax = axes(f);
    hold(ax,'on');

    lgd    = gobjects(0,1);
    lgdTxt = strings(0,1);
    anyPlotted = false;

    for ci = 1:nCases
        cd = allCollapse.data{ci};
        if isempty(cd) || isempty(cd.fftFreq_Hz) || cd.nQualified < 4
            continue;
        end

        freqKHz = cd.fftFreq_Hz(:) / 1000;
        power   = cd.fftPower(:);
        col     = cmap(ci,:);

        % Remove non-positive values to avoid log-scale issues
        validIdx = power > 0 & isfinite(power);

        % Skip DC bin (index 1, freq = 0)
        validIdx(1) = false;

        if ~any(validIdx), continue; end

        hLine = plot(ax, freqKHz(validIdx), power(validIdx), '-', ...
            'Color', col, 'LineWidth', 1.6);
        anyPlotted = true;

        lgd(end+1,1)    = hLine; %#ok<AGROW>
        lgdTxt(end+1,1) = sprintf('k/d=%.4g  (n=%d)', ...
            allCollapse.kD(ci), cd.nQualified); %#ok<AGROW>

        % Mark dominant peaks with filled triangles
        domF = cd.domFreqs_Hz(:);
        domP = cd.domFreqPowers(:);
        if ~isempty(domF)
            domFkHz = domF / 1000;
            validPeak = domP > 0 & isfinite(domP) & isfinite(domFkHz);
            if any(validPeak)
                plot(ax, domFkHz(validPeak), domP(validPeak), '^', ...
                    'Color', col, ...
                    'MarkerFaceColor', col, ...
                    'MarkerEdgeColor', [1 1 1], ...
                    'MarkerSize', 8, ...
                    'LineWidth', 0.8, ...
                    'HandleVisibility','off');
            end
        end
    end

    if ~anyPlotted
        warning('plot_collapse_power_spectrum: no cases with sufficient data to plot.');
        close(f);
        continue;
    end

    xlabel(ax, 'Frequency (kHz)', 'Interpreter','latex');
    ylabel(ax, 'Power (a.u.)', 'Interpreter','latex');
    title(ax,'');
    set(ax, 'YScale','log', 'FontName',fontName);
    grid(ax,'off');
    box(ax,'on');

    if ~isempty(lgd)
        leg = legend(ax, lgd, cellstr(lgdTxt), ...
            'Location','northoutside','NumColumns',3,'Box','off');
    else
        leg = [];
    end

    apply_plot_theme(ax, char(theme));
    style_legend_for_theme(leg, char(theme));
    set(ax, 'YScale','log');  % re-apply after theme

    outBase = fullfile(figDir, "CollapsePowerSpectrum_" + theme);
    save_fig_dual_safe(f, outBase, plotOpts);
    if ~isfield(plotOpts,'keepFiguresOpen') || ~plotOpts.keepFiguresOpen
        close(f);
    end
end
fprintf('Saved collapse power spectrum to: %s\n', figDir);
end


% =========================================================================
function style_legend_for_theme(leg, theme)
if isempty(leg) || ~isgraphics(leg), return; end
if strcmp(theme,'poster')
    leg.TextColor = [1 1 1];
    leg.Color     = [0 0 0];
else
    leg.TextColor = [0 0 0];
    leg.Color     = 'none';
end
end
