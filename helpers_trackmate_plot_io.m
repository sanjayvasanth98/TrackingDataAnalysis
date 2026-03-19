function helpers_trackmate_plot_io()
%HELPERS_TRACKMATE_PLOT_IO Backward-compatible no-op helper entrypoint.
%
% This repository now exposes plotting I/O helpers as standalone functions:
%   - apply_plot_theme.m
%   - save_fig_dual_safe.m
%   - save_fig_dual.m (compatibility wrapper)
%
% Keeping this function avoids breaking older scripts that call
% `helpers_trackmate_plot_io` for initialization.
end
