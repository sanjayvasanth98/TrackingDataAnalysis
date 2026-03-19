function save_fig_dual(figHandle, outBase, opts)
%SAVE_FIG_DUAL Backward-compatible wrapper around save_fig_dual_safe.
%
% Older plotting functions call save_fig_dual(fig, outBase). This wrapper
% keeps those calls working while centralizing save behavior.

if nargin < 3
    save_fig_dual_safe(figHandle, outBase);
else
    save_fig_dual_safe(figHandle, outBase, opts);
end
end
