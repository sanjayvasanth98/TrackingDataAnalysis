# TrackingDataAnalysis Edit Context (2026-03-25)

## Scope
This note captures the current analysis logic and the major code edits applied in this working branch so the pipeline can be resumed or tuned quickly.

## Core Workflow Decisions
- Bulk flow direction is treated as `left_to_right`.
- Strict recirculation is the authoritative primary set for metrics and blue tracks.
- A secondary microbubble-start rescue set is included for visualization and activation accounting.
- Inception plotting uses strict activation events, plus microbubble rescue activation events.

## Current Key Parameters (Local + ARC)
- `flowOpts.bulkDirection = "left_to_right"`
- `flowOpts.minNetDxCounterflow_mm = 0.08`
- `flowOpts.minNegativeStepFraction = 0.65`
- `flowOpts.maxPositiveStepFraction = 0.30`
- `flowOpts.requireRightOrigin = true`
- `flowOpts.rightOriginFrac = 0.25`
- `flowOpts.netLeftMinConsecutiveNegSteps = 2`
- `flowOpts.netLeftBandRescueEnabled = true`
- `flowOpts.netLeftBandRescueX_mm = [0.5 1.5]`
- `qcOpts.minTrackSpots = 2`
- `qcOpts.maxTrackGaps = 1`
- `qcOpts.rejectSplitMergeComplex = false`
- `activationOpts.areaJumpFactor = 6.0`
- `activationOpts.requiredPostFrames = 3`
- `activationOpts.microbubbleStartAreaRange_px2 = [1 120]`
- `activationOpts.microbubbleRequireOutsideStrictPrimary = true`

## New Origin Exclusion Gate
- Any track whose start point is inside:
- `x in [0, 0.5] mm`
- `y in [0, 1.2] mm`
- is excluded from consideration.
- This exclusion happens early in `trackmate_case_metrics` and therefore affects strict, legacy left-moving, and microbubble rescue eligibility.
- Diagnostics selection (`resolve_diagnostic_track_indices`) also drops excluded-origin tracks so they do not appear in GIF/static outputs.
- Console strict-gate summary now reports `originBox` rejects (`nRejectedOriginWindow`).

## Strict (Blue) vs Microbubble (Green) Semantics
- Blue tracks: strict primary track IDs (`strictPrimaryTrackIds`).
- Green tracks: non-strict tracks in the microbubble-start range (`1–120 px^2`).
- Activation events include strict + microbubble rescue events.
- Static diagnostics legend uses:
- `Left moving`
- `Microbubble start 1-120`
- `Activation events`

## Dynamic Color Switching for Green Tracks in GIFs
- Implemented in:
- `save_diagnostic_track_gifs.m`
- `save_video_overlay_gif_from_avi.m`
- Rule for microbubble (green) tracks:
- If last 3 consecutive step displacements are leftward, switch rendering to blue.
- If last 3 consecutive step displacements are rightward, switch back to green.
- This is a visual state machine per track using `recent_motion_signal(xTail, 3)`.
- Strict tracks remain blue at all times.

## Inception / Activation Location Selection
- In `choose_inception_activation_xy` (local and ARC):
- Primary source is `strictActivationEvent_xy`.
- Then microbubble activation points are appended (`microbubbleActivationEvent_xy_nonLeft`).
- Unique rows are retained.

## Caching and Reproducibility
- Cache keys include policy tags from parser/QC/flow/activation settings.
- The new origin-exclusion policy is included in `build_cache_policy_tag`.
- Changing these knobs creates a distinct cache key and avoids stale-policy reuse.

## Main Files Touched
- `main_batch_trackmate_local.m`
- `main_batch_trackmate_arc.m`
- `trackmate_case_metrics.m`
- `save_diagnostic_track_gifs.m`
- `save_video_overlay_gif_from_avi.m`
- `plot_verification_tracks_for_case.m`
- `resolve_diagnostic_track_indices.m`

## Operational Notes
- Area thresholds for microbubble logic are in `px^2` (not converted to physical area).
- MATLAB runtime was not available in this shell environment, so changes were verified by code inspection and grep-level consistency checks, not by execution.
