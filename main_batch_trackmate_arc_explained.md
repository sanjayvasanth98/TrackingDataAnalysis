# `main_batch_trackmate_arc.m` Explained

This document explains the current logic in:

- `TrackingDataAnalysis/main_batch_trackmate_arc.m`
- and the core functions it calls:
  - `analyze_trackmate_xml_arc.m`
  - `trackmate_case_metrics.m`
  - `plot_verification_tracks_for_case.m`
  - `save_diagnostic_track_gifs.m`
  - plotting/CSV helpers

## 1) What this script does

`main_batch_trackmate_arc.m` is the batch driver for your TrackMate XML workflow.

For each case (5um, 12um, ...), it:

1. Loads parser cache if available.
2. Parses XML tracks/spots (or reuses cached parsed output).
3. Computes unified recirculation/activation metrics.
4. Writes summary/debug CSV files.
5. Generates diagnostic and publication plots.

At the end it writes:

- `Summary_tracks.csv` (primary per-case physics table)
- `track_gate_summary.csv` (strict legacy gating audit)
- `framewise_counts_<case>_Re_<...>_kDh_<...>.csv` (per-frame denominator/numerator traces)
- figure outputs in `Figures_PNG_SVG/...`

## 2) High-level data flow

`main_batch_trackmate_arc` -> `analyze_trackmate_xml_arc` -> `trackmate_case_metrics` -> CSV + plots + GIF diagnostics

The key consistency design is:

- One shared track/event catalog is built in `trackmate_case_metrics`.
- The same catalog feeds:
  - A/I computation
  - static diagnostic plot
  - GIF diagnostic plot

So blue tracks/red activation markers and CSV counts use the same classification source.

## 3) Script sections in `main_batch_trackmate_arc.m`

## 3.1 Startup and run mode

- `runMode = "arc"` by default.
- `isArc` controls headless plotting behavior:
  - ARC: figures hidden
  - local: figures visible
- `maxTracksToParse` lets local runs use a smaller parse cap while ARC uses all tracks (`Inf`).

## 3.2 Cache setup

- Parser outputs are cached in `all_cases_cache.mat`.
- Cache entries are keyed by:
  - XML path + file signature (size + timestamp)
  - calibration (`pixelSize`, `dt`)
  - parse cap (`maxTracks`)
  - policy tag (parser + gating/activation settings)

This prevents stale outputs when parser or metric policy changes.

## 3.3 Policy knobs

Four option groups control behavior:

- `parserOpts`
  - `parserVersion`, `parseTrackedSpotsOnly`, `parseFilteredTracksOnly`
- `flowOpts`
  - bulk direction and counterflow thresholds
- `qcOpts`
  - min spots, topology limits, optional wall band
- `activationOpts`
  - sustained-growth detection thresholds (2x jump with post-jump persistence + burst fallback)

These are used both for computation and cache-policy hashing.

## 3.4 Case definitions

Each case includes:

- `name`, `Re`, `kDh`
- `xmlFile`, `pixelSize`, `dt`
- optional `diagnosticTrackIds` (empty means all tracks for GIF diagnostics)

## 3.5 Per-case loop

For each case:

1. Build cache key.
2. If compatible cache exists, reuse parsed `out`.
3. Otherwise parse XML using `analyze_trackmate_xml_arc`.
4. Compute metrics using `trackmate_case_metrics(out, qcOpts, flowOpts, activationOpts)`.
5. Append one row to summary table.
6. Append strict legacy gate-audit row.
7. Write framewise per-case CSV.
8. Render static track diagnostics and GIF diagnostics (if enabled).

Then it closes figures and continues.

## 3.6 Final outputs and plots

After all cases:

- Sort rows by `Re`, `kDh`.
- Write:
  - `Summary_tracks.csv`
  - `track_gate_summary.csv`
- Plot:
  - A/I vs kDh (`plot_ai_vs_kdh_re`)
  - Inception map (`plot_inception_locations_by_re`)
  - Tau vs kDh (`plot_tau_vs_kdh_re`)
  - Upstream size distributions (`plot_upstream_size_distribution_by_re`)
- Save updated cache if anything was reparsed.

## 4) Parser logic (`analyze_trackmate_xml_arc.m`)

## 4.1 Input behavior

Key parser switches:

- `parseFilteredTracksOnly=true`:
  - parse only track IDs in TrackMate `<FilteredTracks>` if present
  - if missing/empty, fallback to all tracks with warning
- `parseTrackedSpotsOnly=true`:
  - parse only spot IDs referenced by retained edges

## 4.2 Track/edge extraction

For each retained `<Track>`:

- read track stats (`NUMBER_SPOTS`, `NUMBER_GAPS`, `NUMBER_SPLITS`, `NUMBER_MERGES`, `NUMBER_COMPLEX`, etc.)
- read each `<Edge>` source/target and edge attributes

Outputs:

- `out.tracks` table (track-level metadata)
- `out.edges` table (graph links)

## 4.3 Spot extraction

- Parses `<Spot>` nodes and keeps only required tracked spots when enabled.
- Reads area and geometry features (`AREA`, `PERIMETER`, `CIRCULARITY`, etc.).

Output:

- `out.spots` table

## 4.4 Trajectory reconstruction

For each track ID in edges:

1. Infer start node (`source` not appearing in `target`).
2. Walk source->target links to build ordered spot sequence.
3. Map spot IDs to spot table rows.
4. Build per-track vectors:
   - frame/time/x/y
   - physical coordinates (mm) via `pixelSize`
   - `ds`, `dt`, speed

Output:

- `out.trajectories` struct array

## 4.5 XML robustness

`read_trackmate_xml_document`:

- preflights file readability
- tries `xmlread`
- falls back to Java DOM parser
- detects Java heap OOM and raises explicit guidance (increase JVM heap on ARC)
- adds readable error context and file head preview for malformed/bad XML issues

## 4.6 Parser metadata

`out.meta` stores parser version and parse policy flags.  
Main script checks this before reusing cache entries.

## 5) Metrics logic (`trackmate_case_metrics.m`)

This file produces both:

- **Primary** framewise A/I outputs (used for plots + Summary CSV)
- **Strict legacy** gated outputs (retained for audit/comparison)

## 5.1 Shared per-track prep and catalog

For each trajectory it builds `prep(k)`:

- validity checks (finite x/y/t, monotonic time, min spots)
- left-moving check: `x_end < x_start`
- activation detection from area series

Then it stores `trackCatalog(k)` with:

- track ID, frame/x/y
- basic-valid flag
- left-moving flag
- activation flag + event location/frame/index

This catalog is the single source for metrics and diagnostics.

## 5.2 Primary denominator/numerator definitions

Primary left-moving set:

- pass basic validity
- net left-moving (`x_end < x_start`)

Primary denominator is **frame exposure**:

`leftMovingTrackFrameExposure = sum_t N_leftMovingVisible(t)`

Primary numerator is activation events:

`activationEventsTotal = sum_t N_activationEvents(t)`

with one activation event per track at first sustained-growth jump.

Primary metric:

`A_over_I = activationEventsTotal / leftMovingTrackFrameExposure`

Wilson CI is computed for this ratio.

## 5.3 Sustained-growth activation detector

`find_sustained_growth_activation` applies:

- trigger jump:
  - `A(i+1)/A(i) >= areaJumpFactor` (default 2.0)
  - `A(i+1)/A_pre >= areaJumpFactor`, where `A_pre = median(pre-window)`
- persistence:
  - require `requiredPostFrames` post samples (default 2)
  - `median(post/A_pre) >= postMedianFactor` and `max(post/A_pre) >= postMaxFactor`
- burst fallback:
  - if only one post sample, accept only if stronger thresholds pass
- no post evidence -> reject candidate

## 5.4 Strict legacy branch (secondary)

The legacy gate uses additional physics/topology filters:

- topology gate (`NUMBER_GAPS`, splits/merges/complex)
- counterflow step-fraction gate
- right-origin gate (for left->right bulk flow)
- optional wall-band gate

Outputs are stored with `_strictLegacy` suffix.

## 6) Diagnostics consistency

## 6.1 GIF (`save_diagnostic_track_gifs`)

- Uses `metrics.trackCatalog` directly.
- Blue tracks = `isLeftMoving`.
- Red markers = activation events from `metrics.activationEvent_xy/frame`.
- Title values (`left-moving`, `act@frame`, cumulative exposure/events, A/I) come from framewise arrays in `metrics`.

## 6.2 Static track diagnostics (`plot_verification_tracks_for_case`)

- Draws `metrics.leftMovingTrack_xy`.
- Marks `metrics.activationEvent_xy`.
- Legend text explicitly says these are denominator/numerator sets for framewise A/I.

## 7) CSV outputs explained

## 7.1 `Summary_tracks.csv` (primary table)

Key groups:

- case metadata: `Case, Re, kDh`
- track counts:
  - `nTracksTotal`
  - `nValidTracks`
  - `nLeftMovingTracks`
  - `nActivatedLeftMovingTracks`
- framewise A/I terms:
  - `leftMovingTrackFrameExposure`
  - `activationEventsTotal`
  - `meanLeftMovingPerFrame`
  - `peakLeftMovingPerFrame`
  - `A_over_I` + CI
- strict legacy audit:
  - `nInjected_strictLegacy`
  - `nActivated_strictLegacy`
  - `A_over_I_strictLegacy` + CI
- time/location/size summary stats:
  - tau mean/median/p90/std/count
  - injection and activation location moments
  - upstream equivalent-diameter stats

## 7.2 `track_gate_summary.csv`

Per-case counts for strict legacy rejection reasons:

- too short
- nonfinite
- non-monotonic time
- topology
- flow
- origin
- no activation
- wall-band

Plus strict injected/activated counts and origin-threshold values.

## 7.3 `framewise_counts_<case>...csv`

Frame-by-frame debug columns:

- `frame`
- `nLeftMovingVisible`
- `nActivationEvents`
- `cumExposure`
- `cumActivationEvents`

Useful to directly verify:

`A_over_I = cumActivationEvents(end) / cumExposure(end)`

## 8) Why this structure is robust

- Parser policy and metric policy are encoded into cache keys.
- XML parser has fallback path + OOM-specific diagnostics.
- One shared classification source avoids GIF/static/CSV mismatch.
- Legacy strict outputs remain available for side-by-side audit without changing primary A/I semantics.

