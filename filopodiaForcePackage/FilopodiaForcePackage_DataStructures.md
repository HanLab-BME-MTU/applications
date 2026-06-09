# FilopodiaForcePackage — data structures & wrapper contracts

Channels (this MovieData): **talin-GFP** = structure + tip/base; **bead (red, 40 nm)** = TFM, already processed by `TFMPackage`. All optical processing is centered on the talin channel; force is read from the existing `TFMPackage`.

Each wrapper has the signature `wrapper(movieData, paramsIn)`, resolves its own `Process` via `movieData.getProcess`/`getPackage`, reads inputs from upstream `outFilePaths_`, writes outputs to `funParams.OutputDirectory`, and finally calls `setOutFilePaths`.

---

## Process 1 — `segmentMovieFilopodia`
**Reads:** raw talin channel.
**Writes (per frame):** `bodyMask` (logical), and the steerable maps `res`, `theta`, `nms`, `scaleMap` from `multiscaleSteerableDetector(img, SteerableOrder, SigmaArray)`.
Cell body = intensity threshold + cleanup. Shaft enhancement = hysteresis threshold on `res`; `theta` is kept for orientation-guided linking in P2; `scaleMap` separates thick body from thin shaft.

## Process 2 — `detectMovieFilopodia`
**Reads:** P1 outputs + talin channel.
**Pipeline:** `pointSourceDetection` on talin → bright puncta. Classify by `bodyMask` (distal → tip, body-edge → base, interior → FA, ignored). Global one-to-one tip↔base assignment (LAP) under `MaxTipBaseDist`. Shaft = minimum-cost path between paired endpoints on cost `1./max(res,PathCostFloor)`, biased to agree with `theta` within `OrientTolerance`. Endpoints anchor the path so dim/broken shafts still resolve.
**Writes:** `filoInfo` (struct array, one element per filopodium per frame):

| field | type | meaning |
|---|---|---|
| `tipPos` | 1x2 | sub-pixel tip [x y] |
| `tipAmp` | 1x1 | talin amplitude at tip (psd) |
| `tipPval` | 1x1 | detection p-value |
| `basePos` | 1x2 | base [x y] |
| `baseAmp` | 1x1 | talin amplitude at base |
| `baseAdhesionId` | 1x1 | FA label at base (optional) |
| `centerline` | Nx2 | traced path base→tip [x y] |
| `arc` | Nx1 | cumulative arclength along centerline |
| `length` | 1x1 | total length L (px) |
| `shaftMeanInt` | 1x1 | mean talin along shaft |
| `frame` | 1x1 | frame index |

## Process 3 — `trackMovieFilopodia`
**Reads:** P2 `filoInfo` (all frames).
**Pipeline:** LAP linking of tip (and base) across frames, gap closing, min-length filter. Per track: `L(t)` = geodesic base→tip; `velocity` = smoothed `dL/dt`; `state` from velocity sign/threshold; `fluctFreq` from `L(t)` spectrum.
**Writes:** `filoTracks` (struct array, one per tracked filopodium):

| field | type | meaning |
|---|---|---|
| `trackId` | 1x1 | unique id |
| `frames` | 1xT | frames present |
| `tipTrack` | Tx2 | tip [x y] over time |
| `baseTrack` | Tx2 | base [x y] over time |
| `L` | 1xT | length time series |
| `velocity` | 1x(T-1) | signed dL/dt (+ protrusion, − retraction) |
| `state` | 1xT | 'protrusion'\|'retraction'\|'pause' |
| `lifetime` | 1x1 | frames (or seconds via metadata) |
| `fluctFreq` | 1x1 | dominant fluctuation frequency |
| `centerlineRef` | 1xT | index back into `filoInfo` per frame |

## Process 4 — `windowMovieFilopodia`
**Reads:** P3 tracks (+ centerlines via `filoInfo`).
**Pipeline:** for each track/frame, resample centerline by arclength and lay windows (`normalized` = constant `NumWindows`; `fixed` = `WindowLength` px). Region by normalized `s`: base if `s<=BaseFraction`, tip if `s>=1-TipFraction`, else shaft. Tip window re-anchored to the current distal endpoint each frame.
**Writes:** `filoWindows` (per track, per frame, per window):

| field | type | meaning |
|---|---|---|
| `trackId` | 1x1 | parent track |
| `frame` | 1x1 | frame |
| `winIndex` | 1x1 | 0..NumWindows-1 (base→tip) |
| `region` | char | 'base'\|'shaft'\|'tip' |
| `sCenter` | 1x1 | arclength (and normalized) of window center |
| `tangent` | 1x2 | unit tangent (axial direction) |
| `pixIdxInt` | Mx1 | pixel footprint for intensity (narrow) |
| `pixIdxForce` | Px1 | pixel footprint for force (wider, TFM scale) |

## Process 5 — `sampleMovieFilopodia`
**Reads:** P4 windows; talin channel; traction field from
`movieData.getPackage(... 'TFMPackage').getProcess('ForceFieldCalculationProcess')`.
**Pipeline:** apply green↔red channel registration; for each window sample traction magnitude and project onto `tangent` → axial / lateral; sample talin intensity (`SampleStat`).
**Writes:** `filoSamples` indexed `[trackId, frame, winIndex]`:

| field | meaning |
|---|---|
| `forceMag` | mean |T| in window |
| `forceAxial` | traction · tangent |
| `forceLateral` | traction · normal |
| `talinInt` | talin intensity (SampleStat) |
| `nPix` | pixels averaged |

## Process 6 — `computeMovieFilopodiaStats`
**Reads:** P3 tracks, P4 windows, P5 samples.
**Writes:** `filoStats`: per-frame and mean count; length / velocity / fluctuation-frequency distributions; lifetime; region-resolved force & talin distributions (tip vs shaft vs base); correlations listed in `funParams.Correlations` (e.g. tip talin vs tip axial force; protrusion velocity vs tip force); filopodia density vs body perimeter/area. Optional CSV + figures.

---

### Notes
- `funName_`/`funParams_` are set after the superclass call. If your installed `Process` base accepts `(owner, name, funName, funParams)`, move them into the super call.
- `GUI()` returns `@noSettingsProcessGUI` as a placeholder; replace with a custom settings GUI when ready.
- TFM (bead) spatial resolution, not the Airyscan optics, limits how finely tip/shaft/base force can be separated; `ForceLateralWidth` exists so force windows can match the regularization length scale while intensity windows stay narrow.
