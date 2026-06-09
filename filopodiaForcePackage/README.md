# FilopodiaForcePackage

lccb/uTrack-style MATLAB package for filopodia mechanics from **talin-GFP**
(structure + tip/base) and a **bead TFM** channel. Mirrors the conventions of
`FocalAdhesionPackage` / `WindowingPackage` (4-arg `Process` super calls,
`getDefault<X>Params` injectors in the Package, `funName_` wrappers driven by
`parseProcessParams` + `setInFilePaths`/`setOutFilePaths`).

## Process chain
| # | Process | Wrapper | Base class | Status |
|---|---|---|---|---|
| 1 | FilopodiaSegmentationProcess | segmentMovieFilopodia | ImageAnalysisProcess | implemented |
| 2 | FilopodiaDetectionProcess | detectMovieFilopodia | DetectionProcess | implemented |
| 3 | FilopodiaTrackingProcess | trackMovieFilopodia | DataProcessingProcess | wired stub |
| 4 | FilopodiaWindowingProcess | windowMovieFilopodia | DataProcessingProcess | wired stub |
| 5 | FilopodiaSamplingProcess | sampleMovieFilopodia | DataProcessingProcess | wired stub |
| 6 | FilopodiaStatisticsProcess | computeMovieFilopodiaStats | DataProcessingProcess | wired stub |

Helpers (implemented): `traceFilopodiaShaft.m` (orientation-guided min-cost
shaft trace between anchored tip/base), `pairFilopodiaTipBase.m` (global
tip<->base LAP using `lap.m`).

## Force source (Process 5)
Traction is read cross-package from the existing `TFMPackage`:
`ForceFieldCalculationProcess` output `tMapUnshifted` (the same output the
WindowingPackage samples). No re-computation of TFM here.

## Notes
- `GUI()` returns `@abstractProcessGUI` as a placeholder; build custom GUIs later.
- Put this folder on the path alongside the FA/uTrack package (needs
  `multiscaleSteerableDetector`, `pointSourceDetection`, `lap`,
  `createDistanceMatrix`, `thresholdRosin/Otsu`, `mkClrDir`, `parseProcessParams`).
- Register on a MovieData (MD) e.g.:
  `MD.addPackage(FilopodiaForcePackage(MD)); pkg = MD.getPackage(MD.getPackageIndex('FilopodiaForcePackage'));`
  then create each process from `getDefaultProcessConstructors` and `.run`.

See FilopodiaForcePackage_DataStructures.md for the on-disk struct fields.
