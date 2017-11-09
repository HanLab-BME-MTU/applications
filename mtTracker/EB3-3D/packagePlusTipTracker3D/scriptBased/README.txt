%%%%%%%% plusTipTracker 3D Alpha 1: INSTRUCTIONS AND RELEASE NOTE %%%%%%%%

%% USAGE: 

The four scripts in 'user-script/' should be used to complete the analysis 
process. They are self-described. 

%% HOW TO

MANAGE MOVIES
In order to access both pixels and metadata, each movie is associated to a
data structure called MovieData. Each MovieData structure is saved in a file 
named 'movieData.mat'. The software takes as input a list of movieData 
encapsulated in a structure called MovieList.  

The script <user-script/createMovieManagementFile.m> automatically generates those data 
structures and save their associated file, provided that the file tree follows a
proper structure (see script description). 

The script <user-script/loadMoviesManagementFile.m> demonstrates how to load 
those files in Matlab to be used with the software. 

DETECT AND TRACK
Two approaches are currently proposed for EB* detection. 'PointSourceAutoSigmaLM'
is a parameter-free technique based on image noise analysis, it suitable 
for raw  (but deskewed) data. 'bandPassWatershed' is the 3D version of the 
original plusTipTracker detector (Matov et al 2010), more suitable for
deconvolve data and requires one parameter. 

For each cells, the resulting files are stored in <CellName>/analysis/EB3/<methodName>/. 
The file 'detection.mat' contains the coordinates values. The folder 
'amiraVertex' contains the time series to be loaded in Amira for display. 
 
Tracking use the exact same set of input parameter as in (Applegate et al 2012), 
the only differences are internal. 

For each cells, the resulting files are stored in <CellName>/analysis/EB3/<methodName>/plustipTrackerio. 
The folder 'amiraTrack' contains the time series to be loaded in Amira for display.
In the folder 'track', the file 'trackResults.mat' contains the legacy output 
(track structure, kalman filter estimates), the file 'trackNewFormat.mat' 
cointains a more human readable structure for tracks. 

The script <user-script/scriptEB3DetectAndTracking.m> performs detection and
tracking automatically. It gives a detailed description of the input parameters.

VISUALIZE
Detection and tracks visualization can be performed in Amira 6. 

Original stacks and overlays must be loaded using the <open Time Series Data> dialog.
The pixel size must be inputed manually while opening original stacks. 

The module "volren" must be used with original stack.The module 
"SpatialGraphView" must be used with tracks and detections. 

ANALYZE RESULTS
The script <scriptAnalyseResults> demonstrate how to compare lifetime 
distribution, mean tracks length and mean track speed, accross various 
conditions and cells. 

%% CURRENT RELEASE NOTE 

::Alpha 1::

Features:
- Detection: two techniques are proposed one automatique for raw data
- Tracking: porting of plusTipTracker in 3D
- Visualization: automatic export of detection and tracking to the Amira 6 format

Known issues:
- In-depth evaluation of tracking results are yet to be performed. 
- Loading of Amira .am file has to be done manually in Amira.
- Track analysis is simplistic
- Tested on Linux only

Next Release: 
- Pole detection
- port of growth analysis from plusTipTracker 2D
- Amira template


%% PAST RELEASE DESCRIPTION
