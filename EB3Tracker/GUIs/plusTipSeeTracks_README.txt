plusTipSeeTracks README

---------------------------------------------------------------------------

Overview of supported visualization modes:

* Track Overlays: for single projects
    Overlay tracks on a frame of the movie, with the option to 
    select individual tracks for more information about them

* Track Movies: for single projects
    Make a movie of either all tracks within a region, within a frame range 
    or one or more individual tracks.

* Speed Movies: for single projects
    Make a movie where comets are color-coded by speed (microns/min)

* Sub-ROIs: for single or multiple projects
    Create sub-regions-of-interest and automatically extract growth
    sub-tracks from them.

* Quadrant Scatter Plots: for single projects or groups
    Make a scatter plot of two parameters (i.e. growth 
    speed and growth lifetime) divided into quadrants depending on values 
    or percentiles. The colors of the four quadrants correspond to tracks 
    overlaid on an image in subsequent figures.

---------------------------------------------------------------------------

Track Type Quick Reference for Track Overlays and Track Movies
1 - growth                               (red solid)
2 - forward gap  (pause)                 (cyan dotted)
3 - backward gap (shrinkage)             (yellow dotted)
4 - unclassified gap                     (magenta dotted)
5 - forward gap reclassified as growth   (green solid)
6 - backward gap reclassified as pause   (blue dotted)

---------------------------------------------------------------------------
OPTIONAL SETTINGS:
---------------------------------------------------------------------------
1. Select Project(s)

Note: once you choose a single project, it can be used for multiple tasks 
(e.g. track overlay followed by sub-ROI selection).  You can choose a new 
one at any time, or click "Reset" to start over.  

Choose multiple projects to run sub-ROI selection or quadrant scatter plot 
analysis in batch mode.

Groups of selected projects can be created using the "Create Group(s)" 
button.  The output is a file called "groupList," which is used by the the 
Quadrant Scatter Plot tool (batch mode) and some stand-alone functions.

If your data is arranged in a data hierarchy such that projects from
different groups are stored at the same level, you may generate groups
automatically by checking "Auto group from hierarchy."  You will be prompted
to select which levels of the directory tree should be used to create unique
group names.  If this option is unchecked, you will be prompted to choose
groups of projects and name them.  Avoid using spaces and hyphens in the
group names.

The "Create Group(s)" button calls plusTipPickGroups.m.

TROUBLESHOOTING:

* If you have created a roi_x directory but have not run tracking and 
post-processing, it will not appear in the list. 

* Track overlays, track movies, and speed movies can only work with one 
project at a time.  Sub-ROI and quadrant scatter plot analysis can work with
one or many projects at a time.

* If no projects are found, check to make sure there are no spaces anywhere 
in the directory path or file names.

* If you get the message "Select any directory above input directory", the 
root of your Matlab current directory does not match the root directory
where your project is stored.  Point to the relevant server location.

---------------------------------------------------------------------------
2. Select Saved ROI

Click the button if you want to load a saved roiYX.mat file, which contains 
the coordinates of a region you have previously selected.
These are saved during project setup with plusTipGetTracks and also for 
each movie that is generated.  Once you load a ROI, it can be 
used for multiple tasks (e.g. track overlay followed by movie making).  
You can choose a new one at any time, or click "Reset" to start over with 
no ROI.  If no ROI is chosen, the whole image will be used for track 
overlays, you will be prompted to select a new ROI for track movies.

---------------------------------------------------------------------------
3. Choose Frame Range

Default is all frames. For track overlays, partial tracks will be shown if 
they exist partially outside the frame range. For quadrant scatter plots, 
partial tracks are excluded from the analysis if "Remove tracks at start/end" 
is checked.

---------------------------------------------------------------------------
4. Select Output Directory

Select the directory in which to store movies.  Note that overlays are not 
automatically saved.

---------------------------------------------------------------------------
VISUALIZATION MODES
---------------------------------------------------------------------------
"Track Overlays"

All tracks within the frame range will appear as an overlay on an image
chosen by the user (e.g. first frame of frame range).
 
If "Select Tracks" is checked, the user will be prompted to click one or 
more times on the image. Information about the tracks will appear in the 
Matlab command window as follows:

Track: trackNumber   Frame: frame closest to where user selected
    [trackNumber, start frame, end frame, speed (microns/min), track type, lifetime (frames), displacement (pixels)]

Track Types:
1 - growth                               (red solid)
2 - forward gap  (pause)                 (cyan dotted)
3 - backward gap (shrinkage)             (yellow dotted)
4 - unclassified gap                     (magenta dotted)
5 - forward gap reclassified as growth   (green solid)
6 - backward gap reclassified as pause   (blue dotted)

The track numbers selected will then show up in a new text window below the 
"Plot tracks" button. These are useful if, for example, you want to quickly
make a movie of the track you selected.  You may also see all compound track 
profiles by loading the projData.mat file from the 'meta' folder and viewing:
projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix

The "Plot Tracks" button calls plusTipPlotTracks.m.

---------------------------------------------------------------------------
"Track Movies"

The "Detected Comet Display Options" radio buttons and the "Display Tracks"
check box control how the detected comets and tracks are displayed in the 
movie:

    * All comets, current frame only: displays ALL the detected comets from
    a given frame in that frame only.

    * All comets, all frames: displays ALL the detected comets (ie including
    those that did not get incorporated into a track), color-coded by frame.  
    This option is useful for checking whether a tracking mistake might be 
    due to a missed detection or to a wrong link, for example.

    * Comets in tracks only, all frames: displays only the comets used in the
    tracks, color-coded by frame, such that comets in a track appear shortly 
    before and after a track.

    * None: use this option if you want to make a movie of the raw data 
    or if you only want to show the track without the detected comets.

The "Individual Track Numbers" text box can be used to make movies of
individual tracks.  The track numbers correspond to those found in the first
column of projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix,
the matrix containing the tracking results after post-processing.

You may find it useful to select tracks using the Track Overlays tool and 
copy and paste the track numbers into this text box. Or, load projData 
manually and look for interesting tracks to plot.

Please note that the individual track movies are still bounded by the frame
range given and the frames in which the track exists, i.e. if the frame range
chosen in Step 3 is 10-20, and the individual track of interest goes from 
frame 15 to frame 30, the movie will only contain frames 15-20.

If the "Individual Track Numbers" text box is empty, all tracks will
be shown for the ROI.

The "Dual panel with raw images" function creates a movie where the raw
image is shown on the left and the detection and/or track overlay is shown
on the right.  
 
The "Save as AVI" check box determines whether the movie will be saved as 
.MOV (default) or .AVI.  The AVI option crashes in some versions of Linux, 
so it is advised to leave this box unchecked when working in Linux.
 
The "Make Track Movie" button calls plusTipTrackMovie.m.

---------------------------------------------------------------------------
"Speed Movies"

"Speed Limit" is the maximum speed used in the jet color map (e.g. an input 
of 20 will map all speeds faster than 20 to 20 and range from dark blue at 0
to deep red at 20). The default option (max) uses the whole range.  

Circles   - growth
Triangles - forward gap (fgap)
Squares   - backward gap (bgap)

The "Save as AVI" check box determines whether the movie will be saved as 
.MOV (default) or .AVI.  The AVI option crashes in some versions of Linux, 
so it is advised to leave this box unchecked when working in Linux.

The "Make Speed Movie" button calls plusTipSpeedMovie.m.

---------------------------------------------------------------------------
"Sub-ROIs"

Sub-regions of the cell can be selected in manual or automatic mode. In 
manual mode, the user may select a variable number of ROIs. If the regions
overlap when selected, they will be automatically adjusted so no overlap
occurs during track extraction.  

In automatic mode, the cell is split into a central and a peripheral sub-ROI.
The peripheral region can be further sub-divided by checking the "Also divide
periphery into quadrants" option.  The thickness of the peripheral band is
chosen by the user in microns or as a fraction of the largest distance from
the cell edge to the center of mass.  Thus, automatic sub-ROI selection 
creates 2, 4, or 5 sub_x folders.

If it is desirable to exclude tracks from some region of the cell (e.g. from
a previously-selected sub-ROI), check the "Choose exclude regions" option
and either load a mask or draw the region(s) for exclusion when prompted.

Next, define how long a track must exist in the ROI to be included.
This duration is given either as a fraction of the track's lifetime or as 
some number of seconds (see dropdown menu).

To begin, select one or more projects and press "Select Sub-ROIs."  
A 'subROIs' folder will be created under the roi_x directory and will 
contain info for all sub-ROIs.

Previously-created sub-ROI projects (sub_x) may be included during project 
selection; for these, new sub-regions cannot be selected, but tracks will 
be re-extracted according to the lifetime fraction/seconds.

Sub-ROI 'meta' folders will contain data for GROWTH PHASES ONLY pulled from 
the original ROI's data.

Sub-ROI projects can now be selected for track overlays, movie making, etc.

---------------------------------------------------------------------------
"Quadrant Scatter Plots"

Use the Quadrant Scatter option to color-code tracks falling within 
specified ranges of various parameters.  Select parameters to be plotted
from the x- and y-axis drop-down menus, such as growth speed and growth 
lifetime.  Adjust the data values or percentiles for each parameter 
independently and provide min/max limits (if desired) for each.  Data 
outside this range will be excluded from the analysis.

Because the values on the x- and y- axes must be paired, only certain 
combinations of parameters work.  The track type (e.g. "fgap") must be
the same for x- and y- axes.

If "Remove tracks at start/end" is checked, any track not entirely
contained within the frame range will be excluded.  (Lifetime measurements
can be biased especially in short movies where most long tracks will exist 
at the beginning or end, thereby getting discarded.) 
If this option is unchecked, any track which ends before the frame range 
begins or begins after the frame range ends will be excluded.

If projects from different groups should be compared, use the "Batch process
on groups" option and select the appropriate groupList (see Step 1 above).

Seven figures for each project will appear: 
- a scatter plot 
- five images with tracks overlaid (four colors separately and together)
- a percentage bar plot

For the track overlays, the colors of the tracks correspond to the color 
map of the scatter plot. For example, if we take the 50th percentile each 
for growth speed and growth lifetime we will see four populations in four 
colors: fast and short-lived, slow and short-lived, fast and long-lived, 
and slow and long-lived.  The four populations will appear separately in 
four images and merged together in a fifth image.  The percentage bar plot
will show the relative proportion of the four populations.

If running in batch mode, summary percentage bars and the raw data of the 
four colors will be saved for each group.  The percentage bars will be
stacked in the order of the group names (grp1, grp2 etc.), and the data
will be stored in "btwGrpQuadStats" file.  To speed up processing during 
batch mode, choose the "Make summary plots only" option to bypass making
track overlays.

Note: It is also possible to divide the population of tracks based on one 
parameter into three groups.  For example, if we choose growth speed for
both the x- and y- axes, and select the 25th and 50th percentiles, 
respectively, we will see three populations in three colors: tracks in Q1, 
tracks in Q4, and tracks in both Q2 and Q3.  In this case one figure will 
simply show the raw image.

