+TIP Detection, Tracking, and Post-Processing README

This README explains how to do the following:
- Set up new projects
- Select one or more projects to analyze together
- Run detection, tracking, and post-processing

---------------------------------------------------------------------------
PROJECT SETUP

A new project begins with a tiff series in a folder called "images."  As
analysis commences, folders called roi_1, roi_2, etc. will be created.  
Multiple projects will only be necessary for a given image series if
there are multiple cells in the field of view (analysis at the sub-cellular
scale is possible with post-processing functions).  Even in that case, it 
is generally ok to analyze them together unless you want to look at cell-cell
heterogeneity.

For a new project, click "Set Up New Project(s)." 
You will be prompted to select a top-level directory.  This can be any
directory as long as there are one or more directories called "images"
containing the tif series you want to analyze.  

If "Draw ROI" is checked, you will be able to draw up to 9 regions-of-interest
(ROIs) per movie.  The first frame of each movie will appear in succession
and allow you to draw a polygons.  Left-clicking to add points and close the
polygon by right-clicking on the first point and selecting "Create mask."  
If "Draw ROI" is unchecked, you will be prompted to click in the center
of the cell on the first frame.  The whole image will then be considered 
to be the ROI for future steps. (The click simply helps the comet detector 
estimate the local background.)

---------------------------------------------------------------------------
PROJECT SELECTION

This step allows you to choose one or more projects (i.e. roi_1 folders) 
for analysis.  This allows you to analyze many movies at once, provided they 
have the same parameters.

If "Select Subset" is checked, a window will pop up asking for one or more
search strings.  These are strings of characters that can be used to narrow
down the number of projects you have to scroll through.  For example, you 
may know that you have some projects that contain "tau" (somewhere in the file
path), so you may input "tau" into the search string list.  Then only those matching
the query will pop up.  If "Select Subset" is unchecked, this step is bypassed
and all the projects will appear in the list.

When the next window pops up with the list of projects, choose one or more 
and use the arrow to move them from the left to the right.

---------------------------------------------------------------------------
TRACKING PARAMETERS

Maximum gap length (frames): max number of frames that can separate two
growth trajectories for them to be considered candidates for gap closing. 
8-10 is generally a good range.

Minimum track length (frames): minimum number of consecutive frames over
which a feature must persist in order to be retained as a track. 3 is
the default.

Search radius range (pixels): min/max distance away from a feature's projected
position (based on its movement in the past) to consider when looking for
candidate features in the next frame.  This will depend on the average speed,
the frame rate, and how much heterogeneity there is between frame-to-frame  
displacements. You will probably want to try several ranges to find the 
optimal values for your data, though the tracker is not super sensitive to 
these values.  (e.g. for 2sec data, try 10-15; for 0.8sec data, try 3-5)

Maximum angle (degrees): only tracks which start within a cone of +/- this 
angle projecting out from the end of a track will be considered as candidates
for pause/out-of-focus events.

Max shrinkage factor: how fast shrinkage is expected to occur relative to
growth.  (e.g. if a MT grows at 10microns/min, and we assume it cannot 
shrink faster than 15microns/min, then this parameter is 1.5.)

Max perp end-start distance (pixels): maximum distance between a track start 
and the nearest point along the current track.  This is used to narrow down
the candidates for backward (shrinkage) linking.  Candidates for shrinkage
should be colinear with the track, even if the track has curvature.
