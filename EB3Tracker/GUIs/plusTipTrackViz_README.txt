+TIP Tracking Visualization README

After detection, tracking and analysis, you can visualize tracks in 3 ways:

- Overlay tracks on a frame of the movie (4a)
- Make a movie of EITHER all tracks within a region and within a frame range 
                  OR     one or more individual tracks (4b)
- Make a movie where tracks are color-coded by velocity (microns/min) (4c)

---------------------------------------------------------------------------

Track Type Quick Reference - see details below
1 - growth                               (red)
2 - forward gap  (pause OR out of focus) (cyan)
3 - backward gap (shrinkage)             (yellow)
4 - unclassified gap                     (magenta)

---------------------------------------------------------------------------

STEPS:

1. Select Project Data (OPTIONAL) 
Click the button and navigate to the projData file corresponding to the 
project of interest (located in the roi_x/meta folder).  Once you choose a 
project, it can be used for multiple tasks (e.g. 4a followed by 4b).  You 
can choose a new one anytime, or click "Reset" to start over. If you do not
choose the project first, you can still run 4a and 4b, but you will have to
select the project for each one separately.

---------------------------------------------------------------------------
2. Select Saved ROI (OPTIONAL) 
Click the button if you want to use a saved roiYX or roiMask.tif file.
These are saved during the intitial ROI selection and also for each movie
that is generated.  Once you load one ROI, it can be used for multiple
tasks (e.g. 4a followed by 4b).  You can choose a new one anytime, or click
"Reset" to start over with no ROI.  If no ROI is chosen, 4a will use the 
whole image, while 4b and 4c will ask you to select a new ROI.

---------------------------------------------------------------------------
3. Choose Frame Range (OPTIONAL)
Default is all frames. The frame range chosen here will be applied to 4a, 
4b, and 4c.

---------------------------------------------------------------------------
4a. Overlay tracks on image
All tracks within the frame range will appear as an overlay on an image
chosen by the user (e.g. first frame of frame range). 
If "Select Tracks" is checked, the user will be prompted to click one or 
more times on the image. Information about the tracks will appear in the 
Matlab command window as follows:

Track: trackNumber   Frame: frame closest to where user selected
    trackNumber startFrame endFrame velocity trackType

where velocity is in microns/min and trackType is:
1 - growth                               (red)
2 - forward gap  (pause OR out of focus) (cyan)
3 - backward gap (shrinkage)             (yellow)
4 - unclassified gap                     (magenta)

The track numbers selected will then show up in a new text window below the 
"Plot tracks" button.

---------------------------------------------------------------------------
4b. Make movie: track overlays
The "Show Tracks" and "Show Detection" check boxes control whether or not
the tracks and detection are overlaid on the movie. The "Save as AVI" check
box determines whether the movie will be saved as .MOV (default) or .AVI.
In Linux the AVI option tends to crash for some reason, so leave unchecked.
If the "Individual Track Numbers" text box is empty, all tracks will be 
shown; otherwise multiple movies will be made for each track number. The
user may find it useful to select tracks in 4a and copy and paste the track
numbers into this text box.

---------------------------------------------------------------------------
4c. Make movie: speeds
Speed Limit is the max speed used in the jet colormap (e.g. an input of 20
will map all speeds lower than -20 microns/min to -20 and all speeds higher
than 20 to 20, ranging from dark blue to deep red.  Default ("max") will use
the whole range.  Negative speeds correspond to shrinkage velocities, 
because they go in the opposite direction.  See 4b for note on AVI movies.




