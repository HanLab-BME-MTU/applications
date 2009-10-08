plusTipSmartSettings README

This README explains how to choose optimal tracking parameters by doing a 
parameter sweep.  

plusTipSmartSettings allows the user to choose a range of values for
each tracking parameter that should be tested, leaving all others constant
at their user-set defaults.

After project setup and feature detection (see plusTipGetTracks), select a
project containing a representative movie for which you would like to 
determine the optimal tracking parameters.

Fill in the frame rate and pixel size.

For each tracking parameter, there is a field called "Default" and a
field called "Range to Test."  The function will test each parameter 
independently, using the range specified, and use the default values for 
all other parameters.  Leave the range blank if you do not want to test a
parameter.

A new folder called "paramTest" will appear under the project folder. In it
will be figure (fig) folders that contain two kinds of .tif images:

    * gapLifetimes_xxxxx.tif: histogram of the lifetime of all gaps (fgaps 
    and bgaps) in frames. The max gap length provides the upper bound.

    * linkingDistance_xx.tif: histogram of the frame-to-frame linking distance,
    which is how far the actual linked feature is from where it was 
    projected to be, given the previous time points of the track.  The max 
    search radius provides the upper bound.

A Matlab cell array called "dataIndivTest" is also saved in the fig folder,
which contains the statistics calculated during post-processing for each
value in the range to test.

In addition, a cell array called "allData" is saved in the paramTest folder.
This contains the "dataIndivTest" in one convenient matrix that can be
imported into an Excel spreadsheet for further analysis.

Optimal parameter selection requires visual inspection of the histograms and
the values recorded in allData.  Good min/max radii will produce linking 
distance histograms where the exponential does not cut off abruptly
at the max radius value.  Similarly, a good max gap length will produce a 
gapLifetimes histogram that tapers off at the max value.  This value can
also be visual inspection of tracks containing gaps.  To choose settings 
for the other parameters, consider plotting how the MT dynamics parameters 
in dataIndivTest change as the tracking parameter of interest changes.

The defaults and ranges given were used to validate parameters for a movie
with a frame rate of 0.8sec and a real-space pixel size of 110nm.