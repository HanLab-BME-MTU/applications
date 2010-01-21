plusTipParamSweepGUI README

This README explains how to choose optimal tracking parameters by doing a 
parameter sweep.  

plusTipParamSweepGUI allows the user to choose a range of values for
each tracking parameter that should be tested, leaving all others constant
at their user-set defaults.

After project setup and particle detection (see plusTipGetTracks), select a
project containing a representative movie for which you would like to 
determine the optimal tracking parameters.

Fill in the frame rate and pixel size.

For each tracking parameter, there is a field called "Default" and a
field called "Range to Test."  The function will test each parameter 
independently, using the range specified, and use the default values for 
all other parameters.  Leave the range blank if you do not want to test a
parameter.  The format for a range is s:i:e, where s is the starting value,
i is the increment to increase by during each round, and e is the ending
value.

A new folder called "paramTest" will appear under the project folder. 
In it will be figure (fig) folders, one for each parameter and one called 
figs_allData.  The latter contains summary plots for how a number of tracking
statistics vary as a function of parameter value.  The individual parameter
folders (e.g. figs_fluctRad) contain two kinds of .tif images:

    * gapLifetimes_xxxxx.tif: histogram of the lifetime of all gaps (fgaps 
    and bgaps) in frames. The maximum gap length provides the upper bound.

    * linkingDistance_xx.tif: histogram of the frame-to-frame linking distance,
    which is how far the actual linked particle is from where it was 
    projected to be, given the previous time points of the track.  The max 
    search radius provides the upper bound.

A Matlab cell array called "dataIndivTest" is also saved in the fig folder,
which contains the statistics calculated during post-processing for each
value in the range to test.

In addition, a cell array called "allData" is saved in the paramTest folder.
This contains the "dataIndivTest" in one convenient matrix that can be
copied into an Excel spreadsheet for further analysis.

Optimal parameter selection requires visual inspection of the allData
figures and the two kinds of histograms, e.g. good min/max radii will produce 
a linking distance histogram where the exponential does not cut off abruptly
at the max radius value.  Visual inspection of tracks containing gaps is a
good way to check whether the parameter settings are appropriate (see
plusTipSeeTracks).  In general, the max gap length and the max shrinkage
factor are the most sensitive parameters.

