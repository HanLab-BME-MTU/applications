% TRACKPARTITION
% Version 1.0 22-07-2016
% Files
%   constantSizeMasking      - Creates an image mask based on detection results
%   cutTracks                - Using the intersection information from trackPartitionInner, cut
%   gapInterp                - Interpolates particle position during gaps in a track (NaN's)
%   gaussianMasking          - Creates an image mask based on Gaussian detection
%   mergeShortTracks         - Merge together partition tracks that are shorter than a
%   plotTracksPart           - plots tracks in 2D highlighting the different diffusion
%   plotTracksPartGUIWrapper - Get channels
%   trackPartition           - Performs track partitioning, i.e. dividing particle tracks
%   trackPartitionInit       - Prepares mask and track struct for use by
%   trackPartitionInner      - Does the meat of the computation for trackPartition
%   TrackPartitionProcess    - Process for trackPartition
%   trackPartitionWrapper    - Pass variables and parameters from TrackPartitionProcess to its
%   trackScalar              - Multiply track coordinates by a scalar
%
% GUI Help
%
% Tabs should be used from left to right:
% 1. Load Files
%     	*Click 'Load MovieList' at the top and select a Movie List or Combined Movie List .mat file containing movies with tracks to be partitioned. These movies must have already been run with TrackingProcess.
% 	*Click 'Load MovieList' at the bottom and select a Movie List or Combined Movie List .mat file containing movies with detected structures that will be used to generate masks. These movies must have already been run with SubResolutionProcess.
% 	*The movies in this Movie List or Combined Movie List will populate the boxes on the left. To create a "job", select corresponding track and mask movies in these boxes. For the track movie, select the channel containing the tracks to be partitioned. For the mask movie, select the channel containing the detected structures to be made into masks. Click "Add pair to job list" to add the job to the box on the right.
% 	*If the movies in the Track MD and Mask MD boxes are in corresponding order, multiple jobs can be added automatically using "Add all pairs". The currently selected channels will be applied to all jobs.
% 	*Asterisks(*) next to a job indicate that the track movie has not been analyzed before. Plus signs (+) indicate that this analysis has been done before. These symbols will change to tildes (~) after the job is run.
% ----------------
% 2. Set Params
% Set parameters for analysis here. If TrackPartitionProcess has been run on the track movie before, the previous parameters will be shown here. Otherwise, default parameters will be shown.
% 	*Min. Track Length: Minimum number of frames a particle must remain outside or inside a masked location for its track to be partitioned; e.g. if a particle enters a masked location and exits only after x frames, those x frames will still be considered as "outside".
% 	*Min. Mask Diameter: Minimum diameter of each mask in nm. Useful when the typical size of the masked structure is known. 
% 	*Gaussian Threshold: Set to some nonzero number to incorporate Gaussian fitting information into masking. Each structure will be masked where its fitted Gaussian exceeds this threshold; i.e. set to smaller values (such as 0.001) to create larger masks. Set to 0 to ignore Gaussian information and make every mask the same size, dictated by Min. Mask Diameter.
% 	*Upscale: Because masks can be very small (1 or 2 pixels in the original image resolution) and motion can be sub-pixel scale, upscaling masks and tracks can improve accuracy. This parameter scales up the number of pixels in the image by the specified integer as well as the masks and track coordinates accordingly. Upscaling past 3x may consume ungodly amounts of memory.
% 	*Analysis Start and End Frames:  Defaults to first and last frames in the movies.
% 	*Perform Diffusion Analysis: Run trackDiffusionAnalysis1 on the resulting partitioned tracks.
% ----------------
% 3. Run
% Click Run.
% 
% "Turbo" knob determines number of parpool workers to use. Currently, there are some memory issues with using too many.
% 
% 
% ----------------
% 4. Plot
% Plot the partitioned tracks of the selected job. Diffusion Analysis must have been run.
% 
% Colors are the same as in plotTracksDiffAnalysis2D (this is basically a ripoff of that function). "Inside" tracks are shown thicker than "outside" tracks. 
% 
% 	*Start and End Frame: Defaults to first and last frames.
% 	*Show Image: Show first frame of the track and mask movies as an RGB image.
% 	*Use White Background: Makes background of image white instead of black. Useful for when plotting only inside tracks and they're tiny and impossible to spot. (What is this, a plot for ants??)
% 	*Show Confinement Radii: Same as in plotTracksDiffAnalysis2D. 
% 	*Simplify Linear Groups: Same as in plotTracksDiffAnalysis2D. "All" makes all linear groups red. "+Random & Superdiffusive" makes all linear groups and random & superdiffusive groups red.
% 	*Plot Subset: Show "All" tracks or a subset of only "Inside" or "Outside" tracks. 
% 




                    