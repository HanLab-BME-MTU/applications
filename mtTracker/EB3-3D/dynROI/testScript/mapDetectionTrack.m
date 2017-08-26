%% testing ROI mapping for detection and tracks to deprecate the pole oriented code. 
MD=MovieData.loadMatFile('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/cell1_12_half_volume_double_time/movieData.mat');
tmp=load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/cell1_12_half_volume_double_time/EB3/detection/detectionLabRef.mat');
tmpPole=load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/cell1_12_half_volume_double_time/poles/simplex_scale_003/poleDetection.mat');

[~,tracks]=detectPoles(MD,'isoOutput',true);

[Refs,ROIs]=buildRefsFromTracks(tracks,tracks,'buildROI',true);
%%
detections=mapDetectionsToROI(tmp.detectionsLabRef,ROIs{1,2},0,'distType','vertexDistOtsu');

%% 
tmp=load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/cell1_12_half_volume_double_time/EB3/track/tracksLabRef.mat');
tmpAugmented=load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/cell1_12_half_volume_double_time/EB3/track/augmentedSpindleRef.mat');

%%
mappedTracks=mapTracksTo1DManifold(ROIs{1,2},tmp.tracksLabRef,0,'position','start','distType','vertexDistOtsu');
EB3Inlier=tmp.tracksLabRef(arrayfun(@(e) e.inliers(1)==1,tmpAugmented.EB3Tracks));

amiraWriteTracks('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/cell1_12_half_volume_double_time/sandbox/tracks/tracks.am',tmp.tracksLabRef,'cumulativeOnly',true)
amiraWriteTracks('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/cell1_12_half_volume_double_time/sandbox/tracksInlier/tracks.am',mappedTracks,'cumulativeOnly',true)
amiraWriteTracks('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/cell1_12_half_volume_double_time/sandbox/tracksInlier/tracksOrig.am',EB3Inlier,'cumulativeOnly',true)