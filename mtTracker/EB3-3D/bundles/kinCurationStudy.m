load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/cell1_12_half_volume_double_time/movieData.mat')

%% load spindle ref data and randomization 
processSpindleRef=MD.getProcess(1);
processRandKin=MD.getProcess(2);

kinTracks=load(poleRefProcess.outFilePaths_{2}); kinTracks=kinTracks.kinTracks;
EB3TracksInliers=load(processSpindleRef.outFilePaths_{3}); EB3TracksInliers=EB3TracksInliers.EB3TracksInliers;
kinTracksInliers=load(processSpindleRef.outFilePaths_{4}); kinTracksInliers=kinTracksInliers.kinTracksInliers;

randKinTracks=load(processRandKin.outFilePaths_{1});  randKinTracks=randKinTracks.randKinTracks;

randKinTracksInlier=randKinTracks(logical(arrayfun(@(t) t.inliers(1),kinTracks)));

%%
% It turns out that the response to bundling study is very dissimetric, the
% first few kin show strong differences but the not the overall responses.
% Let s test various part of the kinetochore.
kinIndxBoundaries=1:100:length(kinTracksInliers);
parfor b=1:length(kinIndxBoundaries)-1 
    MTBundling(MD,kinTracksInliers,EB3TracksInliers,randKinTracksInlier, 'kinRange',kinIndxBoundaries(1):kinIndxBoundaries(b+1))
end