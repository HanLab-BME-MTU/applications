% to run the function:
colocalizationAdhesionsWithTFM('/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/IntraVsExtraForce/Youbean/130110 Cell7 2Frames/ROIAnalysis',8,false);
colocalizationAdhesionsWithTFM('/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/IntraVsExtraForce/Youbean/130131 cell3paxTIRF 2Frames/ROIAnalysis','L1 2nd',16,false,4000);
colocalizationAdhesionsWithTFM('/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/IntraVsExtraForce/Youbean/130131 cell3paxTIRF 2Frames/ROIAnalysis','L2 0th',16,false,4000);

%%
pathForTheMovieDataFile = '/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM/130429_cell11_100f/ROI';
outputPath = 'tracksSeparation';
band = 16;
showAllTracks = false;
plotEachTrack = false;
tmaxEach = [];
tmax=1600;
[tracksNA,forceFC,forceFA,forceBGband,forceBG] = colocalizationAdhesionsWithTFM( pathForTheMovieDataFile,outputPath,band,tmax,showAllTracks,plotEachTrack,tmaxEach);
%% with L2optimal - 1e-7
pathForTheMovieDataFile = '/Users/joshua2/Documents/PostdocResearch/Traction Force/Manuscript/Draft1.2/Data/130429_cell11_100f/ROI';
outputPath = 'L2optimal';
band = 16;
showAllTracks = true;
plotEachTrack = false;
tmaxEach = [];
tmax=1600;
[tracksNAL2opt,forceFCL2opt,forceFAL2opt,forceBGbandL2opt,forceBGL2opt] = ...
    colocalizationAdhesionsWithTFM( pathForTheMovieDataFile,outputPath,band,tmax,showAllTracks,plotEachTrack,tmaxEach);
%% for L1 again - to be consistent with tracks in L2
outputPath = 'L1l';
band = 16;
showAllTracks = true;
plotEachTrack = false;
tmaxEach = [];
tmax=1600;
[tracksNAL1,forceFCL1,forceFAL1,forceBGbandL1,forceBGL1] = ...
    colocalizationAdhesionsWithTFM( pathForTheMovieDataFile,outputPath,band,tmax,showAllTracks,plotEachTrack,tmaxEach);
[tracksNAL1, tracksNAfailingL1,tracksNAmaturingL1] = ...
    separateMatureAdhesionTracks(tracksNAL1,'/Users/joshua2/Documents/PostdocResearch/Traction Force/Manuscript/Draft1.2/Data/130429_cell11_100f/ROI/Colocalization/L1l');
mean(arrayfun(@(x) x.forceMag(x.emergingFrame),tracksNAL1(arrayfun(@(x) x.emerging, tracksNAL1))))
mean(arrayfun(@(x) x.mean,forceBGL1))
mean(arrayfun(@(x) x.mean,forceBGbandL1))
mean(forceFCL1)
mean(forceFAL1)
outputPath = '/Users/joshua2/Documents/PostdocResearch/Traction Force/Manuscript/Draft1.2/Data/130429_cell11_100f/ROI/Colocalization/L1l/data';
[slopeMaturingMeanL1,slopeFailingMeanL1,InitForceMaturingMeanL1,InitForceFailingMeanL1] = ...
    postColocalizationAdhesionTFM(tracksNAmaturingL1, tracksNAfailingL1, outputPath);
%% for beforeFAKi
pathForTheMovieDataFile = '/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM/130511 Cell2 beforeFAKi/ROItop';
outputPath = 'tracksSeparation';
band = 16;
showAllTracks = false;
plotEachTrack = false;
tmaxEach = [];
tmax=4500;
[tracksNA,forceFC,forceFA,forceBGband,forceBG] = colocalizationAdhesionsWithTFM( pathForTheMovieDataFile,band,tmax,showAllTracks,plotEachTrack,tmaxEach,'outputPath',outputPath);

