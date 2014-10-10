% to run the function:
colocalizationAdhesionsWithTFM('/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/IntraVsExtraForce/Youbean/130110 Cell7 2Frames/ROIAnalysis',8,false);
colocalizationAdhesionsWithTFM('/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/IntraVsExtraForce/Youbean/130131 cell3paxTIRF 2Frames/ROIAnalysis','L1 2nd',16,false,4000);
colocalizationAdhesionsWithTFM('/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/IntraVsExtraForce/Youbean/130131 cell3paxTIRF 2Frames/ROIAnalysis','L2 0th',16,false,4000);

%%
pathForTheMovieDataFile = '/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM/130429_cell11_100f/ROI';
outputPath = '141009Denser';
band = 16;
showAllTracks = true;
tmaxEach = false;
tmax=1600;
[tracksNA,forceFC,forceFA,forceBGband,forceBG] = colocalizationAdhesionsWithTFM( pathForTheMovieDataFile,outputPath,band,tmax,showAllTracks,plotEachTrack,tmaxEach);