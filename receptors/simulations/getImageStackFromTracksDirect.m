function getImageStackFromTracksDirect(tracksSim,saveDir)

%load('target_rD10_ap0p5_lR0p3/out1/compTracksALT1.mat');
% loadStruct = load([pathToTrackMAT,'/out',int2str(runIndx),'/compTracksALT',int2str(runIndx),'.mat']);

% compTracksALTdef_sparse = loadStruct.compTracksALT.defaultFormatTracks;
%{
     numTracks = length(compTracksALTdef_sparse);

     for indx=1:numTracks
         compTracksALTdef_full(indx).tracksFeatIndxCG = full(compTracksALTdef_sparse(indx).tracksFeatIndxCG);
         compTracksALTdef_full(indx).tracksCoordAmpCG = full(compTracksALTdef_sparse(indx).tracksCoordAmpCG);
         compTracksALTdef_full(indx).aggregState = full(compTracksALTdef_sparse(indx).aggregState);
         compTracksALTdef_full(indx).seqOfEvents = compTracksALTdef_sparse(indx).seqOfEvents;
     end
%}

% [trackedFeatureInfo,trackedFeatureIndx,trackStartRow,numSegments,aggregStateMat] =...
%     convStruct2MatIgnoreMS(compTracksALTdef_sparse,0);


%check directory for saving
if nargin < 2 || isempty(saveDir)
    saveDir = [];
end

%convert tracks structure to matrix
[trackedFeatureInfo,trackedFeatureIndx,trackStartRow,numSegments,aggregStateMat] =...
    convStruct2MatIgnoreMS(tracksSim,0);

%Use a portion of the simulation if desired (truncating)
%i.e, only the first 1000 iterations (10 seconds, at time step = 0.01)
endIter = 8 * 1000;
trackedFeatureInfo = trackedFeatureInfo(:,1:endIter);

%Scale up intensities
trackedFeatureInfo(:,4:8:end) = trackedFeatureInfo(:,4:8:end) * 1000;

%pixel in microns
pixelSize = 0.09;
trackedFeatureInfo(:,1:8:end) = trackedFeatureInfo(:,1:8:end)/pixelSize;
trackedFeatureInfo(:,2:8:end) = trackedFeatureInfo(:,2:8:end)/pixelSize;

%remove 0 columns
trackedFeatureInfo(:,3:8:end) = [];
trackedFeatureInfo(:,4:7:end) = [];
trackedFeatureInfo(:,4:6:end) = [];
trackedFeatureInfo(:,4:5:end) = [];
trackedFeatureInfo(:,4:4:end) = [];

%Replace zeros with NaNs
trackedFeatureInfo(trackedFeatureInfo == 0) = NaN;

%Subsample whole simulation every 10 iterations
xCoord = trackedFeatureInfo(:,1:30:end);
yCoord = trackedFeatureInfo(:,2:30:end);
amps = trackedFeatureInfo(:,3:30:end);
[numTracks,numIters] = size(xCoord);

%alleviate boundary effects on detection: add 10 pixels to
%particle positions; a few lines down image size will be expanded as
%well
xCoord = xCoord + 10;
yCoord = yCoord + 10;

%reformat track information following function input
trackInfo = reshape([xCoord;yCoord;amps],numTracks,3*numIters);

%Set parameters
bgav = 5000;
bgnoise = 50;
sigma = 1;

%Convert 25 micron square to pixels (1 pixel is 90 nm x 90 nm)
%also add 10 pixels to each image size to alleviate boundary effects
%on detection
imsize(1:2) = round(25/pixelSize) + 20;

rad = 3;
saveVar = 1;
% saveFolder = ['outB',int2str(runIndx)];

makeAiryImageFromMPM(trackInfo,bgav,bgnoise,sigma,imsize,rad,saveVar,[],saveDir);

end



