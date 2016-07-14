dataDir = '/project/biophysics/jaqaman_lab/vegf_tsp1/knguyen/test-data/testIllustration/';
addpath(genpath(dataDir))
movie = 'm-07';
load([dataDir,movie,'/TrackingPackage/GaussianMixtureModels/detections_for_channel_2/Channel_2_detection_result.mat'])
load([dataDir,movie,'/TrackingPackage/tracks/Channel_1_tracking_result.mat'])
load([dataDir,movie,filesep,movie,'.mat'])
%%
tic; [mask,tracks] = trackPartitionInit(tracksFinal,movieInfo,MD,0.001,100,3); toc
%% test
tic; tracksPartTemp = trackPartition(tracks,mask,3,1,200); toc
%% Un-upscale tracks
tracksPart = trackScalar(tracksPartTemp,1/3);
if ~exist([dataDir,movie,'/TrackPartition/'],'dir')
    mkdir([dataDir,movie,'/TrackPartition/'])
end
save([dataDir,movie,'/TrackPartition/tracksPart.mat'],'tracksPart');
%% longest inside track
longest = 1;
totalLength = 0;
totalInsideTracks = 0;
for iTrack = 1:size(tracksPart,1)
    if tracksPart(iTrack).isInside == 1
        length = size(tracksPart(iTrack).tracksFeatIndxCG,2);
        if length > longest
            longest = length;
        end
        totalLength = totalLength+length;
        totalInsideTracks = totalInsideTracks+1;
    end
end
meanLength = totalLength/totalInsideTracks;
%%
tic; diffAnalysisRes = trackDiffusionAnalysis1(tracksPart); toc
save([dataDir,movie,'/TrackPartition/diffAnalysisRes.mat'],'diffAnalysisRes');
%% read image
m_01 = imread([dataDir,filesep,movie,'.tif']);
channel1 = m_01(:,1:512);
channel2 = m_01(:,513:end);
channel1 = mat2gray(channel1);
channel2 = mat2gray(channel2);

% black background
% I = zeros([size(channel1),3]);
% I(:,:,1) = channel1;
% I(:,:,2) = channel2;

% white background
I = ones([size(channel1),3]);
I(:,:,2) = ones(size(channel1))-channel1; 
I(:,:,3) = ones(size(channel1))-channel1-channel2; 
I(:,:,1) = ones(size(channel1))-channel2;
imshow(I)
%%
plotTracksPart(tracksPart,diffAnalysisRes,[],[],I,[],[],[],'inside');
