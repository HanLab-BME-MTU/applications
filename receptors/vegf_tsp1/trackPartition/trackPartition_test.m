dataDir = '/project/biophysics/jaqaman_lab/vegf_tsp1/knguyen/test-data/testIllustration';
addpath(genpath(dataDir))
load([dataDir,'/m-01/TrackingPackage/GaussianMixtureModels/detections_for_channel_2/Channel_2_detection_result.mat'])
load([dataDir,'/m-01/TrackingPackage/tracks/Channel_1_tracking_result.mat'])
load([dataDir,'/m-01/m-01.mat'])
%%
tic; [mask,tracks] = trackPartitionInit(tracksFinal,movieInfo,MD,0,100,3); toc
%% test
tic; tracksPart = trackPartition(tracks,mask,1,200,3); toc
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
%% read image
m_01 = imread([dataDir,'/m-01.tif']);
channel1 = m_01(:,1:512);
channel2 = m_01(:,513:end);
channel1 = mat2gray(channel1);
channel2 = mat2gray(channel2);
I = zeros([size(channel1),3]);
I(:,:,1) = channel1;
I(:,:,2) = channel2;
imtool(I)
%%
plotTracksPart(tracksPart,diffAnalysisRes,[],[],I,[],[],[]);
