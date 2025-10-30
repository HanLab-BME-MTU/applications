% figure5_maturingNAinR8: Figure 5 of NA recruitment paper
% Showing the overlay of NAs that mature to FCs and FAs

%% Load MD
% MDstruct=load('/project/bioinformatics/Danuser_lab/P01adhesion/analysis/Sangyoon/NA_RecruitmentProject/Kevin/20200819_imcdTalinR8/TalinR8_012/TalinR8_012.mat');
MDstruct=load('/project/bioinformatics/Danuser_lab/P01adhesion/analysis/Sangyoon/NA_RecruitmentProject/Kevin/20200814_IMCDKO_TalinWT_5kPa/20200814_imcdTalinWT/TalinWT_013/TalinWT_013.mat');
% MDstruct=load('/storage/disk2/Kevin/2017-05-17/ChoK1_TalinShRNA_TalinWT_20170517_1520_006/movieData.mat');
% MDstruct=load('/storage/disk2/Kevin/2017-04-24/ChoK1_TalinShRNA_R8Mut_1520_005/movieData.mat');
MD = MDstruct.MD;

%% Get the tracks
[tracksNA,adhAnalProc]=getTracksNAFromMD(MD);

%% Get the indices
[indMature,indMatureNAtoFA,indMatureNAtoFC,indMatureFCtoFA,...
 indFail,indFailFC,indFailFA, indStableNA,indStableFC,indStableFA,...
 pNAtoFC,pNAtoFA,pFCtoFA]=findAdhesionMaturation(tracksNA);

%% Show over the image
iFrameInterest=round(0.8*MD.nFrames_);
[imgStack,tMap] = getAnyStacks(MD);
%% Show segmentation
iiformat = '%.3d';
p=adhAnalProc.funParams_;
labelTifPath = [p.OutputDirectory filesep 'labelTifs'];
labelAdhesion = imread(strcat(labelTifPath,'/label',num2str(iFrameInterest,iiformat),'.tif'));
adhBound = bwboundaries(labelAdhesion>0,4,'noholes'); % strongly assumes each has only one boundary

% Get the adhesion ID information
% s = struct2table(tracksNA);
refineFAID_cell = arrayfun(@(x) x.refineFAID,tracksNA,'unif',false);

% refineFAID_cell is numTracks x
% refineID_for_everyFramesInvolved. So for each track (each
% row), I'll make each raw a full frame entries although it
% is a bit memory intensive
maxFrame = max(cellfun(@length,refineFAID_cell));
insuffRows = cellfun(@(x) length(x)<maxFrame,refineFAID_cell);
for k=find(insuffRows')
    refineFAID_cell{k} = [refineFAID_cell{k} ...
                NaN(1,maxFrame-length(refineFAID_cell{k}))];
end
refineFAID = cell2mat(refineFAID_cell);
%% Show image
mapFig = figure; imshow(imgStack(:,:,iFrameInterest),[]), hold on
idCurrent=arrayfun(@(x) iFrameInterest>=x.startingFrameExtra & iFrameInterest<=x.endingFrameExtra,tracksNA);

% htrackLine = arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color','b'),tracksNA(indMatureNAtoFC & idCurrent),'UniformOutput',false);
markerSize = 4;
% htrackCircles = arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','MarkerSize',markerSize,'Color','b'),...
%     tracksNA(indMatureNAtoFC & idCurrent),'UniformOutput',false);

%% Draw the adhesion boundaries
% validAdhState = refineFAID(indMatureNAtoFC & idCurrent,iFrameInterest); %cellfun(@(x) x(iFrame),refineFAID(validState));
% cellfun(@(x) plot(x(:,2),x(:,1),'Color',[0.2 0.2 1]),adhBound(validAdhState(~isnan(validAdhState) & validAdhState>0)),'UniformOutput',false);

colorNotMuture=[.2 0.2 .8]; %includes all not-maturing, stable adhesions
colorMatureNA = [.9 0.2 .2];
colorMatureFC = [.7 0.6 0.2];
colorFA = [.7 0.4 0.2];
% For NAs maturing to FAs
arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colorMatureNA),tracksNA(indMatureNAtoFA & idCurrent),'UniformOutput',false);
htrackG1 = arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','MarkerSize',markerSize,'Color',colorMatureNA),...
    tracksNA(indMatureNAtoFA & idCurrent),'UniformOutput',false);
validAdhState = refineFAID(indMatureNAtoFA & idCurrent,iFrameInterest); %cellfun(@(x) x(iFrame),refineFAID(validState));
cellfun(@(x) plot(x(:,2),x(:,1),'Color',colorMatureNA),adhBound(validAdhState(~isnan(validAdhState) & validAdhState>0)),'UniformOutput',false);

% For NAs failing
arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colorNotMuture),tracksNA(indFail & idCurrent),'UniformOutput',false);
htrackG2 = arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','MarkerSize',markerSize,'Color',colorNotMuture),...
    tracksNA(indFail & idCurrent),'UniformOutput',false);
validAdhState = refineFAID(indFail & idCurrent,iFrameInterest); %cellfun(@(x) x(iFrame),refineFAID(validState));
cellfun(@(x) plot(x(:,2),x(:,1),'Color',colorNotMuture),adhBound(validAdhState(~isnan(validAdhState) & validAdhState>0)),'UniformOutput',false);

% For FCs failing
arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colorNotMuture),tracksNA(indFailFC & idCurrent),'UniformOutput',false);
htrackG3 = arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','MarkerSize',markerSize,'Color',colorNotMuture),...
    tracksNA(indFailFC & idCurrent),'UniformOutput',false);
validAdhState = refineFAID(indFailFC & idCurrent,iFrameInterest); %cellfun(@(x) x(iFrame),refineFAID(validState));
cellfun(@(x) plot(x(:,2),x(:,1),'Color',colorNotMuture),adhBound(validAdhState(~isnan(validAdhState) & validAdhState>0)),'UniformOutput',false);

% For FAs failing
arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colorNotMuture),tracksNA(indFailFA & idCurrent),'UniformOutput',false);
htrackG4 = arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','MarkerSize',markerSize,'Color',colorNotMuture),...
    tracksNA(indFailFA & idCurrent),'UniformOutput',false);
validAdhState = refineFAID(indFailFA & idCurrent,iFrameInterest); %cellfun(@(x) x(iFrame),refineFAID(validState));
cellfun(@(x) plot(x(:,2),x(:,1),'Color',colorNotMuture),adhBound(validAdhState(~isnan(validAdhState) & validAdhState>0)),'UniformOutput',false);

% For NAs stable
arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colorNotMuture),tracksNA(indStableNA & idCurrent),'UniformOutput',false);
htrackG5 = arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','MarkerSize',markerSize,'Color',colorNotMuture),...
    tracksNA(indStableNA & idCurrent),'UniformOutput',false);
validAdhState = refineFAID(indStableNA & idCurrent,iFrameInterest); %cellfun(@(x) x(iFrame),refineFAID(validState));
cellfun(@(x) plot(x(:,2),x(:,1),'Color',colorNotMuture),adhBound(validAdhState(~isnan(validAdhState) & validAdhState>0)),'UniformOutput',false);

% For FCs stable
arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colorNotMuture),tracksNA(indStableFC & idCurrent),'UniformOutput',false);
htrackG6 = arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','MarkerSize',markerSize,'Color',colorNotMuture),...
    tracksNA(indStableFC & idCurrent),'UniformOutput',false);
validAdhState = refineFAID(indStableFC & idCurrent,iFrameInterest); %cellfun(@(x) x(iFrame),refineFAID(validState));
cellfun(@(x) plot(x(:,2),x(:,1),'Color',colorNotMuture),adhBound(validAdhState(~isnan(validAdhState) & validAdhState>0)),'UniformOutput',false);

% For FAs stable
arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colorFA),tracksNA(indStableFA & idCurrent),'UniformOutput',false);
htrackG7 = arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','MarkerSize',markerSize,'Color',colorFA),...
    tracksNA(indStableFA & idCurrent),'UniformOutput',false);
validAdhState = refineFAID(indStableFA & idCurrent,iFrameInterest); %cellfun(@(x) x(iFrame),refineFAID(validState));
cellfun(@(x) plot(x(:,2),x(:,1),'Color',colorFA),adhBound(validAdhState(~isnan(validAdhState) & validAdhState>0)),'UniformOutput',false);

% For NAs maturing to FCs
arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colorMatureNA),tracksNA(indMature & idCurrent),'UniformOutput',false);
htrackG8 = arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','MarkerSize',markerSize,'Color',colorMatureNA),...
    tracksNA(indMature & idCurrent),'UniformOutput',false);
validAdhState = refineFAID(indMature & idCurrent,iFrameInterest); %cellfun(@(x) x(iFrame),refineFAID(validState));
cellfun(@(x) plot(x(:,2),x(:,1),'Color',colorMatureNA),adhBound(validAdhState(~isnan(validAdhState) & validAdhState>0)),'UniformOutput',false);

% For FCs maturing to FAs
arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colorMatureFC),tracksNA(indMatureFCtoFA & idCurrent),'UniformOutput',false);
htrackG9 = arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','MarkerSize',markerSize,'Color',colorMatureFC),...
    tracksNA(indMatureFCtoFA & idCurrent),'UniformOutput',false);
validAdhState = refineFAID(indMatureFCtoFA & idCurrent,iFrameInterest); %cellfun(@(x) x(iFrame),refineFAID(validState));
htrackSegLine = cellfun(@(x) plot(x(:,2),x(:,1),'Color',colorMatureFC),adhBound(validAdhState(~isnan(validAdhState) & validAdhState>0)),'UniformOutput',false);

%% Labels
classDescription={'NAs maturing', 'FCs maturing to FAs','NAs failing','FCs failing',...
    'FAs failing','NAs not-maturing','FCs not-maturing','FAs stable'};%,
%     'NAs maturing to FCs'};
lgdHandle=legend([htrackG1{1} htrackG9{1} htrackG2{1} htrackG3{1} htrackG4{1} htrackG5{1} htrackG6{1} htrackG7{1}],... %htrackG8{1}],...
    classDescription,'TextColor','k','Location','best');
lgdHandle.Color='k'; lgdHandle.TextColor='w';

%% Save the figure
savefig(mapFig, [MD.getPath filesep 'FocalAdhesionPackage' filesep 'maturingNAsFCs.fig'],'compact')

%% Show time series
% Pick one good one - from maturing ones
[idMaNA,groupIDsAll]=pickAdhesionTracksInteractive(tracksNA(indMatureNAtoFA), imgStack, 'movieData',MD,'tMap',tMap);


