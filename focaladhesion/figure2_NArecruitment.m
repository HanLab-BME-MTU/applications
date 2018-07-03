%% Folder Setting
% box plot for initial rise time lag
% First, gather all the data
paxFolder{1} = '/storage/disk3/NA_recruitment/analysis/Sangyoon/NA_recruitment/Alexia/2015-01-28/Fibro3T3_paxillin/Colocalization/analysis1';
paxFolder{2} = '/storage/disk3/NA_recruitment/analysis/Sangyoon/NA_recruitment/Alexia/2015-07-17/Paxillin1/Colocalization/analysis2';
paxFolder{3} = '/storage/disk3/NA_recruitment/analysis/Sangyoon/NA_recruitment/Alexia/2015-07-24/Paxillin1/Colocalization/analysis2';
paxFolder{4} = '/storage/disk3/NA_recruitment/analysis/Sangyoon/NA_recruitment/Alexia/2015-07-24/Paxillin2/Colocalization/analysis2';
% paxFolder{5} =
% '/storage/disk3/NA_recruitment/analysis/Sangyoon/NA_recruitment/Alexia/2014-07-Alexia/5kPa/ROI/Colocalization/AveragedPaxillin2';
% % Potentially...

vinFolder{1} = '/storage/disk3/NA_recruitment/analysis/Sangyoon/NA_recruitment/Alexia/2015-01-28/CHOK1-vinculin/ROI/Colocalization/analysis2';
vinFolder{2} = '/storage/disk3/NA_recruitment/analysis/Sangyoon/NA_recruitment/Alexia/31-01-2015/Vinc2/Colocalization/analysis1';
vinFolder{3} = '/storage/disk3/NA_recruitment/analysis/Sangyoon/NA_recruitment/Alexia/31-01-2015/Vinc3/Colocalization/analysis2';
vinFolder{4} = '/storage/disk3/NA_recruitment/analysis/Sangyoon/NA_recruitment/Alexia/2015-04-29/cell2/Int1sec/Colocalization/analysis3';
vinFolder{5} = '/storage/disk3/NA_recruitment/analysis/Sangyoon/NA_recruitment/Alexia/2015-07-17/Vinculin5/Colocalization/analysis1';
vinFolder{6} = '/storage/disk3/NA_recruitment/analysis/Sangyoon/NA_recruitment/Alexia/2015-07-24/Vinculin1/Colocalization/analysis2';

talFolder{1} = '/storage/disk3/NA_recruitment/analysis/Sangyoon/NA_recruitment/Alexia/2015-07-24/Talin1/Colocalization/analysisBleachCorrected';
% talFolder{2} = '/storage/disk2/Kevin/2017-02/Talin1_ChoK1_5kPa_1520_01/Colocalization/analysis2'; % Too early spread cell, no FAs
% talFolder{3} = '/storage/disk2/Kevin/2017-02/Talin1_ChoK1_5kPa_1520_03/Colocalization/analysis2'; % Not protruding much
talFolder{2} = '/storage/disk2/Kevin/2017-02/Talin1_ChoK1_5kPa_1520_05_2017_02_09/Colocalization/analysis2'; % Good protrusion, but no G2 observed
% talFolder{5} = '/storage/disk2/Kevin/2017-02/Talin1_ChoK1_5kPa_1520_06/Colocalization/analysis2'; % No active protrusion. Cell is relaxing...
% talFolder{6} = '/storage/disk2/Kevin/2017-02/Talin1_ChoK1_5kPa_1520_06_2017_02_09/Colocalization/analysis2'; % Very little protrusion, no G2
% talFolder{7} = '/storage/disk2/Kevin/2017-02/Talin1_ChoK1_5kPa_1520_07/Colocalization/analysis2'; % Also mostly relaxing
% talFolder{8} = '/storage/disk2/Kevin/2017-02/Talin1_ChoK1_5kPa_1520_07_2017_02_09/Colocalization/analysis2'; % Not healthy looking
talFolder{3} = '/storage/disk2/Kevin/2017-02/Talin1_ChoK1_5kPa_1520_09_2017_02_09/Colocalization/analysis2'; % Good protrusion, no G2
% talFolder{10} = '/storage/disk2/Kevin/2017-02/Talin1_ChoK1_5kPa_1520_08_2017_02_09/Colocalization/analysis2'; % Not healthy looking
% talFolder{4} = '/storage/disk2/Kevin/2017-05-17/ChoK1_TalinShRNA_TalinWT_20170517_1520_008'; % One successful G2
talFolder{4} = '/storage/disk2/Kevin/2017-06-29/ChoK1_shRNA_WT_Rescue_FACS_5kPa_006/Colocalization/analysis1'; % Many G2, almost no G1
talFolder{5} = '/storage/disk2/Kevin/2017-06-29/ChoK1_shRNA_WT_Rescue_FACS_5kPa_011/Colocalization/analysis1'; % Many G1, cloud-like
% talFolder{5} = '/storage/disk2/Kevin/2017-06-29/ChoK1_shRNA_WT_Rescue_FACS_5kPa_012/Colocalization/analysis1'; % Many G1, nice protrusion, not good SNR

% talFolder{2} = '/project/bioinformatics/Danuser_lab/P01adhesion/analysis/Sangyoon/NA_RecruitmentProject/Kevin/2015-07-22_16kPa/Cell2/Colocalization/analysis3';
% talFolder{3} = '/project/bioinformatics/Danuser_lab/P01adhesion/analysis/Sangyoon/NA_RecruitmentProject/Kevin/2015-07-22_16kPa/Cell4/Colocalization/analysis3';
% talFolder{4} = '/project/bioinformatics/Danuser_lab/P01adhesion/analysis/Sangyoon/NA_RecruitmentProject/Kevin/2015-07-22_16kPa/Cell5/Colocalization/analysis3';
% talFolder{5} = '/project/bioinformatics/Danuser_lab/P01adhesion/analysis/Sangyoon/NA_RecruitmentProject/Kevin/2015-07-22_16kPa/Cell8/Colocalization/analysis3';
% talFolder{6} = '/project/bioinformatics/Danuser_lab/P01adhesion/analysis/Sangyoon/NA_RecruitmentProject/Kevin/2015-07-22_16kPa/Cell10/Colocalization/analysis2';
% talFolder{7} = '/project/bioinformatics/Danuser_lab/P01adhesion/analysis/Sangyoon/NA_RecruitmentProject/Kevin/2015-07-22_16kPa/Cell6/Colocalization/analysis3';
numTalMovies=numel(talFolder);
numVinMovies=numel(vinFolder);
numPaxMovies=numel(paxFolder);
%% Do some post-processing - temporary
% progressText(0,'Pax Post-processing')
% for ii=1:4
%     tempPostAnalysisAfterAdhesionClassification(paxFolder{ii})
%     progressText(ii/4,'Pax Post-processing')
% end
% progressText(0,'Vin Post-processing')
% for ii=1:6
%     tempPostAnalysisAfterAdhesionClassification(vinFolder{ii})
%     progressText(ii/6,'Vin Post-processing')
% end
% progressText(0,'Tal Post-processing')
% for ii=1:1
%     tempPostAnalysisAfterAdhesionClassification(talFolder{ii})
%     progressText(ii/1,'Tal Post-processing')
% end
%% Fig 2 High-throughput collection of NA tracks with classification
% Showing overall look over each protein
iRepTal=1; iRepVin=2; iRepPax=2; % or iRepTal=4;
close all
f1path = '/home/sjhan/Documents/Manuscript/NascentAdhesionForce/Figures/Fig1';
if ~exist(f1path,'dir')
    mkdir(f1path)
end
% positioning
pos = get(0, 'DefaultFigurePosition');
% Get all MDs
MDpathTal = fileparts(fileparts(talFolder{iRepTal}));
MDtal = MovieData.load([MDpathTal filesep 'movieData.mat']);
MDpathVin = fileparts(fileparts(vinFolder{iRepVin}));
MDvin = MovieData.load([MDpathVin filesep 'movieData.mat']);
MDpathPax = fileparts(fileparts(paxFolder{iRepPax}));
MDpax = MovieData.load([MDpathPax filesep 'movieData.mat']);
unifiedWidthNanometer = 450*MDtal.pixelSize_;
[minPixSize,iMinPixSize] = min([MDtal.pixelSize_ MDvin.pixelSize_ MDpax.pixelSize_]);
% It is Vinculin that has the smallest pixel size, but it seems strange --
% celll in vinculin exp looks so small compared to the cell in talin exp.
%% Fig 2 time-course in G1 adhesions and time-lag quantification and Cross-variance calculation results
%% Fig 2 a------------------------ Talin--------data loading------------------------------------------
talForceStack = load([talFolder{iRepTal} filesep 'fMap' filesep 'tMap.mat'],'tMap');
talForceStack = talForceStack.tMap;
talImgStack = load([talFolder{iRepTal} filesep 'pax' filesep 'paxImgStack.mat'],'paxImgStack');
talImgStack = talImgStack.paxImgStack;
forceFieldTal = load([MDpathTal filesep 'TFMPackage' filesep 'forceField' filesep 'forceField.mat']);
forceFieldTal = forceFieldTal.forceField;
displFieldTal = load([MDpathTal filesep 'TFMPackage' filesep 'correctedDisplacementField' filesep 'displField.mat']);
displFieldTal = displFieldTal.displField;
%% 2a G1 with adhesion channel and force channel
pos(3:4) = [600 600];
unifiedWidthNanometer = 300*MDtal.pixelSize_;
widthPixelTalin = unifiedWidthNanometer/MDtal.pixelSize_;
CurrentFrameTal=67;

figure('Position', pos, 'PaperPositionMode', 'auto','Color','w');
axes('Position', [0 5/6 1/6 1/6]); avgWidth=0; % avgWidth
imshow(imcomplement(mean(talImgStack(:,:,CurrentFrameTal-avgWidth:CurrentFrameTal+avgWidth),3)),[-900 -100])
% G1 tracks
if ~exist('tracksNAtal','var')
    tic
    tracksNAtal = load([talFolder{iRepTal} filesep 'data' filesep 'tracksNA.mat'],'tracksNA');
    tracksNAtal = tracksNAtal.tracksNA; toc
end
if ~exist('idGroupsTal','var')
    idGroupsTal1struct = load([talFolder{iRepTal} filesep 'data' filesep 'idsClassified.mat'],'idGroup1','idGroup2','idGroup3',...
    'idGroup4','idGroup5','idGroup6','idGroup7','idGroup8','idGroup9');
    idGroupsTal{1} = idGroupsTal1struct.idGroup1;
    idGroupsTal{2} = idGroupsTal1struct.idGroup2;
    idGroupsTal{3} = idGroupsTal1struct.idGroup3;
    idGroupsTal{4} = idGroupsTal1struct.idGroup4;
    idGroupsTal{5} = idGroupsTal1struct.idGroup5;
    idGroupsTal{6} = idGroupsTal1struct.idGroup6;
    idGroupsTal{7} = idGroupsTal1struct.idGroup7;
    idGroupsTal{8} = idGroupsTal1struct.idGroup8;
    idGroupsTal{9} = idGroupsTal1struct.idGroup9;
end

% % CurrentFrame = 197;
% idCurrent = arrayfun(@(x) x.startingFrameExtra<=CurrentFrameTal & x.endingFrameExtra>=CurrentFrameTal,tracksNAtal);
% markerSize=2;
% % plot(arrayfun(@(x) x.xCoord(CurrentFrame),tracksNAtal(idCurrent & (idGroupsTal{1} | idGroupsTal{7}))),arrayfun(@(x) x.yCoord(CurrentFrame),tracksNAtal(idCurrent & (idGroupsTal{1} | idGroupsTal{7}))),'go','MarkerSize',markerSize)
% plot(arrayfun(@(x) x.xCoord(CurrentFrameTal),tracksNAtal(idCurrent & (idGroupsTal{1}))),arrayfun(@(x) x.yCoord(CurrentFrameTal),tracksNAtal(idCurrent & (idGroupsTal{1}))),'go','MarkerSize',markerSize)
% xmat = cell2mat(arrayfun(@(x) x.xCoord(1:CurrentFrameTal),tracksNAtal(idCurrent & idGroupsTal{1}),'UniformOutput',false));
% ymat = cell2mat(arrayfun(@(x) x.yCoord(1:CurrentFrameTal),tracksNAtal(idCurrent & idGroupsTal{1}),'UniformOutput',false));
% plot(xmat',ymat','g','linewidth',0.25)
% plot(arrayfun(@(x) x.xCoord(CurrentFrameTal),tracksNAtal(idCurrent & (idGroupsTal{2}))),arrayfun(@(x) x.yCoord(CurrentFrameTal),tracksNAtal(idCurrent & (idGroupsTal{2}))),'bo','MarkerSize',markerSize)
% 
% axes('Position', [0 4/6 1/6 1/6]); % 
% imshow(mean(talForceStack(:,:,CurrentFrameTal-1:CurrentFrameTal+1),3),[20 400]), colormap(gca,mycmap),hold on
% set(gca,'XLim',[64 64+widthPixelTalin],'YLim',[200 200+widthPixelTalin])
% % plot(arrayfun(@(x) x.xCoord(CurrentFrameTal),tracksNAtal(idCurrent & (idGroupsTal{1} | idGroupsTal{7}))),arrayfun(@(x) x.yCoord(CurrentFrameTal),tracksNAtal(idCurrent & (idGroupsTal{1} | idGroupsTal{7}))),'go','MarkerSize',markerSize)
% plot(arrayfun(@(x) x.xCoord(CurrentFrameTal),tracksNAtal(idCurrent & (idGroupsTal{1}))),arrayfun(@(x) x.yCoord(CurrentFrameTal),tracksNAtal(idCurrent & (idGroupsTal{1}))),'go','MarkerSize',markerSize)
% plot(xmat',ymat','g','linewidth',0.25)
% plot(arrayfun(@(x) x.xCoord(CurrentFrameTal),tracksNAtal(idCurrent & (idGroupsTal{2}))),arrayfun(@(x) x.yCoord(CurrentFrameTal),tracksNAtal(idCurrent & (idGroupsTal{2}))),'wo','MarkerSize',markerSize)
% %% ------------------------ Talin---------showing-----------------------------------------
% widthPixelTalin = unifiedWidthNanometer/MDtal.pixelSize_;
% avgWidth=1;
% CurrentFrameTal=300;
% CurrentFrameTal = min(size(talImgStack,3)-avgWidth,CurrentFrameTal-avgWidth);
% pos(3:4) = [600 600];
% figure('Position', pos, 'PaperPositionMode', 'auto','Color','w');
% % subplot(2,1,1)
% % force
% 
% axes('Position', [30/pos(3) 240/pos(4) 160/pos(3) 160/pos(4)]); % talin
% % axes('Position', [0.03 .35 .3 .30]); % talin
% % imshow(imcomplement(mean(talImgStack(:,:,195:200),3)),[-900 -100])
% imshow(imcomplement(mean(talImgStack(:,:,CurrentFrameTal-avgWidth:CurrentFrameTal+avgWidth),3)),[-800 -100])
% set(gca,'XLim',[30 30+widthPixelTalin],'YLim',[30 30+widthPixelTalin])
% % We want to keep the same zoom scale for all three images
% % set(gca,'XLim',[64 64+widthPixelTalin],'YLim',[120 120+widthPixelTalin])
% hold on
% line([30+20 30+20+round(5000/MDtal.pixelSize_)],[30+widthPixelTalin-20 30+widthPixelTalin-20],'LineWidth',2,'Color','k')
% % % line([64+10 64+10+round(5000/MDtal.pixelSize_)],[120+widthPixelTalin-20 120+widthPixelTalin-20],'LineWidth',2,'Color','k')
% text(30+20, 30+widthPixelTalin-20-30,'5 um','Color','k','Fontsize',7)
%% Force-transmitting-normalized - talin
forceTransmittingTal = arrayfun(@(x) ~isempty(x.forceTransmitting) & x.forceTransmitting,tracksNAtal,'UniformOutput',false);
% forceTransmittingTal = arrayfun(@(x) ~isempty(x.forceTransmitting) & x.forceTransmitting,tracksG1,'UniformOutput',false);
idEmptyFTtal=cellfun(@isempty,forceTransmittingTal);
forceTransmittingTal(idEmptyFTtal)={false};
forceTransmittingTal = cell2mat(forceTransmittingTal);
% axes('Position', [0.21 .45 .11 .12]); % only force-transmitting
% axes('Position', [530/600 530/600 65/600 65/600]); % only force-transmitting
idTalG1 = find((idGroupsTal{1}) & forceTransmittingTal);
lifeTimeTalG1=arrayfun(@(x) x.lifeTime, tracksNAtal(idTalG1));

% presentTalG1=arrayfun(@(x) x.presence(CurrentFrameTal), tracksNAtal(idTalG1));
% midLT_talG1=quantile(lifeTimeTalG1,0.9);
% % find the id whose liftime is cloeset to midLT)talG1
% [LTTalG1sorted,idLTTalG1sorted] = sort(abs(lifeTimeTalG1-midLT_talG1));
% presentTalG1sorted = presentTalG1(idLTTalG1sorted);
% repID_talG1 = idLTTalG1sorted(find(presentTalG1sorted,3));

midLT_talG1=quantile(lifeTimeTalG1,0.9);
repID_talG1 = find(lifeTimeTalG1==midLT_talG1,1);
curTrack = tracksNAtal(idTalG1(repID_talG1)); %3040
% Show this first on top of the image
sF  = curTrack.startingFrameExtraExtra; eF = curTrack.endingFrameExtraExtra;
% hold on,plot(curTrack.xCoord(sF:CurrentFrameTal),curTrack.yCoord(sF:CurrentFrameTal),'g')
arrowWidth=6;
% hold on,arrow3([curTrack.xCoord(CurrentFrameTal)+arrowWidth,curTrack.yCoord(CurrentFrameTal)]+1.5*arrowWidth,...
%     [curTrack.xCoord(CurrentFrameTal),curTrack.yCoord(CurrentFrameTal)],'g0',arrowWidth,1.4*arrowWidth)
hold on,plot(curTrack.xCoord(CurrentFrameTal),curTrack.yCoord(CurrentFrameTal),'go','MarkerSize',2*arrowWidth,'LineWidth',2)
widthPixelTal=60;
% set(gca,'XLim',[234 234+widthPixelTalin],'YLim',[402 402+widthPixelTalin])
% line([234+5 234+5+round(1000/MDtal.pixelSize_)],[402+widthPixelTalin-5 402+widthPixelTalin-5],'LineWidth',2,'Color','k')
% text(234+5, 402+widthPixelTalin-5-5,'1 um','Color','k','Fontsize',7)
set(gca,'XLim',[curTrack.xCoord(CurrentFrameTal)-widthPixelTal/2 curTrack.xCoord(CurrentFrameTal)+widthPixelTal/2],...
    'YLim',[curTrack.yCoord(CurrentFrameTal)-widthPixelTal/2 curTrack.yCoord(CurrentFrameTal)+widthPixelTal/2])
line([curTrack.xCoord(CurrentFrameTal)-widthPixelTal/2+5 curTrack.xCoord(CurrentFrameTal)-widthPixelTal/2+5+round(1000/MDvin.pixelSize_)],...
    [curTrack.yCoord(CurrentFrameTal)+widthPixelTal/2-5 curTrack.yCoord(CurrentFrameTal)+widthPixelTal/2-5],'LineWidth',2,'Color','k')
text(curTrack.xCoord(CurrentFrameTal)-widthPixelTal/2+5, curTrack.yCoord(CurrentFrameTal)+widthPixelTal/2-5-5,'1 um','Color','k','Fontsize',7)
text(curTrack.xCoord(CurrentFrameTal)+widthPixelTal/2-13, curTrack.yCoord(CurrentFrameTal)+widthPixelTal/2-5,...
    [num2str(round((CurrentFrameTal-sF)*MDtal.timeInterval_)) '"'],'Color','k','Fontsize',7)

% Traction map
axes('Position', [1/6 5/6 1/6 1/6]); % talin
fMaxTal=300; avgWidth=0;
imshow(mean(talForceStack(:,:,CurrentFrameTal-avgWidth:CurrentFrameTal+avgWidth),3),[0 fMaxTal])
load('/storage/disk3/NA_recruitment/analysis/Sangyoon/NA_recruitment/MyColormaps','mycmap')
colormap(gca,mycmap); 
hold on
overlayTractionVectors(forceFieldTal(CurrentFrameTal),displFieldTal(CurrentFrameTal),0.1,3,'Color',[0.7 0.7 0.7],'ShiftField',false);
hc=colorbar('east'); hc.Position(3)=0.01; hc.Ticks=[0 fMaxTal];
% set(gca,'XLim',[234 234+widthPixelTalin],'YLim',[402 402+widthPixelTalin])
set(gca,'XLim',[curTrack.xCoord(CurrentFrameTal)-widthPixelTal/2 curTrack.xCoord(CurrentFrameTal)+widthPixelTal/2],...
    'YLim',[curTrack.yCoord(CurrentFrameTal)-widthPixelTal/2 curTrack.yCoord(CurrentFrameTal)+widthPixelTal/2])
hold on
plot(curTrack.xCoord(CurrentFrameTal),curTrack.yCoord(CurrentFrameTal),'go','MarkerSize',2*arrowWidth,'LineWidth',2)
line([234+5 234+5+round(1000/MDtal.pixelSize_)],[402+widthPixelTalin-5 402+widthPixelTalin-5],'LineWidth',2,'Color','w')
title('Talin')

tRange = ((sF:eF) - sF)*MDtal.timeInterval_;
axes('Position', [40/600 430/600 145/600 65/600]); % 
plot(tRange, curTrack.ampTotal(sF:eF),'k-')
xlabel('Time (sec)'); ylabel('F.I. (a.u.)')
set(gca,'FontSize',7)
xlim([tRange(1),tRange(end)]); ylim auto; 
hold on
yl = ylim; line([curTrack.firstIncreseTimeInt-sF*MDtal.timeInterval_ curTrack.firstIncreseTimeInt-sF*MDtal.timeInterval_],...
    [yl(1) yl(2)],'Color',[0.5 .5 .5],'LineStyle','--')
% plotIntensityForce(tracksNAtal((idGroupsTal{1}) & forceTransmittingTal),[],false,false,'UseCurrentAxis',true,...
%     'Source',{'ampTotal'},'plotCohorts',true,'tInterval',MDtal.timeInterval_,'prePostFrames',10,'yNormalization',false)
% plotIntensityForce(tracksNAtal((idGroupsTal{1} | idGroupsTal{7}) & forceTransmittingTal),[],false,false,'UseCurrentAxis',true,...
%     'Source',{'ampTotal'},'plotCohorts',true,'tInterval',MDtal.timeInterval_,'prePostFrames',10,'yNormalization',true)
% xlim([-10 150])
% ylim([251 460])
axes('Position', [40/600 330/600 145/600 65/600]); % 
plot(tRange, curTrack.forceMag(sF:eF),'k-')
xlabel('Time (sec)'); ylabel('Traction (Pa)')
set(gca,'FontSize',7)
xlim([tRange(1),tRange(end)]); ylim auto; 
hold on
yl = ylim; line([curTrack.firstIncreseTimeForce-sF*MDtal.timeInterval_ curTrack.firstIncreseTimeForce-sF*MDtal.timeInterval_],...
    [yl(1) yl(2)],'Color',[0.5 .5 .5],'LineStyle','--')

% axes('Position', [130/600 430/600 65/600 65/600]); % only force-transmitting
% plotIntensityForce(tracksNAtal((idGroupsTal{1}) & forceTransmittingTal),[],false,false,'UseCurrentAxis',true,...
%     'Source',{'forceMag'},'plotCohorts',true,'tInterval',MDtal.timeInterval_,'prePostFrames',10,'yNormalization',false)
% % plotIntensityForce(tracksNAtal((idGroupsTal{1} | idGroupsTal{7}) & forceTransmittingTal),[],false,false,'UseCurrentAxis',true,...
% %     'Source',{'forceMag'},'plotCohorts',true,'tInterval',MDtal.timeInterval_,'prePostFrames',10,'yNormalization',true)
% xlim([-10 150])
% ylim([20 100])

% Bcc calculation and display
curBcc = crossVariance(tracksNAtal(idTalG1(repID_talG1)).ampTotal,tracksNAtal(idTalG1(repID_talG1)).forceMag,11);
axes('Position', [40/600 230/600 145/600 65/600]); % 
plot(tRange,curBcc(sF:eF),'k'); xlabel('Time (sec)'); ylabel('Cross Variance, Bcc')
set(gca,'FontSize',7)
xlim([tRange(1),tRange(end)]); ylim auto 
hold on
yl = ylim; line([(curTrack.startingFrameExtra-sF)*MDtal.timeInterval_ (curTrack.startingFrameExtra-sF)*MDtal.timeInterval_],...
    [yl(1) yl(2)],'Color',[0.1 .2 .8],'LineStyle','--')
[~,peakBcc] = nanmax(curBcc);
line([(peakBcc-sF)*MDtal.timeInterval_ (peakBcc-sF)*MDtal.timeInterval_],...
    [yl(1) yl(2)],'Color',[0.8 .6 .2],'LineStyle','--')
%% vinculin loading
clear talForceStack forceFieldTal talImgStack 
vinForceStack = load([vinFolder{iRepVin} filesep 'fMap' filesep 'tMap.mat'],'tMap');
vinForceStack = vinForceStack.tMap;
vinImgStack = load([vinFolder{iRepVin} filesep 'pax' filesep 'paxImgStack.mat'],'paxImgStack');
vinImgStack = vinImgStack.paxImgStack;
forceFieldVin = load([MDpathVin filesep 'TFMPackage' filesep 'forceField' filesep 'forceField.mat']);
forceFieldVin = forceFieldVin.forceField;
displFieldVin = load([MDpathVin filesep 'TFMPackage' filesep 'correctedDisplacementField' filesep 'displField.mat']);
displFieldVin = displFieldVin.displField;

%% ------------------------ Vinculin--------------------------------------------------
% axes('Position', [0.35 .65 .3 .35]); % vinculin
% subplot(2,1,1)
MDvin5.pixelSize_=108;
% force
% widthPixelVinculin = unifiedWidthNanometer/MDvin5.pixelSize_;
% Vin image
axes('Position', [1/3 5/6 1/6 1/6]); % 
CurrentFrameVin = 100;
imshow(imcomplement(mean(vinImgStack(:,:,CurrentFrameVin-2:CurrentFrameVin+2),3)),[-29000 -3000])

if ~exist('tracksNAvin','var')
    tic
    tracksNAvin = load([vinFolder{iRepVin} filesep 'data' filesep 'tracksNA.mat'],'tracksNA');
    tracksNAvin = tracksNAvin.tracksNA; toc
end
% Example of intensity, force, edge position and adhesion movement in g1
if ~exist('idGroupsVin5struct','var')
    idGroupsVin5struct = load([vinFolder{iRepVin} filesep 'data' filesep 'idsClassified.mat'],'idGroup1filtered','idGroup2filtered','idGroup3filtered',...
    'idGroup4filtered','idGroup5filtered','idGroup6filtered','idGroup7filtered','idGroup8filtered','idGroup9filtered',...
    'idGroup1','idGroup2','idGroup3',...
    'idGroup4','idGroup5','idGroup6','idGroup7','idGroup8','idGroup9');
    idGroupsVin5{1} = idGroupsVin5struct.idGroup1;
    idGroupsVin5{2} = idGroupsVin5struct.idGroup2;
    idGroupsVin5{3} = idGroupsVin5struct.idGroup3;
    idGroupsVin5{4} = idGroupsVin5struct.idGroup4;
    idGroupsVin5{5} = idGroupsVin5struct.idGroup5;
    idGroupsVin5{6} = idGroupsVin5struct.idGroup6;
    idGroupsVin5{7} = idGroupsVin5struct.idGroup7;
    idGroupsVin5{8} = idGroupsVin5struct.idGroup8;
    idGroupsVin5{9} = idGroupsVin5struct.idGroup9;
%     end
end

forceTransmittingVin = arrayfun(@(x) ~isempty(x.forceTransmitting) & x.forceTransmitting,tracksNAvin,'UniformOutput',false);
idEmptyFTvin=cellfun(@isempty,forceTransmittingVin);
forceTransmittingVin(idEmptyFTvin)={false};
forceTransmittingVin = cell2mat(forceTransmittingVin);

idVinG1 = find((idGroupsVin5{1}) & forceTransmittingVin);
lifeTimeVinG1=arrayfun(@(x) x.lifeTime, tracksNAvin(idVinG1));
presentVinG1=arrayfun(@(x) x.presence(CurrentFrameVin), tracksNAvin(idVinG1));
midLT_vinG1=quantile(lifeTimeVinG1,0.9);
% find the id whose liftime is cloeset to midLT)vinG1
[LTVinG1sorted,idLTVinG1sorted] = sort(abs(lifeTimeVinG1-midLT_vinG1));
presentVinG1sorted = presentVinG1(idLTVinG1sorted);
repID_vinG1 = idLTVinG1sorted(find(presentVinG1sorted,3));
repID_vinG1 = repID_vinG1(1);
curTrack = tracksNAvin(idVinG1(repID_vinG1));

sF  = curTrack.startingFrameExtraExtra; eF = curTrack.endingFrameExtraExtra;
arrowWidth=6;
hold on,plot(curTrack.xCoord(CurrentFrameVin),curTrack.yCoord(CurrentFrameVin),'go','MarkerSize',2*arrowWidth,'LineWidth',2)
widthPixelVin=60;
set(gca,'XLim',[curTrack.xCoord(CurrentFrameVin)-widthPixelVin/2 curTrack.xCoord(CurrentFrameVin)+widthPixelVin/2],...
    'YLim',[curTrack.yCoord(CurrentFrameVin)-widthPixelVin/2 curTrack.yCoord(CurrentFrameVin)+widthPixelVin/2])
line([curTrack.xCoord(CurrentFrameVin)-widthPixelVin/2+5 curTrack.xCoord(CurrentFrameVin)-widthPixelVin/2+5+round(1000/MDvin.pixelSize_)],...
    [curTrack.yCoord(CurrentFrameVin)+widthPixelVin/2-5 curTrack.yCoord(CurrentFrameVin)+widthPixelVin/2-5],'LineWidth',2,'Color','k')
text(curTrack.xCoord(CurrentFrameVin)-widthPixelVin/2+5, curTrack.yCoord(CurrentFrameVin)+widthPixelVin/2-5-5,'1 um','Color','k','Fontsize',7)
text(curTrack.xCoord(CurrentFrameVin)+widthPixelVin/2-13, curTrack.yCoord(CurrentFrameVin)+widthPixelVin/2-5,...
    [num2str(round((CurrentFrameVin-sF)*MDvin.timeInterval_)) '"'],'Color','k','Fontsize',7)

% Traction map
axes('Position', [3/6 5/6 1/6 1/6]); % talin
fMaxVin=300; avgWidth=1;
imshow(mean(vinForceStack(:,:,CurrentFrameVin-avgWidth:CurrentFrameVin+avgWidth),3),[0 fMaxVin])
colormap(gca,mycmap); 
hold on
overlayTractionVectors(forceFieldVin(CurrentFrameVin),forceFieldVin(CurrentFrameVin),0.1,3,'Color',[0.7 0.7 0.7],'ShiftField',false);
hc=colorbar('east'); hc.Position(3)=0.01; hc.Ticks=[0 fMaxVin];
set(gca,'XLim',[curTrack.xCoord(CurrentFrameVin)-widthPixelVin/2 curTrack.xCoord(CurrentFrameVin)+widthPixelVin/2],...
    'YLim',[curTrack.yCoord(CurrentFrameVin)-widthPixelVin/2 curTrack.yCoord(CurrentFrameVin)+widthPixelVin/2])
hold on
plot(curTrack.xCoord(CurrentFrameVin),curTrack.yCoord(CurrentFrameVin),'go','MarkerSize',2*arrowWidth,'LineWidth',2)
line([curTrack.xCoord(CurrentFrameVin)-widthPixelVin/2+5 curTrack.xCoord(CurrentFrameVin)-widthPixelVin/2+5+round(1000/MDvin.pixelSize_)],...
    [curTrack.yCoord(CurrentFrameVin)+widthPixelVin/2-5 curTrack.yCoord(CurrentFrameVin)+widthPixelVin/2-5],'LineWidth',2,'Color','w')

tRange = ((sF:eF) - sF)*MDvin.timeInterval_;
axes('Position', [240/600 430/600 145/600 65/600]); % 
plot(tRange, curTrack.ampTotal(sF:eF),'k-')
xlabel('Time (sec)'); ylabel('F.I. (a.u.)')
set(gca,'FontSize',7)
xlim([tRange(1),tRange(end)]); ylim auto; 
hold on
%firstIncreseTimeInt Vin
yl = ylim; line([curTrack.firstIncreseTimeInt-sF*MDvin.timeInterval_ curTrack.firstIncreseTimeInt-sF*MDvin.timeInterval_],...
    [yl(1) yl(2)],'Color',[0.5 .5 .5],'LineStyle','--')

axes('Position', [240/600 330/600 145/600 65/600]); % 
plot(tRange, curTrack.forceMag(sF:eF),'k-')
xlabel('Time (sec)'); ylabel('Traction (Pa)')
set(gca,'FontSize',7)
xlim([tRange(1),tRange(end)]); ylim auto; 
hold on
yl = ylim; line([curTrack.firstIncreseTimeForce-sF*MDvin.timeInterval_ curTrack.firstIncreseTimeForce-sF*MDvin.timeInterval_],...
    [yl(1) yl(2)],'Color',[0.5 .5 .5],'LineStyle','--')

% Bcc calculation and display
curBcc = crossVariance(curTrack.ampTotal,curTrack.forceMag,11);
axes('Position', [240/600 230/600 145/600 65/600]); % 
plot(tRange,curBcc(sF:eF),'k'); xlabel('Time (sec)'); ylabel('Cross Variance, Bcc')
set(gca,'FontSize',7)
xlim([tRange(1),tRange(end)]); ylim auto 
hold on
yl = ylim; line([(curTrack.startingFrameExtra-sF)*MDvin.timeInterval_ (curTrack.startingFrameExtra-sF)*MDvin.timeInterval_],...
    [yl(1) yl(2)],'Color',[0.1 .2 .8],'LineStyle','--')
[~,peakBcc] = nanmax(curBcc);
line([(peakBcc-sF)*MDvin.timeInterval_ (peakBcc-sF)*MDvin.timeInterval_],...
    [yl(1) yl(2)],'Color',[0.8 .6 .2],'LineStyle','--')

% axes('Position', [230/600 330/600 65/600 65/600]); % 
% % plotIntensityForce(tracksNAvin((idGroupsVin5{1} | idGroupsVin5{7}) & forceTransmittingVin),[],false,false,'UseCurrentAxis',true,...
% %     'Source',{'ampTotal'},'plotCohorts',true,'tInterval',MDvin.timeInterval_,'prePostFrames',10,'yNormalization',true)
% plotIntensityForce(tracksNAvin((idGroupsVin5{1}) & forceTransmittingVin),[],false,false,'UseCurrentAxis',true,...
%     'Source',{'ampTotal'},'plotCohorts',true,'tInterval',MDvin.timeInterval_,'prePostFrames',10,'yNormalization',false)
% xlim([-5 220])
% ylim([230 680])
% % axes('Position', [230/600 230/600 65/600 65/600]); % only force-transmitting
% axes('Position', [130/600 230/600 65/600 65/600]); % 
% % plotIntensityForce(tracksNAvin((idGroupsVin5{1} | idGroupsVin5{7}) & forceTransmittingVin),[],false,false,'UseCurrentAxis',true,...
% %     'Source',{'forceMag'},'plotCohorts',true,'tInterval',MDvin.timeInterval_,'prePostFrames',10,'yNormalization',true)
% plotIntensityForce(tracksNAvin((idGroupsVin5{1}) & forceTransmittingVin),[],false,false,'UseCurrentAxis',true,...
%     'Source',{'forceMag'},'plotCohorts',true,'tInterval',MDvin.timeInterval_,'prePostFrames',10,'yNormalization',false)
% xlim([-5 220])
% ylim([5 305])

%% pax loading
paxForceStack = load([paxFolder{iRepPax} filesep 'fMap' filesep 'tMap.mat'],'tMap');
paxForceStack = paxForceStack.tMap;
paxImgStack = load([paxFolder{iRepPax} filesep 'pax' filesep 'paxImgStack.mat'],'paxImgStack');
paxImgStack = paxImgStack.paxImgStack;
forceFieldPax = load([MDpathPax filesep 'TFMPackage' filesep 'forceField' filesep 'forceField.mat']);
forceFieldPax = forceFieldPax.forceField;
displFieldPax = load([MDpathPax filesep 'TFMPackage' filesep 'correctedDisplacementField' filesep 'displField.mat']);
displFieldPax = displFieldPax.displField;
%% Paxillin
% axes('Position', [0.66 .8 .15 .17]); % 
if ~exist('tracksNApax','var')
    tic
    tracksNApax = load([paxFolder{iRepPax} filesep 'data' filesep 'tracksNA.mat'],'tracksNA');
    tracksNApax = tracksNApax.tracksNA; toc
end
if ~exist('idGroupsPax','var')
    idGroupsPaxStruct = load([paxFolder{iRepPax} filesep 'data' filesep 'idsClassified.mat'],'idGroup1','idGroup2','idGroup3',...
    'idGroup4','idGroup5','idGroup6','idGroup7','idGroup8','idGroup9');
    idGroupsPax{1} = idGroupsPaxStruct.idGroup1;
    idGroupsPax{2} = idGroupsPaxStruct.idGroup2;
    idGroupsPax{3} = idGroupsPaxStruct.idGroup3;
    idGroupsPax{4} = idGroupsPaxStruct.idGroup4;
    idGroupsPax{5} = idGroupsPaxStruct.idGroup5;
    idGroupsPax{6} = idGroupsPaxStruct.idGroup6;
    idGroupsPax{7} = idGroupsPaxStruct.idGroup7;
    idGroupsPax{8} = idGroupsPaxStruct.idGroup8;
    idGroupsPax{9} = idGroupsPaxStruct.idGroup9;
end
% Extra investigation
tracksNApax=readIntensityFromTracks(tracksNApax,paxImgStack,1,'extraLength',30,'movieData',MDpax);
tracksNApax=readIntensityFromTracks(tracksNApax,paxForceStack,2,'movieData',MDpax);
[firstIncreseTimeIntAgainstForce,forceTransmittingAll...
    ,firstIncreseTimeIntAll,firstIncreseTimeForce,bkgMaxIntAll,bkgMaxForce] =...
    calculateFirstIncreaseTimeTracks(tracksNApax,.5,10,MDpax.timeInterval_);

for pp=1:numel(tracksNApax)
    tracksNApax(pp).forceTransmitting=forceTransmittingAll(pp);
    tracksNApax(pp).firstIncreseTimeInt=firstIncreseTimeIntAll(pp);
    tracksNApax(pp).firstIncreseTimeForce=firstIncreseTimeForce(pp);
    tracksNApax(pp).firstIncreseTimeIntAgainstForce=firstIncreseTimeIntAgainstForce(pp);
    tracksNApax(pp).bkgMaxInt=bkgMaxIntAll(pp);
    tracksNApax(pp).bkgMaxForce=bkgMaxForce(pp);
end

%% pax displaying
quantilePax=0.56; %figure('Position', pos, 'PaperPositionMode', 'auto','Color','w');

CurrentFramePax = 176;
axes('Position', [4/6 5/6 1/6 1/6]); % 
imshow(imcomplement(mean(paxImgStack(:,:,CurrentFramePax-1:CurrentFramePax+1),3)),[-800 -160])
% widthPixelPaxilln = unifiedWidthNanometer/MDpax.pixelSize_;
% set(gca,'XLim',[30 30+widthPixelPaxilln],'YLim',[30 30+widthPixelPaxilln])
% hold on
% idCurrent = arrayfun(@(x) x.startingFrameExtra<=CurrentFramePax & x.endingFrameExtra>=CurrentFramePax,tracksNApax);
% plot(arrayfun(@(x) x.xCoord(CurrentFramePax),tracksNApax(idCurrent & (idGroupsPax{1}))),arrayfun(@(x) x.yCoord(CurrentFramePax),tracksNApax(idCurrent & (idGroupsPax{1}))),'go','MarkerSize',markerSize)
% plot(arrayfun(@(x) x.xCoord(CurrentFramePax),tracksNApax(idCurrent & (idGroupsPax{2}))),arrayfun(@(x) x.yCoord(CurrentFramePax),tracksNApax(idCurrent & (idGroupsPax{2}))),'bo','MarkerSize',markerSize)
% plot(arrayfun(@(x) x.xCoord(CurrentFrame),tracksNApax(idCurrent & (idGroupsPax{1} |  idGroupsPax{7}))),arrayfun(@(x) x.yCoord(CurrentFrame),tracksNApax(idCurrent & (idGroupsPax{1} |  idGroupsPax{7}))),'go','MarkerSize',markerSize)

forceTransmittingPax = arrayfun(@(x) ~isempty(x.forceTransmitting) & x.forceTransmitting,tracksNApax,'UniformOutput',false);
idEmptyFTpax=cellfun(@isempty,forceTransmittingPax);
forceTransmittingPax(idEmptyFTpax)={false};
forceTransmittingPax = cell2mat(forceTransmittingPax);
% axes('Position', [130/600 130/600 65/600 65/600]); % 
% plotIntensityForce(tracksNApax((idGroupsPax{1}) & forceTransmittingPax),[],false,false,'UseCurrentAxis',true,...
%     'Source',{'ampTotal'},'plotCohorts',true,'tInterval',MDpax.timeInterval_,'prePostFrames',10,'yNormalization',false)
% axes('Position', [230/600 30/600 65/600 65/600]); % only force-transmitting
% xlim([-5 140])
% ylim([300 750])
% axes('Position', [130/600 30/600 65/600 65/600]); % 
% plotIntensityForce(tracksNApax((idGroupsPax{1}) & forceTransmittingPax),[],false,false,'UseCurrentAxis',true,...
%     'Source',{'forceMag'},'plotCohorts',true,'tInterval',MDpax.timeInterval_,'prePostFrames',10,'yNormalization',false)
% xlim([-5 130])
% ylim([10 150])

% Find the long-enough representative track
idPaxG1 = find((idGroupsPax{1}) & forceTransmittingPax);
lifeTimePaxG1=arrayfun(@(x) x.lifeTime, tracksNApax(idPaxG1));
presentPaxG1=arrayfun(@(x) x.presence(CurrentFramePax), tracksNApax(idPaxG1));
midLT_paxG1=quantile(lifeTimePaxG1,quantilePax);
% find the id whose liftime is cloeset to midLT)paxG1
[LTPaxG1sorted,idLTPaxG1sorted] = sort(abs(lifeTimePaxG1-midLT_paxG1));
presentPaxG1sorted = presentPaxG1(idLTPaxG1sorted);
repID_paxG1 = idLTPaxG1sorted(find(presentPaxG1sorted,4));
repID_paxG1 = repID_paxG1(1);
absoluteIDG1 = idPaxG1(repID_paxG1);
% absoluteIDG1 = 182; %754 are G2
curTrack = tracksNApax(absoluteIDG1);
disp(['absoluteIDG1= ' num2str(absoluteIDG1)])    

% Overlaying the track on top of the image
sF  = curTrack.startingFrameExtraExtra; eF = curTrack.endingFrameExtraExtra;
arrowWidth=6;
hold on,plot(curTrack.xCoord(CurrentFramePax),curTrack.yCoord(CurrentFramePax),'go','MarkerSize',2*arrowWidth,'LineWidth',2)
widthPixelPax=60;
set(gca,'XLim',[curTrack.xCoord(CurrentFramePax)-widthPixelPax/2 curTrack.xCoord(CurrentFramePax)+widthPixelPax/2],...
    'YLim',[curTrack.yCoord(CurrentFramePax)-widthPixelPax/2 curTrack.yCoord(CurrentFramePax)+widthPixelPax/2])
line([curTrack.xCoord(CurrentFramePax)-widthPixelPax/2+5 curTrack.xCoord(CurrentFramePax)-widthPixelPax/2+5+round(1000/MDtal.pixelSize_)],...
    [curTrack.yCoord(CurrentFramePax)+widthPixelPax/2-5 curTrack.yCoord(CurrentFramePax)+widthPixelPax/2-5],'LineWidth',2,'Color','k')
text(curTrack.xCoord(CurrentFramePax)-widthPixelPax/2+5, curTrack.yCoord(CurrentFramePax)+widthPixelPax/2-5-5,'1 um','Color','k','Fontsize',7)
text(curTrack.xCoord(CurrentFramePax)+widthPixelPax/2-13, curTrack.yCoord(CurrentFramePax)+widthPixelPax/2-5,...
    [num2str(round((CurrentFramePax-sF)*MDtal.timeInterval_)) '"'],'Color','k','Fontsize',7)

% Traction map
axes('Position', [5/6 5/6 1/6 1/6]); % talin
fMaxPax=300; avgWidth=1;
imshow(mean(paxForceStack(:,:,CurrentFramePax-avgWidth:CurrentFramePax+avgWidth),3),[0 fMaxPax])
colormap(gca,mycmap); 
hold on
overlayTractionVectors(forceFieldPax(CurrentFramePax),forceFieldPax(CurrentFramePax),0.1,3,'Color',[0.7 0.7 0.7],'ShiftField',false);
hc=colorbar('east'); hc.Position(3)=0.01; hc.Ticks=[0 fMaxPax];
set(gca,'XLim',[curTrack.xCoord(CurrentFramePax)-widthPixelPax/2 curTrack.xCoord(CurrentFramePax)+widthPixelPax/2],...
    'YLim',[curTrack.yCoord(CurrentFramePax)-widthPixelPax/2 curTrack.yCoord(CurrentFramePax)+widthPixelPax/2])
hold on
plot(curTrack.xCoord(CurrentFramePax),curTrack.yCoord(CurrentFramePax),'go','MarkerSize',2*arrowWidth,'LineWidth',2)
line([curTrack.xCoord(CurrentFramePax)-widthPixelPax/2+5 curTrack.xCoord(CurrentFramePax)-widthPixelPax/2+5+round(1000/MDpax.pixelSize_)],...
    [curTrack.yCoord(CurrentFramePax)+widthPixelPax/2-5 curTrack.yCoord(CurrentFramePax)+widthPixelPax/2-5],'LineWidth',2,'Color','w')

tRange = ((sF:eF) - sF)*MDpax.timeInterval_;
axes('Position', [440/600 430/600 145/600 65/600]); % 
plot(tRange, curTrack.ampTotal(sF:eF),'k-')
xlabel('Time (sec)'); ylabel('F.I. (a.u.)')
set(gca,'FontSize',7)
xlim([tRange(1),tRange(end)]); ylim auto; 
hold on
%firstIncreseTimeInt Pax
yl = ylim; line([curTrack.firstIncreseTimeInt-sF*MDpax.timeInterval_ curTrack.firstIncreseTimeInt-sF*MDpax.timeInterval_],...
    [yl(1) yl(2)],'Color',[0.5 .5 .5],'LineStyle','--')

axes('Position', [440/600 330/600 145/600 65/600]); % 
plot(tRange, curTrack.forceMag(sF:eF),'k-')
xlabel('Time (sec)'); ylabel('Traction (Pa)')
set(gca,'FontSize',7)
xlim([tRange(1),tRange(end)]); ylim auto; 
hold on
yl = ylim; line([curTrack.firstIncreseTimeForce-sF*MDpax.timeInterval_ curTrack.firstIncreseTimeForce-sF*MDpax.timeInterval_],...
    [yl(1) yl(2)],'Color',[0.5 .5 .5],'LineStyle','--')

% Bcc calculation and display
curBcc = crossVariance(curTrack.ampTotal,curTrack.forceMag,11);
axes('Position', [440/600 230/600 145/600 65/600]); % 
plot(tRange,curBcc(sF:eF),'k'); xlabel('Time (sec)'); ylabel('Cross Variance, Bcc')
set(gca,'FontSize',7)
xlim([tRange(1),tRange(end)]); ylim auto 
hold on
yl = ylim; line([(curTrack.startingFrameExtra-sF)*MDpax.timeInterval_ (curTrack.startingFrameExtra-sF)*MDpax.timeInterval_],...
    [yl(1) yl(2)],'Color',[0.1 .2 .8],'LineStyle','--')
[~,peakBcc] = nanmax(curBcc);
line([(peakBcc-sF)*MDpax.timeInterval_ (peakBcc-sF)*MDpax.timeInterval_],...
    [yl(1) yl(2)],'Color',[0.8 .6 .2],'LineStyle','--')

% % axes('Position', [0.66 .6 .15 .17]); % 
% axes('Position', [0 0 1/6 1/6]); % 
% imshow(mean(paxForceStack(:,:,295:300),3),[20 1000]), colormap(gca,mycmap),hold on
% set(gca,'XLim',[30 30+widthPixelPaxilln],'YLim',[30 30+widthPixelPaxilln])
% plot(arrayfun(@(x) x.xCoord(CurrentFrame),tracksNApax(idCurrent & (idGroupsPax{1}))),arrayfun(@(x) x.yCoord(CurrentFrame),tracksNApax(idCurrent & (idGroupsPax{1}))),'go','MarkerSize',markerSize)
% plot(arrayfun(@(x) x.xCoord(CurrentFrame),tracksNApax(idCurrent & (idGroupsPax{2}))),arrayfun(@(x) x.yCoord(CurrentFrame),tracksNApax(idCurrent & (idGroupsPax{2}))),'wo','MarkerSize',markerSize)
% % plot(arrayfun(@(x) x.xCoord(CurrentFrame),tracksNApax(idCurrent & (idGroupsPax{1} |  idGroupsPax{7}))),arrayfun(@(x) x.yCoord(CurrentFrame),tracksNApax(idCurrent & (idGroupsPax{1} |  idGroupsPax{7}))),'go','MarkerSize',markerSize)
%% Statistics about G1-talin
% Let's go first in naive way
% Talin first
tLagBccPosG1All=cell(numTalMovies,1);
tLagBccNegG1All=cell(numTalMovies,1);
BccInLaterG1All=cell(numTalMovies,1);
BccInPriorG1All=cell(numTalMovies,1);
tLagBccPosG2All=cell(numTalMovies,1);
tLagBccNegG2All=cell(numTalMovies,1);
BccInLaterG2All=cell(numTalMovies,1);
BccInPriorG2All=cell(numTalMovies,1);

totalInitTLagTal=[];
totalPeakTLagTal=[];
splineParamInit=0.95;
preDetecFactor=1.5/10; 
%% running
for ii=1:numTalMovies
    if exist([talFolder{ii} filesep 'data' filesep 'tracksG2real.mat'],'file')
        tracksG2 = load([talFolder{ii} filesep 'data' filesep 'tracksG2real.mat'],'tracksG2');
        curTracksNAtalG2f=tracksG2.tracksG2;
        tracksG1 = load([talFolder{ii} filesep 'data' filesep 'tracksG1real.mat'],'tracksG1');
        curTracksNAtalG1=tracksG1.tracksG1;
        curMDpath = fileparts(fileparts(talFolder{ii}));
        curMD = load([curMDpath filesep 'movieData.mat']);
        curMD = curMD.MD;
    else

        % try to filter G2 based on starting amplitude, peaks and slopes...
        curTracksNAtalG2 = load([talFolder{ii} filesep 'data' filesep 'tracksG2.mat'],'tracksG2');
        curTracksNAtalG2 = curTracksNAtalG2.tracksG2; 
        % Get rid of overlapping tracks
        curTracksNAtalG2 = filterOverlappingTracks(curTracksNAtalG2);
        curForceTransmittingG2 = arrayfun(@(x) ~isempty(x.firstIncreseTimeIntAgainstForce), curTracksNAtalG2);
        curStartingAmpG2 = arrayfun(@(x) x.ampTotal(x.startingFrameExtra), curTracksNAtalG2);

        curTracksNAtalG1 = load([talFolder{ii} filesep 'data' filesep 'tracksG1.mat'],'tracksG1');
        curTracksNAtalG1 = curTracksNAtalG1.tracksG1; 
        curTimeDelay = arrayfun(@(x) x.firstIncreseTimeIntAgainstForce, curTracksNAtalG1,'UniformOutput',false);
        curStartingAmpG1 = arrayfun(@(x) x.ampTotal(x.startingFrameExtra), curTracksNAtalG1);
        thresAmp = mean(curStartingAmpG1)+3*std(curStartingAmpG1);

        % filtering for G2
        stateFA = arrayfun(@(x) any(cellfun(@(y) strcmp(y,'FC') || strcmp(y,'FA'), x.state)), curTracksNAtalG2);
    %         curTimeDelay = arrayfun(@(x) x.firstIncreseTimeIntAgainstForce, curTracksNAtalG2(realG2tal));

        % Go through each adhesion track and make sure
        curMDpath = fileparts(fileparts(talFolder{ii}));
        curMD = load([curMDpath filesep 'movieData.mat']);
        curMD = curMD.MD;
        tMap = load([talFolder{ii} filesep 'fMap' filesep 'tMap.mat'],'tMap');
        tMap = tMap.tMap;
        imgMap = load([talFolder{ii} filesep 'pax' filesep 'paxImgStack.mat'],'paxImgStack');
        imgMap = imgMap.paxImgStack;
        tInterval = curMD.timeInterval_;

        % Get rid of overlapping tracks
        curTracksNAtalG1 = filterOverlappingTracks(curTracksNAtalG1);
        % Get rid of tracks that are obviously inside the cell
        distG1 = arrayfun(@(x) mean(x.distToEdge),curTracksNAtalG1);
        thresDist = 10; % um
        curTracksNAtalG1 = curTracksNAtalG1(distG1<thresDist);

        [firstIncreseTimeIntAgainstForce,forceTransmittingAll...
        ,firstIncreseTimeIntAll,firstIncreseTimeForce,bkgMaxIntAll,bkgMaxForce] =...
            calculateFirstIncreaseTimeTracks(curTracksNAtalG1,splineParamInit,preDetecFactor,tInterval);
        curTracksNAtalG1 = curTracksNAtalG1(forceTransmittingAll);
        firstIncreseTimeIntAll=firstIncreseTimeIntAll(forceTransmittingAll);
        firstIncreseTimeForce=firstIncreseTimeForce(forceTransmittingAll);
        firstIncreseTimeIntAgainstForce=firstIncreseTimeIntAgainstForce(forceTransmittingAll);
        bkgMaxIntAll=bkgMaxIntAll(forceTransmittingAll);
        bkgMaxForce=bkgMaxForce(forceTransmittingAll);
        forceTransmittingAll=forceTransmittingAll(forceTransmittingAll);

        % Further G1 critera
        % 1. Based on amplitude 
        curAmpSlopeGroup = arrayfun(@(x) x.ampSlope, curTracksNAtalG1);
        curEarlyAmpSlopeGroup = arrayfun(@(x) x.earlyAmpSlope, curTracksNAtalG1);
        % We decided to regard amplitude with flat slope as noise.
        %             figure, plot(curAmpSlopeGroup,curEarlyAmpSlopeGroup,'*')
        indFlatAmp = curAmpSlopeGroup<=0 & curEarlyAmpSlopeGroup<=0;
        % 2. Based on forceMag
        curVeryEarlyAmpSlopeGroup = NaN(numel(curTracksNAtalG1),1);
        curForceSlopeGroup = NaN(numel(curTracksNAtalG1),1);
        curForceEarlySlopeGroup = NaN(numel(curTracksNAtalG1),1);
        curAmpOverallSlopeGroup = NaN(numel(curTracksNAtalG1),1);
        curAmpOverallSlopeUpErrGroup = NaN(numel(curTracksNAtalG1),1);
        curAmpOverallSlopeDownErrGroup = NaN(numel(curTracksNAtalG1),1);
        curForceLateSlopeGroup = NaN(numel(curTracksNAtalG1),1);

        periodFrames = 30;
        rr=0;
        for pp=1:numel(curTracksNAtalG1)
            rr=rr+1;
            curTrack = curTracksNAtalG1(pp);
            [~,curForceSlopeGroup(rr)] = regression((1:curTrack.lifeTime+1),curTrack.forceMag(curTrack.startingFrameExtra:curTrack.endingFrameExtra));
            curEndFrame = min(curTrack.startingFrameExtra+periodFrames-1,curTrack.endingFrame);
            curEarlyPeriod = curEndFrame - curTrack.startingFrameExtra+1;
            [~,curForceEarlySlopeGroup(rr)] = regression((1:curEarlyPeriod),curTrack.forceMag(curTrack.startingFrameExtra:curEndFrame));
            sF=curTrack.startingFrameExtraExtra; eF=curTrack.endingFrameExtraExtra;
            curEndFrame = min(sF+periodFrames-1,curTrack.endingFrame);
            curEarlyPeriod = curEndFrame - sF+1;
            [~,curVeryEarlyAmpSlopeGroup(rr)] = regression((1:curEarlyPeriod),curTrack.ampTotal(sF:curEndFrame));

            % There is amplitude that has rising phase but just stops in the middle
            % of course. (It means the the particle loses it's point-ness while
            % increasing the intensity. This does not help in characterizing true
            % G1 behavior
            % I'll filter these out with lateAmpSlope
            curStartFrame = max(eF-periodFrames,curTrack.startingFrame);
            curLatePeriod = curTrack.endingFrameExtraExtra - curStartFrame + 1;
            [~,curForceLateSlopeGroup(rr)] = regression((1:curLatePeriod),curTrack.ampTotal(curStartFrame:eF));
        end
        indFlatForce = curForceSlopeGroup<=0 & curForceEarlySlopeGroup<=0;
        validTrackID = (~indFlatForce & ~indFlatAmp) & curVeryEarlyAmpSlopeGroup>0 ...
            & curForceLateSlopeGroup<0;

        tracksG1fromG1_1=curTracksNAtalG1(validTrackID);

        %Potential G2 in G1
        ampSlopeG1 = arrayfun(@(x) x.ampSlope,curTracksNAtalG1);
        lifeTimeG1 = arrayfun(@(x) x.lifeTime,curTracksNAtalG1);
        ampEndingG1 = arrayfun(@(x) x.ampTotal(x.endingFrameExtra),curTracksNAtalG1);
        ampStartingG1 = arrayfun(@(x) x.ampTotal(x.startingFrameExtra),curTracksNAtalG1);
        idxIncreasingAmpG2inG1 = ampSlopeG1>0 & ampEndingG1>ampStartingG1;
        % idxLowInitForceG2= initForceG2<500;
        % maturing NA should have at least 3 min of lifetime
        thresLT_G2_frames_inG1 = 3*60/tInterval;
        idxLongLifeTimeG2inG1=lifeTimeG1>thresLT_G2_frames_inG1;
        idGroup2fromG1 = idxIncreasingAmpG2inG1 & idxLongLifeTimeG2inG1; %& idxLowInitForceG2 
        tracksG2fromG1_1 = curTracksNAtalG1(idGroup2fromG1);

        realTypesFromG1 = zeros(size(tracksG1fromG1_1));
        for kk=1:numel(tracksG1fromG1_1)
            tracksG1fromG1_1(kk)=readIntensityFromTracks(tracksG1fromG1_1(kk),imgMap,1,'extraLength',30,'movieData',curMD);
            tracksG1fromG1_1(kk)=readIntensityFromTracks(tracksG1fromG1_1(kk),tMap,2,'movieData',curMD);
            [firstIncreseTimeIntAgainstForce,forceTransmittingAll...
            ,firstIncreseTimeIntAll,firstIncreseTimeForce,bkgMaxIntAll,bkgMaxForce] =...
                calculateFirstIncreaseTimeTracks(tracksG1fromG1_1(kk),splineParamInit,preDetecFactor,tInterval);

            tracksG1fromG1_1(kk).forceTransmitting=forceTransmittingAll;
            tracksG1fromG1_1(kk).firstIncreseTimeInt=firstIncreseTimeIntAll;
            tracksG1fromG1_1(kk).firstIncreseTimeForce=firstIncreseTimeForce;
            tracksG1fromG1_1(kk).firstIncreseTimeIntAgainstForce=firstIncreseTimeIntAgainstForce;
            tracksG1fromG1_1(kk).bkgMaxInt=bkgMaxIntAll;
            tracksG1fromG1_1(kk).bkgMaxForce=bkgMaxForce;

            hSum=showSingleAdhesionTrackSummary(curMD,tracksG1fromG1_1(kk),imgMap,tMap,kk);
            hBcc=figure; plot(crossVariance(tracksG1fromG1_1(kk).ampTotal,tracksG1fromG1_1(kk).forceMag,13))
            [ cur_tLagBccPosG2,cur_tLagBccNegG2,cur_BccInLaterG2,cur_BccInPriorG2 ] ...
                = calculateBccPeakLag( tracksG1fromG1_1(kk).ampTotal, tracksG1fromG1_1(kk).forceMag , tracksG1fromG1_1(kk).startingFrameExtra,curMD.timeInterval_);        
            disp(['cur_tLagBccPosG2:' num2str(cur_tLagBccPosG2) ' cur_BccInLaterG2:' num2str(cur_BccInLaterG2) ' cur_tLagBccNegG2:' num2str(cur_tLagBccNegG2) ' cur_BccInPriorG2:' num2str(cur_BccInPriorG2) ])

            s=(inputdlg('What class is this NA? (1-9): '));
            try
                realTypesFromG1(kk)=str2num(s{1});
            catch
                realTypesFromG1(kk)=3;
            end
            close(hSum); close(hBcc);
            disp(['G1: ' num2str(kk) ' out of ' num2str(numel(tracksG1fromG1_1)) ' done.'])
        end
        tracksG2fromG1_2 = tracksG1fromG1_1(realTypesFromG1==2);
        tracksG1fromG1_2 = tracksG1fromG1_1(realTypesFromG1==1);


        % Now curTracksNAtalG2

        ampSlopeG2 = arrayfun(@(x) x.ampSlope,curTracksNAtalG2);
        initForceG2 = arrayfun(@(x) x.forceMag(x.startingFrame),curTracksNAtalG2);
        lifeTimeG2 = arrayfun(@(x) x.lifeTime,curTracksNAtalG2);
        ampEndingG2 = arrayfun(@(x) x.ampTotal(x.endingFrameExtra),curTracksNAtalG2);
        ampStartingG2 = arrayfun(@(x) x.ampTotal(x.startingFrameExtra),curTracksNAtalG2);
        idxIncreasingAmpG2 = ampSlopeG2>0 & ampEndingG2>ampStartingG2;
        % idxLowInitForceG2= initForceG2<500;
        % maturing NA should have at least 3 min of lifetime
        thresLT_G2_frames = 2*60/tInterval;
        idxLongLifeTimeG2=lifeTimeG2>thresLT_G2_frames;
        realG2tal = curForceTransmittingG2 & curStartingAmpG2<thresAmp & stateFA;
        idGroup2valid = idxIncreasingAmpG2 & idxLongLifeTimeG2 & realG2tal; %& idxLowInitForceG2 

        curTracksNAtalG2f = curTracksNAtalG2(idGroup2valid);

        realTypesFromG2 = zeros(numel(curTracksNAtalG2f),1);
        for kk=1:numel(curTracksNAtalG2f)
            curTracksNAtalG2f(kk)=readIntensityFromTracks(curTracksNAtalG2f(kk),imgMap,1,'extraLength',30,'movieData',curMD);
            curTracksNAtalG2f(kk)=readIntensityFromTracks(curTracksNAtalG2f(kk),tMap,2,'movieData',curMD);
            [firstIncreseTimeIntAgainstForce,forceTransmittingAll...
            ,firstIncreseTimeIntAll,firstIncreseTimeForce,bkgMaxIntAll,bkgMaxForce] =...
                calculateFirstIncreaseTimeTracks(curTracksNAtalG2f(kk),splineParamInit,preDetecFactor,tInterval);
            curTracksNAtalG2f(kk).forceTransmitting=forceTransmittingAll;
            curTracksNAtalG2f(kk).firstIncreseTimeInt=firstIncreseTimeIntAll;
            curTracksNAtalG2f(kk).firstIncreseTimeForce=firstIncreseTimeForce;
            curTracksNAtalG2f(kk).firstIncreseTimeIntAgainstForce=firstIncreseTimeIntAgainstForce;
            curTracksNAtalG2f(kk).bkgMaxInt=bkgMaxIntAll;
            curTracksNAtalG2f(kk).bkgMaxForce=bkgMaxForce;

            hSum=showSingleAdhesionTrackSummary(curMD,curTracksNAtalG2f(kk),imgMap,tMap,kk);
            hBcc=figure; plot(crossVariance(curTracksNAtalG2f(kk).ampTotal,curTracksNAtalG2f(kk).forceMag,13))
            [ cur_tLagBccPosG2,cur_tLagBccNegG2,cur_BccInLaterG2,cur_BccInPriorG2 ] ...
                = calculateBccPeakLag( curTracksNAtalG2f(kk).ampTotal, curTracksNAtalG2f(kk).forceMag , curTracksNAtalG2f(kk).startingFrameExtra,curMD.timeInterval_);        
            disp(['cur_tLagBccPosG2:' num2str(cur_tLagBccPosG2) ' cur_tLagBccNegG2:' num2str(cur_tLagBccNegG2) ' cur_BccInLaterG2:' num2str(cur_BccInLaterG2) ' cur_BccInPriorG2:' num2str(cur_BccInPriorG2) ])
            s=(inputdlg('What class is this NA? (1-9): '));
            try
                realTypesFromG2(kk)=str2num(s{1});
            catch
                realTypesFromG2(kk)=3;
            end
            close(hSum); close(hBcc);
            disp(['G2: ' num2str(kk) ' out of ' num2str(numel(curTracksNAtalG2f)) ' done.'])
        end
        tracksG2fromG2 = curTracksNAtalG2f(realTypesFromG2==2);
    %     tracksG2fromG2 = curTracksNAtalG1(realTypesFromG2==2);
        tracksG1fromG2 = curTracksNAtalG2f(realTypesFromG2==1);
    %     
        curTracksNAtalG2f = [tracksG2fromG1_1; tracksG2fromG1_2; tracksG2fromG2];
        curTracksNAtalG1 = [tracksG1fromG1_2; tracksG1fromG2];
        tracksG2=curTracksNAtalG2f; tracksG1=curTracksNAtalG1;
        save([talFolder{ii} filesep 'data' filesep 'tracksG2real.mat'],'tracksG2');
        save([talFolder{ii} filesep 'data' filesep 'tracksG1real.mat'],'tracksG1');
    end  
    % Get the time measurement
    tInterval = curMD.timeInterval_;
    [curTimeDelayG1] = calculateFirstIncreaseTimeTracks(curTracksNAtalG1,splineParamInit,preDetecFactor,tInterval);
    disp(nanmedian(curTimeDelayG1))
    [curTimeDelayG2] = calculateFirstIncreaseTimeTracks(curTracksNAtalG2f,splineParamInit,preDetecFactor,tInterval);
    disp(nanmedian(curTimeDelayG2))
%         curTimeDelay = arrayfun(@(x) x.firstIncreseTimeIntAgainstForce, curTracksNAtalG2(realG2tal));
%         curTimeDelay = arrayfun(@(x) x.firstIncreseTimeIntAgainstForce, curTracksNAtalG2(curForceTransmittingG2 & curStartingAmpG2<thresAmp));
    InitTimeLagMatureTal{ii} = curTimeDelayG2;
%     tracksNAtalG2{ii}=curTracksNAtalG2(curForceTransmittingG2 & curStartingAmpG2<thresAmp);
    disp([talFolder{ii} ' has been analyzed for t_initForMaturingAdhesion.'])


    tlagInitTal{ii} = curTimeDelayG1; %cur_tlagInitTal;
    totalInitTLagTal = [totalInitTLagTal; curTimeDelayG1]; %cur_tlagInitTal];

    cur_tlagPeakTal = load([talFolder{ii} '/data/timePeaksG1.mat'],'peakTimeIntAgainstForceAll');
    cur_tlagPeakTal = cur_tlagPeakTal.peakTimeIntAgainstForceAll;
    if ii>1
        cur_tlagPeakTal=cur_tlagPeakTal*2;
    end
    tlagPeakTal{ii} = cur_tlagPeakTal;
    totalPeakTLagTal = [totalPeakTLagTal; cur_tlagPeakTal];
    
    % Get the assembly rate and force growth rate
    [~,curEarlyAssmRateG1] = arrayfun(@(x) regression(1:3,x.ampTotal(x.startingFrameExtra:x.startingFrameExtra+2)),curTracksNAtalG1);
    [~,curEarlyAssmRateG2] = arrayfun(@(x) regression(1:5,x.ampTotal(x.startingFrameExtra:x.startingFrameExtra+4)),curTracksNAtalG2f);
    earlyAssmRateG1All{ii} = curEarlyAssmRateG1; %cur_tlagInitPax;
    earlyAssmRateG2All{ii} = curEarlyAssmRateG2; %cur_tlagInitPax;
    
    [~,earlyForceRateG1] = arrayfun(@(x) regression(1:3,x.forceMag(x.startingFrameExtra:x.startingFrameExtra+2)),curTracksNAtalG1);
    [~,earlyForceRateG2] = arrayfun(@(x) regression(1:3,x.forceMag(x.startingFrameExtra:x.startingFrameExtra+2)),curTracksNAtalG2f);
    earlyForceRateG1All{ii} = earlyForceRateG1; %cur_tlagInitPax;
    earlyForceRateG2All{ii} = earlyForceRateG2; %cur_tlagInitPax;
    
    % Bcc measurement - G1
    tLagBccNeg = NaN(numel(curTracksNAtalG1),1);
    BccInPrior = NaN(numel(curTracksNAtalG1),1);
    tLagBccPos = NaN(numel(curTracksNAtalG1),1);
    BccInLater = NaN(numel(curTracksNAtalG1),1);
    for pp=1:numel(curTracksNAtalG1)
        curTrack = curTracksNAtalG1(pp);
        if length(curTrack.forceMag)<length(curTrack.ampTotal)
            tempForce=curTrack.forceMag;
            curTrack.forceMag=curTrack.ampTotal;
            curTrack.forceMag(1:length(tempForce))=tempForce;
        end
        curBcc = crossVariance(curTrack.ampTotal,curTrack.forceMag,11);
        %Will find the peak against curTrack.startingFrameExtra to the
        %front and back
        maxBcc=nanmax(curBcc);
        indCands=locmax1d(curBcc,5);
        % Exclude the locmaxes at the ends
        firstInd = find(~isnan(curBcc),1);
        lastInd = find(~isnan(curBcc),1,'last');
        indCands(ismember(indCands,[firstInd,lastInd]))=[];
        % if the loc max is close to curTrack.startingFrameExtra and the
        % actual value is above 0.8*maxBcc, we'll admit the time as
        % BccShift time
        
        % Sort based on proximity to curTrack.startingFrameExtra
        SFE = curTrack.startingFrameExtra;
        frameDiff = indCands-SFE;
        
        if any(frameDiff<0)
            indNegDiff= indCands(frameDiff<0);
            for curIndCand = indNegDiff'
                if curBcc(curIndCand)>0.5*maxBcc
                    tLagBccNeg(pp) = frameDiff(indCands==curIndCand)*curMD.timeInterval_;
                    BccInPrior(pp) = curBcc(curIndCand);
                    break
                end
            end
        end
        indPosDiff= indCands(frameDiff>=0);

        for curIndCand = indPosDiff'
            if curBcc(curIndCand)>0.7*maxBcc
                tLagBccPos(pp) = frameDiff(indCands==curIndCand)*curMD.timeInterval_;
                BccInLater(pp) = curBcc(curIndCand);
                break
            end
        end
    end
    tLagBccPosG1All{ii} = tLagBccPos; 
    tLagBccNegG1All{ii} = tLagBccNeg; 
    BccInLaterG1All{ii} = BccInLater; 
    BccInPriorG1All{ii} = BccInPrior; 

    % Bcc measurement - G2
    tLagBccNegG2 = NaN(numel(curTracksNAtalG2f),1);
    BccInPriorG2 = NaN(numel(curTracksNAtalG2f),1);
    tLagBccPosG2 = NaN(numel(curTracksNAtalG2f),1);
    BccInLaterG2 = NaN(numel(curTracksNAtalG2f),1);
    for pp=1:numel(curTracksNAtalG2f)
        curTrack = curTracksNAtalG2f(pp);
        if length(curTrack.forceMag)<length(curTrack.ampTotal)
            tempForce=curTrack.forceMag;
            curTrack.forceMag=curTrack.ampTotal;
            curTrack.forceMag(1:length(tempForce))=tempForce;
        end
        curBcc = crossVariance(curTrack.ampTotal,curTrack.forceMag,11);
        %Will find the peak against curTrack.startingFrameExtra to the
        %front and back
        maxBcc=nanmax(curBcc);
        indCands=locmax1d(curBcc,5);
        % Exclude the locmaxes at the ends
        firstInd = find(~isnan(curBcc),1);
        lastInd = find(~isnan(curBcc),1,'last');
        indCands(ismember(indCands,[firstInd,lastInd]))=[];
        % if the loc max is close to curTrack.startingFrameExtra and the
        % actual value is above 0.8*maxBcc, we'll admit the time as
        % BccShift time
        
        % Sort based on proximity to curTrack.startingFrameExtra
        SFE = curTrack.startingFrameExtra;
        frameDiff = indCands-SFE;
        
        if any(frameDiff<0)
            indNegDiff= indCands(frameDiff<0);
            for curIndCand = indNegDiff'
                if curBcc(curIndCand)>0.4*maxBcc
                    tLagBccNegG2(pp) = frameDiff(indCands==curIndCand)*curMD.timeInterval_;
                    BccInPriorG2(pp) = curBcc(curIndCand);
                    break
                end
            end
        end
        indPosDiff= indCands(frameDiff>=0);

        for curIndCand = indPosDiff'
            if curBcc(curIndCand)>0.5*maxBcc
                tLagBccPosG2(pp) = frameDiff(indCands==curIndCand)*curMD.timeInterval_;
                BccInLaterG2(pp) = curBcc(curIndCand);
                break
            end
        end
    end
    tLagBccPosG2All{ii} = tLagBccPosG2; 
    tLagBccNegG2All{ii} = tLagBccNegG2; 
    BccInLaterG2All{ii} = BccInLaterG2; 
    BccInPriorG2All{ii} = BccInPriorG2; 
end
totalInitTLagTal =cell2mat(tlagInitTal(1:numTalMovies)');
totalInitTimeLagMatureTal =cell2mat(InitTimeLagMatureTal(1:numTalMovies)');
% totalInitTLagTal=cell2mat(tlagInitTal);%[tlagInitTal{1}];
totalPeakTLagTal=cell2mat(tlagPeakTal(1:numTalMovies)');

tLagBccPosG1Talin =cell2mat(tLagBccPosG1All(1:numTalMovies));
tLagBccNegG1Talin =cell2mat(tLagBccNegG1All(1:numTalMovies));
tLagBccPosG2Talin =cell2mat(tLagBccPosG2All(1:numTalMovies));
tLagBccNegG2Talin =cell2mat(tLagBccNegG2All(1:numTalMovies));
BccInPriorG1Talin =cell2mat(BccInPriorG1All(1:numTalMovies));
BccInLaterG1Talin =cell2mat(BccInLaterG1All(1:numTalMovies));
BccInPriorG2Talin =cell2mat(BccInPriorG2All(1:numTalMovies));
BccInLaterG2Talin =cell2mat(BccInLaterG2All(1:numTalMovies));

totalEarlyAssmRateG1Tal= cell2mat(earlyAssmRateG1All(1:numTalMovies)');
totalEarlyAssmRateG2Tal= cell2mat(earlyAssmRateG2All(1:numTalMovies)');
totalEarlyForceRateG1Tal= cell2mat(earlyForceRateG1All(1:numTalMovies)');
totalEarlyForceRateG2Tal= cell2mat(earlyForceRateG2All(1:numTalMovies)');

% Get all G1 from talin rep exp
% Calculate Bcc and Peak Bcc time in reference to startingFrameExtra
%% clearing large variables
clear talImgStack talForceStack tMap paxForceStack paxImgStack imgMap forceFiledPax forceFieldVin forceFieldTal displFieldPax displFieldTal displFieldVin curTracksNA*
%% for vinculin 
tLagBccPosG1All=cell(numVinMovies,1);
tLagBccNegG1All=cell(numVinMovies,1);
BccInLaterG1All=cell(numVinMovies,1);
BccInPriorG1All=cell(numVinMovies,1);
tLagBccPosG2All=cell(numVinMovies,1);
tLagBccNegG2All=cell(numVinMovies,1);
BccInLaterG2All=cell(numVinMovies,1);
BccInPriorG2All=cell(numVinMovies,1);
preDetecFactor=1/10; %for vinculin 
%% running
for ii=1:numVinMovies
%     cur_tlagInitVin = load([vinFolder{ii} '/data/timeInitLagsG1.mat'],'firstIncreseTimeIntAgainstForceAll');
%     cur_tlagInitVin = cur_tlagInitVin.firstIncreseTimeIntAgainstForceAll;
%     tlagInitVin{ii} = cur_tlagInitVin;
    % try to filter G2 based on starting amplitude, peaks and slopes...
    if exist([vinFolder{ii} filesep 'data' filesep 'tracksG2real.mat'],'file')
        tracksG2 = load([vinFolder{ii} filesep 'data' filesep 'tracksG2real.mat'],'tracksG2');
        curTracksNAvinG2=tracksG2.tracksG2;
        tracksG1 = load([vinFolder{ii} filesep 'data' filesep 'tracksG1real.mat'],'tracksG1');
        curTracksNAvinG1=tracksG1.tracksG1;
        curMDpath = fileparts(fileparts(vinFolder{ii}));
        curMD = load([curMDpath filesep 'movieData.mat']);
        curMD = curMD.MD;
    else
        curTracksNAvinG2 = load([vinFolder{ii} filesep 'data' filesep 'tracksG2.mat'],'tracksG2');
        curTracksNAvinG2 = curTracksNAvinG2.tracksG2; 
        % Get rid of overlapping tracks
        curTracksNAvinG2 = filterOverlappingTracks(curTracksNAvinG2);
        curForceTransmittingG2 = arrayfun(@(x) ~isempty(x.firstIncreseTimeIntAgainstForce), curTracksNAvinG2);
        curStartingAmpG2 = arrayfun(@(x) x.ampTotal(x.startingFrameExtra), curTracksNAvinG2);

        curTracksNAvinG1 = load([vinFolder{ii} filesep 'data' filesep 'tracksG1.mat'],'tracksG1');
        curTracksNAvinG1 = curTracksNAvinG1.tracksG1; 
        curTimeDelay = arrayfun(@(x) x.firstIncreseTimeIntAgainstForce, curTracksNAvinG1,'UniformOutput',false);
        curStartingAmpG1 = arrayfun(@(x) x.ampTotal(x.startingFrameExtra), curTracksNAvinG1);
        thresAmp = mean(curStartingAmpG1)+3*std(curStartingAmpG1);

        % filtering for G2
        stateFA = arrayfun(@(x) any(cellfun(@(y) strcmp(y,'FC') || strcmp(y,'FA'), x.state)), curTracksNAvinG2);
        realG2vin = curForceTransmittingG2 & curStartingAmpG2<thresAmp & stateFA;
        curTimeDelay = arrayfun(@(x) x.firstIncreseTimeIntAgainstForce, curTracksNAvinG2(realG2vin));

        % Go through each adhesion track and make sure
        curMDpath = fileparts(fileparts(vinFolder{ii}));
        curMD = load([curMDpath filesep 'movieData.mat']);
        curMD = curMD.MD;
        tMap = load([vinFolder{ii} filesep 'fMap' filesep 'tMap.mat'],'tMap');
        tMap = tMap.tMap;
        imgMap = load([vinFolder{ii} filesep 'pax' filesep 'paxImgStack.mat'],'paxImgStack');
        imgMap = imgMap.paxImgStack;
        tInterval = curMD.timeInterval_;

        % Get rid of overlapping tracks
        curTracksNAvinG1 = filterOverlappingTracks(curTracksNAvinG1);
        % Get rid of tracks that are obviously inside the cell
        distG1 = arrayfun(@(x) mean(x.distToEdge),curTracksNAvinG1);
        thresDist = 10; % um
        curTracksNAvinG1 = curTracksNAvinG1(distG1<thresDist);

        [firstIncreseTimeIntAgainstForce,forceTransmittingAll...
        ,firstIncreseTimeIntAll,firstIncreseTimeForce,bkgMaxIntAll,bkgMaxForce] =...
            calculateFirstIncreaseTimeTracks(curTracksNAvinG1,splineParamInit,preDetecFactor,tInterval);
        curTracksNAvinG1 = curTracksNAvinG1(forceTransmittingAll);
        firstIncreseTimeIntAll=firstIncreseTimeIntAll(forceTransmittingAll);
        firstIncreseTimeForce=firstIncreseTimeForce(forceTransmittingAll);
        firstIncreseTimeIntAgainstForce=firstIncreseTimeIntAgainstForce(forceTransmittingAll);
        bkgMaxIntAll=bkgMaxIntAll(forceTransmittingAll);
        bkgMaxForce=bkgMaxForce(forceTransmittingAll);
        forceTransmittingAll=forceTransmittingAll(forceTransmittingAll);

        % Further G1 critera
        % 1. Based on amplitude 
        curAmpSlopeGroup = arrayfun(@(x) x.ampSlope, curTracksNAvinG1);
        curEarlyAmpSlopeGroup = arrayfun(@(x) x.earlyAmpSlope, curTracksNAvinG1);
        % We decided to regard amplitude with flat slope as noise.
        %             figure, plot(curAmpSlopeGroup,curEarlyAmpSlopeGroup,'*')
        indFlatAmp = curAmpSlopeGroup<=0 & curEarlyAmpSlopeGroup<=0;
        % 2. Based on forceMag
        curVeryEarlyAmpSlopeGroup = NaN(numel(curTracksNAvinG1),1);
        curForceSlopeGroup = NaN(numel(curTracksNAvinG1),1);
        curForceEarlySlopeGroup = NaN(numel(curTracksNAvinG1),1);
        curAmpOverallSlopeGroup = NaN(numel(curTracksNAvinG1),1);
        curAmpOverallSlopeUpErrGroup = NaN(numel(curTracksNAvinG1),1);
        curAmpOverallSlopeDownErrGroup = NaN(numel(curTracksNAvinG1),1);
        curForceLateSlopeGroup = NaN(numel(curTracksNAvinG1),1);

        periodFrames = 30;
        rr=0;
        for pp=1:numel(curTracksNAvinG1)
            rr=rr+1;
            curTrack = curTracksNAvinG1(pp);
            [~,curForceSlopeGroup(rr)] = regression((1:curTrack.lifeTime+1),curTrack.forceMag(curTrack.startingFrameExtra:curTrack.endingFrameExtra));
            curEndFrame = min(curTrack.startingFrameExtra+periodFrames-1,curTrack.endingFrame);
            curEarlyPeriod = curEndFrame - curTrack.startingFrameExtra+1;
            [~,curForceEarlySlopeGroup(rr)] = regression((1:curEarlyPeriod),curTrack.forceMag(curTrack.startingFrameExtra:curEndFrame));
            sF=curTrack.startingFrameExtraExtra; eF=curTrack.endingFrameExtraExtra;
            curEndFrame = min(sF+periodFrames-1,curTrack.endingFrame);
            curEarlyPeriod = curEndFrame - sF+1;
            [~,curVeryEarlyAmpSlopeGroup(rr)] = regression((1:curEarlyPeriod),curTrack.ampTotal(sF:curEndFrame));

            % There is amplitude that has rising phase but just stops in the middle
            % of course. (It means the the particle loses it's point-ness while
            % increasing the intensity. This does not help in characterizing true
            % G1 behavior
            % I'll filter these out with lateAmpSlope
            curStartFrame = max(eF-periodFrames,curTrack.startingFrame);
            curLatePeriod = curTrack.endingFrameExtraExtra - curStartFrame + 1;
            [~,curForceLateSlopeGroup(rr)] = regression((1:curLatePeriod),curTrack.ampTotal(curStartFrame:eF));
        end
        indFlatForce = curForceSlopeGroup<=0 & curForceEarlySlopeGroup<=0;
        validTrackID = (~indFlatForce & ~indFlatAmp) & curVeryEarlyAmpSlopeGroup>0 ...
            & curForceLateSlopeGroup<0;

        tracksG1fromG1_1=curTracksNAvinG1(validTrackID);

        %Potential G2 in G1
        ampSlopeG1 = arrayfun(@(x) x.ampSlope,curTracksNAvinG1);
        lifeTimeG1 = arrayfun(@(x) x.lifeTime,curTracksNAvinG1);
        ampEndingG1 = arrayfun(@(x) x.ampTotal(x.endingFrameExtra),curTracksNAvinG1);
        ampStartingG1 = arrayfun(@(x) x.ampTotal(x.startingFrameExtra),curTracksNAvinG1);
        idxIncreasingAmpG2inG1 = ampSlopeG1>0 & ampEndingG1>ampStartingG1;
        % idxLowInitForceG2= initForceG2<500;
        % maturing NA should have at least 3 min of lifetime
        thresLT_G2_frames_inG1 = 3*60/tInterval;
        idxLongLifeTimeG2inG1=lifeTimeG1>thresLT_G2_frames_inG1;
        idGroup2fromG1 = idxIncreasingAmpG2inG1 & idxLongLifeTimeG2inG1; %& idxLowInitForceG2 
        tracksG2fromG1_1 = curTracksNAvinG1(idGroup2fromG1);

        realTypesFromG1 = zeros(size(tracksG1fromG1_1));
        for kk=1:numel(tracksG1fromG1_1)
            tracksG1fromG1_1(kk)=readIntensityFromTracks(tracksG1fromG1_1(kk),imgMap,1,'extraLength',30,'movieData',curMD);
            tracksG1fromG1_1(kk)=readIntensityFromTracks(tracksG1fromG1_1(kk),tMap,2,'movieData',curMD);
            [firstIncreseTimeIntAgainstForce,forceTransmittingAll...
            ,firstIncreseTimeIntAll,firstIncreseTimeForce,bkgMaxIntAll,bkgMaxForce] =...
                calculateFirstIncreaseTimeTracks(tracksG1fromG1_1(kk),splineParamInit,preDetecFactor,tInterval);

            tracksG1fromG1_1(kk).forceTransmitting=forceTransmittingAll;
            tracksG1fromG1_1(kk).firstIncreseTimeInt=firstIncreseTimeIntAll;
            tracksG1fromG1_1(kk).firstIncreseTimeForce=firstIncreseTimeForce;
            tracksG1fromG1_1(kk).firstIncreseTimeIntAgainstForce=firstIncreseTimeIntAgainstForce;
            tracksG1fromG1_1(kk).bkgMaxInt=bkgMaxIntAll;
            tracksG1fromG1_1(kk).bkgMaxForce=bkgMaxForce;

            hSum=showSingleAdhesionTrackSummary(curMD,tracksG1fromG1_1(kk),imgMap,tMap,kk);
            hBcc=figure; plot(crossVariance(tracksG1fromG1_1(kk).ampTotal,tracksG1fromG1_1(kk).forceMag,13))
            [ cur_tLagBccPosG1,cur_tLagBccNegG1,cur_BccInLaterG1,cur_BccInPriorG1 ] ...
                = calculateBccPeakLag( tracksG1fromG1_1(kk).ampTotal, tracksG1fromG1_1(kk).forceMag , tracksG1fromG1_1(kk).startingFrameExtra,curMD.timeInterval_);        
            disp(['tLagPos = ' num2str(cur_tLagBccPosG1) ', BccLater = ' num2str(cur_BccInLaterG1) ' tLagNeg=' num2str(cur_tLagBccNegG1) ' BccPrior:' num2str(cur_BccInPriorG1) ])

            s=(inputdlg('What class is this NA? (1-9): '));
            try
                realTypesFromG1(kk)=str2num(s{1});
            catch
                realTypesFromG1(kk)=3;
            end
            close(hSum); close(hBcc);
            disp(['G1: ' num2str(kk) ' out of ' num2str(numel(tracksG1fromG1_1)) ' done.'])
        end
        tracksG2fromG1_2 = tracksG1fromG1_1(realTypesFromG1==2);
        tracksG1fromG1_2 = tracksG1fromG1_1(realTypesFromG1==1);

        % Now curTracksNAvinG2

        ampSlopeG2 = arrayfun(@(x) x.ampSlope,curTracksNAvinG2);
        initForceG2 = arrayfun(@(x) x.forceMag(x.startingFrame),curTracksNAvinG2);
        lifeTimeG2 = arrayfun(@(x) x.lifeTime,curTracksNAvinG2);
        ampEndingG2 = arrayfun(@(x) x.ampTotal(x.endingFrameExtra),curTracksNAvinG2);
        ampStartingG2 = arrayfun(@(x) x.ampTotal(x.startingFrameExtra),curTracksNAvinG2);
        idxIncreasingAmpG2 = ampSlopeG2>0 & ampEndingG2>ampStartingG2;
        % idxLowInitForceG2= initForceG2<500;
        % maturing NA should have at least 3 min of lifetime
        thresLT_G2_frames = 2*60/tInterval;
        idxLongLifeTimeG2=lifeTimeG2>thresLT_G2_frames;
        realG2 = curForceTransmittingG2 & curStartingAmpG2<thresAmp & stateFA;
        idGroup2valid = idxIncreasingAmpG2 & idxLongLifeTimeG2 & realG2; %& idxLowInitForceG2 

        curTracksNAvinG2f = curTracksNAvinG2(idGroup2valid);

        realTypesFromG2 = zeros(numel(curTracksNAvinG2f),1);
    %     realTypesFromG2 = zeros(size(realG2vin));
        for kk=1:numel(curTracksNAvinG2f)
            curTracksNAvinG2f(kk)=readIntensityFromTracks(curTracksNAvinG2f(kk),imgMap,1,'extraLength',30,'movieData',curMD);
            curTracksNAvinG2f(kk)=readIntensityFromTracks(curTracksNAvinG2f(kk),tMap,2,'movieData',curMD);
            [firstIncreseTimeIntAgainstForce,forceTransmittingAll...
            ,firstIncreseTimeIntAll,firstIncreseTimeForce,bkgMaxIntAll,bkgMaxForce] =...
                calculateFirstIncreaseTimeTracks(curTracksNAvinG2f(kk),splineParamInit,preDetecFactor,tInterval);
            curTracksNAvinG2f(kk).forceTransmitting=forceTransmittingAll;
            curTracksNAvinG2f(kk).firstIncreseTimeInt=firstIncreseTimeIntAll;
            curTracksNAvinG2f(kk).firstIncreseTimeForce=firstIncreseTimeForce;
            curTracksNAvinG2f(kk).firstIncreseTimeIntAgainstForce=firstIncreseTimeIntAgainstForce;
            curTracksNAvinG2f(kk).bkgMaxInt=bkgMaxIntAll;
            curTracksNAvinG2f(kk).bkgMaxForce=bkgMaxForce;

            hSum=showSingleAdhesionTrackSummary(curMD,curTracksNAvinG2f(kk),imgMap,tMap,kk);
            hBcc=figure; plot(crossVariance(curTracksNAvinG2f(kk).ampTotal,curTracksNAvinG2f(kk).forceMag,13))
            [ cur_tLagBccPosG2,cur_tLagBccNegG2,cur_BccInLaterG2,cur_BccInPriorG2 ] ...
                = calculateBccPeakLag( curTracksNAvinG2f(kk).ampTotal, curTracksNAvinG2f(kk).forceMag , curTracksNAvinG2f(kk).startingFrameExtra,curMD.timeInterval_);        
            disp(['tLagPosG2:' num2str(cur_tLagBccPosG2) ', tLagBccNeg: ' num2str(cur_tLagBccNegG2) ', cur_BccInLaterG2: ' num2str(cur_BccInLaterG2) ' cur_BccInPriorG2:' num2str(cur_BccInPriorG2) ])
            s=(inputdlg('What class is this NA? (1-9): '));
            try
                realTypesFromG2(kk)=str2num(s{1});
            catch
                realTypesFromG2(kk)=3;
            end
            close(hSum)
        end
        tracksG2fromG2 = curTracksNAvinG2f(realTypesFromG2==2);
        tracksG1fromG2 = curTracksNAvinG2f(realTypesFromG2==1);

        curTracksNAvinG2 = [tracksG2fromG1_1; tracksG2fromG1_2; tracksG2fromG2]; %[tracksG2fromG1; tracksG2fromG2];
        curTracksNAvinG1 = [tracksG1fromG1_2; tracksG1fromG2];
        tracksG2=curTracksNAvinG2; tracksG1=curTracksNAvinG1;
        save([vinFolder{ii} filesep 'data' filesep 'tracksG2real.mat'],'tracksG2');
        save([vinFolder{ii} filesep 'data' filesep 'tracksG1real.mat'],'tracksG1');
    end
    % Get the time measurement
    tInterval = curMD.timeInterval_;
    [curTimeDelayG2] = calculateFirstIncreaseTimeTracks(curTracksNAvinG2,splineParamInit,preDetecFactor,tInterval);
    disp(nanmedian(curTimeDelayG2))
    [curTimeDelayG1] = calculateFirstIncreaseTimeTracks(curTracksNAvinG1,splineParamInit,preDetecFactor,tInterval);
    disp(nanmedian(curTimeDelayG1))
%         curTimeDelay = arrayfun(@(x) x.firstIncreseTimeIntAgainstForce, curTracksNAvinG2(curForceTransmittingG2 & curStartingAmpG2<thresAmp));
    InitTimeLagMatureVin{ii} = curTimeDelayG2;
    disp([vinFolder{ii} ' has been analyzed for t_initForMaturingAdhesion.'])
    tlagInitVin{ii} = curTimeDelayG1;
    
    cur_tlagPeakVin = load([vinFolder{ii} '/data/timePeaksG1.mat'],'peakTimeIntAgainstForceAll');
    cur_tlagPeakVin = cur_tlagPeakVin.peakTimeIntAgainstForceAll;
    tlagPeakVin{ii} = cur_tlagPeakVin;
    
    % Get the assembly rate and force growth rate
    [~,curEarlyAssmRateG1] = arrayfun(@(x) regression(1:3,x.ampTotal(x.startingFrameExtra:x.startingFrameExtra+2)),curTracksNAvinG1);
    [~,curEarlyAssmRateG2] = arrayfun(@(x) regression(1:5,x.ampTotal(x.startingFrameExtra:x.startingFrameExtra+4)),curTracksNAvinG2);
    earlyAssmRateG1All{ii} = curEarlyAssmRateG1; %cur_tlagInitPax;
    earlyAssmRateG2All{ii} = curEarlyAssmRateG2; %cur_tlagInitPax;
    
    [~,earlyForceRateG1] = arrayfun(@(x) regression(1:3,x.forceMag(x.startingFrameExtra:x.startingFrameExtra+2)),curTracksNAvinG1);
    [~,earlyForceRateG2] = arrayfun(@(x) regression(1:3,x.forceMag(x.startingFrameExtra:x.startingFrameExtra+2)),curTracksNAvinG2);
    earlyForceRateG1All{ii} = earlyForceRateG1; %cur_tlagInitPax;
    earlyForceRateG2All{ii} = earlyForceRateG2; %cur_tlagInitPax;
    
    % Bcc measurement - G1
    tLagBccNeg = NaN(numel(curTracksNAvinG1),1);
    BccInPrior = NaN(numel(curTracksNAvinG1),1);
    tLagBccPos = NaN(numel(curTracksNAvinG1),1);
    BccInLater = NaN(numel(curTracksNAvinG1),1);
    for pp=1:numel(curTracksNAvinG1)
        curTrack = curTracksNAvinG1(pp);
        if length(curTrack.forceMag)<length(curTrack.ampTotal)
            tempForce=curTrack.forceMag;
            curTrack.forceMag=curTrack.ampTotal;
            curTrack.forceMag(1:length(tempForce))=tempForce;
        end
        curBcc = crossVariance(curTrack.ampTotal,curTrack.forceMag,11);
        %Will find the peak against curTrack.startingFrameExtra to the
        %front and back
        maxBcc=nanmax(curBcc);
        indCands=locmax1d(curBcc,5);
        % Exclude the locmaxes at the ends
        firstInd = find(~isnan(curBcc),1);
        lastInd = find(~isnan(curBcc),1,'last');
        indCands(ismember(indCands,[firstInd,lastInd]))=[];
        % if the loc max is close to curTrack.startingFrameExtra and the
        % actual value is above 0.8*maxBcc, we'll admit the time as
        % BccShift time
        
        % Sort based on proximity to curTrack.startingFrameExtra
        SFE = curTrack.startingFrameExtra;
        frameDiff = indCands-SFE;
        
        if any(frameDiff<0)
            indNegDiff= indCands(frameDiff<0);
            for curIndCand = indNegDiff'
                if curBcc(curIndCand)>0.3*maxBcc
                    tLagBccNeg(pp) = frameDiff(indCands==curIndCand)*curMD.timeInterval_;
                    BccInPrior(pp) = curBcc(curIndCand);
                    break
                end
            end
        end
        indPosDiff= indCands(frameDiff>=0);

        for curIndCand = indPosDiff'
            if curBcc(curIndCand)>0.3*maxBcc
                tLagBccPos(pp) = frameDiff(indCands==curIndCand)*curMD.timeInterval_;
                BccInLater(pp) = curBcc(curIndCand);
                break
            end
        end
    end
    tLagBccPosG1All{ii} = tLagBccPos; 
    tLagBccNegG1All{ii} = tLagBccNeg; 
    BccInLaterG1All{ii} = BccInLater; 
    BccInPriorG1All{ii} = BccInPrior; 

    % Bcc measurement - G2
    tLagBccNegG2 = NaN(numel(curTracksNAvinG2),1);
    BccInPriorG2 = NaN(numel(curTracksNAvinG2),1);
    tLagBccPosG2 = NaN(numel(curTracksNAvinG2),1);
    BccInLaterG2 = NaN(numel(curTracksNAvinG2),1);
    for pp=1:numel(curTracksNAvinG2)
        curTrack = curTracksNAvinG2(pp);
        if length(curTrack.forceMag)<length(curTrack.ampTotal)
            tempForce=curTrack.forceMag;
            curTrack.forceMag=curTrack.ampTotal;
            curTrack.forceMag(1:length(tempForce))=tempForce;
        end
        curBcc = crossVariance(curTrack.ampTotal,curTrack.forceMag,11);
        %Will find the peak against curTrack.startingFrameExtra to the
        %front and back
        maxBcc=nanmax(curBcc);
        indCands=locmax1d(curBcc,5);
        % Exclude the locmaxes at the ends
        firstInd = find(~isnan(curBcc),1);
        lastInd = find(~isnan(curBcc),1,'last');
        indCands(ismember(indCands,[firstInd,lastInd]))=[];
        % if the loc max is close to curTrack.startingFrameExtra and the
        % actual value is above 0.8*maxBcc, we'll admit the time as
        % BccShift time
        
        % Sort based on proximity to curTrack.startingFrameExtra
        SFE = curTrack.startingFrameExtra;
        frameDiff = indCands-SFE;
        
        if any(frameDiff<0)
            indNegDiff= indCands(frameDiff<0);
            for curIndCand = indNegDiff'
                if curBcc(curIndCand)>0.3*maxBcc
                    tLagBccNegG2(pp) = frameDiff(indCands==curIndCand)*curMD.timeInterval_;
                    BccInPriorG2(pp) = curBcc(curIndCand);
                    break
                end
            end
        end
        indPosDiff= indCands(frameDiff>=0);

        for curIndCand = indPosDiff'
            if curBcc(curIndCand)>0.3*maxBcc
                tLagBccPosG2(pp) = frameDiff(indCands==curIndCand)*curMD.timeInterval_;
                BccInLaterG2(pp) = curBcc(curIndCand);
                break
            end
        end
    end
    tLagBccPosG2All{ii} = tLagBccPosG2; 
    tLagBccNegG2All{ii} = tLagBccNegG2; 
    BccInLaterG2All{ii} = BccInLaterG2; 
    BccInPriorG2All{ii} = BccInPriorG2; 
    
%     clear curTracksNAvinG2 curTracksNAvinG1 imgMap tMap
end
totalInitTLagVin =cell2mat(tlagInitVin(1:numVinMovies)');
totalInitTimeLagMatureVin =cell2mat(InitTimeLagMatureVin(1:numVinMovies)');
totalPeakTLagVin=cell2mat(tlagPeakVin(1:numVinMovies)');

tLagBccPosG1Vin =cell2mat(tLagBccPosG1All(1:numVinMovies));
tLagBccNegG1Vin =cell2mat(tLagBccNegG1All(1:numVinMovies));
tLagBccPosG2Vin =cell2mat(tLagBccPosG2All(1:numVinMovies));
tLagBccNegG2Vin =cell2mat(tLagBccNegG2All(1:numVinMovies));
BccInPriorG1Vin =cell2mat(BccInPriorG1All(1:numVinMovies));
BccInLaterG1Vin =cell2mat(BccInLaterG1All(1:numVinMovies));
BccInPriorG2Vin =cell2mat(BccInPriorG2All(1:numVinMovies));
BccInLaterG2Vin =cell2mat(BccInLaterG2All(1:numVinMovies));

totalEarlyAssmRateG1Vin=cell2mat(earlyAssmRateG1All(1:numVinMovies)');
totalEarlyAssmRateG2Vin=cell2mat(earlyAssmRateG2All(1:numVinMovies)');
totalEarlyForceRateG1Vin=cell2mat(earlyForceRateG1All(1:numVinMovies)');
totalEarlyForceRateG2Vin=cell2mat(earlyForceRateG2All(1:numVinMovies)');
% clear earlyAssmRateG2All earlyAssmRateG1All earlyForceRateG1All earlyForceRateG2All
%% Initial time lag - paxillin
splineParamInit=0.8;
preDetecFactor=1/10; %for paxillin because diffuse signal before point-source-like signal
tLagBccPosG1All=cell(numPaxMovies,1);
tLagBccNegG1All=cell(numPaxMovies,1);
BccInLaterG1All=cell(numPaxMovies,1);
BccInPriorG1All=cell(numPaxMovies,1);
tLagBccPosG2All=cell(numPaxMovies,1);
tLagBccNegG2All=cell(numPaxMovies,1);
BccInLaterG2All=cell(numPaxMovies,1);
BccInPriorG2All=cell(numPaxMovies,1);
%% running
for ii=1:numPaxMovies
%     cur_tlagInitPax = load([paxFolder{ii} '/data/timeInitLagsG1.mat'],'firstIncreseTimeIntAgainstForceAll');
%     cur_tlagInitPax = cur_tlagInitPax.firstIncreseTimeIntAgainstForceAll;
        % try to filter G2 based on starting amplitude, peaks and slopes...
    if exist([paxFolder{ii} filesep 'data' filesep 'tracksG2real.mat'],'file')
        tracksG2 = load([paxFolder{ii} filesep 'data' filesep 'tracksG2real.mat'],'tracksG2');
        curTracksNApaxG2f=tracksG2.tracksG2;
        tracksG1 = load([paxFolder{ii} filesep 'data' filesep 'tracksG1real.mat'],'tracksG1');
        curTracksNApaxG1=tracksG1.tracksG1;
        curMDpath = fileparts(fileparts(paxFolder{ii}));
        curMD = load([curMDpath filesep 'movieData.mat']);
        curMD = curMD.MD;
    else
        curTracksNApaxG2 = load([paxFolder{ii} filesep 'data' filesep 'tracksG2.mat'],'tracksG2');
        curTracksNApaxG2 = curTracksNApaxG2.tracksG2; 
        % Get rid of overlapping tracks
        curTracksNApaxG2 = filterOverlappingTracks(curTracksNApaxG2);
        curForceTransmittingG2 = arrayfun(@(x) ~isempty(x.firstIncreseTimeIntAgainstForce), curTracksNApaxG2);
        curStartingAmpG2 = arrayfun(@(x) x.ampTotal(x.startingFrameExtra), curTracksNApaxG2);

        curTracksNApaxG1 = load([paxFolder{ii} filesep 'data' filesep 'tracksG1.mat'],'tracksG1');
        curTracksNApaxG1 = curTracksNApaxG1.tracksG1; 
        curTimeDelay = arrayfun(@(x) x.firstIncreseTimeIntAgainstForce, curTracksNApaxG1,'UniformOutput',false);
        curStartingAmpG1 = arrayfun(@(x) x.ampTotal(x.startingFrameExtra), curTracksNApaxG1);
        thresAmp = mean(curStartingAmpG1)+3*std(curStartingAmpG1);

        % filtering for G2
        stateFA = arrayfun(@(x) any(cellfun(@(y) strcmp(y,'FC') || strcmp(y,'FA'), x.state)), curTracksNApaxG2);
        realG2pax = curForceTransmittingG2 & curStartingAmpG2<thresAmp; % & stateFA;
        curTimeDelay = arrayfun(@(x) x.firstIncreseTimeIntAgainstForce, curTracksNApaxG2(realG2pax));

        % Go through each adhesion track and make sure
        curMDpath = fileparts(fileparts(paxFolder{ii}));
        curMD = load([curMDpath filesep 'movieData.mat']);
        curMD = curMD.MD;
        tMap = load([paxFolder{ii} filesep 'fMap' filesep 'tMap.mat'],'tMap');
        tMap = tMap.tMap;
        imgMap = load([paxFolder{ii} filesep 'pax' filesep 'paxImgStack.mat'],'paxImgStack');
        imgMap = imgMap.paxImgStack;
        tInterval = curMD.timeInterval_;
        % Get rid of overlapping tracks
        curTracksNApaxG1 = filterOverlappingTracks(curTracksNApaxG1);
        % Get rid of tracks that are obviously inside the cell
        distG1 = arrayfun(@(x) mean(x.distToEdge),curTracksNApaxG1);
        thresDist = 10; % um
        curTracksNApaxG1 = curTracksNApaxG1(distG1<thresDist);

        [firstIncreseTimeIntAgainstForce,forceTransmittingAll...
        ,firstIncreseTimeIntAll,firstIncreseTimeForce,bkgMaxIntAll,bkgMaxForce] =...
            calculateFirstIncreaseTimeTracks(curTracksNApaxG1,splineParamInit,preDetecFactor,tInterval);
        curTracksNApaxG1 = curTracksNApaxG1(forceTransmittingAll);
        firstIncreseTimeIntAll=firstIncreseTimeIntAll(forceTransmittingAll);
        firstIncreseTimeForce=firstIncreseTimeForce(forceTransmittingAll);
        firstIncreseTimeIntAgainstForce=firstIncreseTimeIntAgainstForce(forceTransmittingAll);
        bkgMaxIntAll=bkgMaxIntAll(forceTransmittingAll);
        bkgMaxForce=bkgMaxForce(forceTransmittingAll);
        forceTransmittingAll=forceTransmittingAll(forceTransmittingAll);

        % Further G1 critera
        % 1. Based on amplitude 
        curAmpSlopeGroup = arrayfun(@(x) x.ampSlope, curTracksNApaxG1);
        curEarlyAmpSlopeGroup = arrayfun(@(x) x.earlyAmpSlope, curTracksNApaxG1);
        % We decided to regard amplitude with flat slope as noise.
        %             figure, plot(curAmpSlopeGroup,curEarlyAmpSlopeGroup,'*')
        indFlatAmp = curAmpSlopeGroup<=0 & curEarlyAmpSlopeGroup<=0;
        % 2. Based on forceMag
        curVeryEarlyAmpSlopeGroup = NaN(numel(curTracksNApaxG1),1);
        curForceSlopeGroup = NaN(numel(curTracksNApaxG1),1);
        curForceEarlySlopeGroup = NaN(numel(curTracksNApaxG1),1);
        curAmpOverallSlopeGroup = NaN(numel(curTracksNApaxG1),1);
        curAmpOverallSlopeUpErrGroup = NaN(numel(curTracksNApaxG1),1);
        curAmpOverallSlopeDownErrGroup = NaN(numel(curTracksNApaxG1),1);
        curForceLateSlopeGroup = NaN(numel(curTracksNApaxG1),1);

        periodFrames = 30;
        rr=0;
        for pp=1:numel(curTracksNApaxG1)
            rr=rr+1;
            curTrack = curTracksNApaxG1(pp);
            [~,curForceSlopeGroup(rr)] = regression((1:curTrack.lifeTime+1),curTrack.forceMag(curTrack.startingFrameExtra:curTrack.endingFrameExtra));
            curEndFrame = min(curTrack.startingFrameExtra+periodFrames-1,curTrack.endingFrame);
            curEarlyPeriod = curEndFrame - curTrack.startingFrameExtra+1;
            [~,curForceEarlySlopeGroup(rr)] = regression((1:curEarlyPeriod),curTrack.forceMag(curTrack.startingFrameExtra:curEndFrame));
            sF=curTrack.startingFrameExtraExtra; eF=curTrack.endingFrameExtraExtra;
            curEndFrame = min(sF+periodFrames-1,curTrack.endingFrame);
            curEarlyPeriod = curEndFrame - sF+1;
            [~,curVeryEarlyAmpSlopeGroup(rr)] = regression((1:curEarlyPeriod),curTrack.ampTotal(sF:curEndFrame));

            % There is amplitude that has rising phase but just stops in the middle
            % of course. (It means the the particle loses it's point-ness while
            % increasing the intensity. This does not help in characterizing true
            % G1 behavior
            % I'll filter these out with lateAmpSlope
            curStartFrame = max(eF-periodFrames,curTrack.startingFrame);
            curLatePeriod = curTrack.endingFrameExtraExtra - curStartFrame + 1;
            [~,curForceLateSlopeGroup(rr)] = regression((1:curLatePeriod),curTrack.ampTotal(curStartFrame:eF));
        end
        indFlatForce = curForceSlopeGroup<=0 & curForceEarlySlopeGroup<=0;
        validTrackID = (~indFlatForce & ~indFlatAmp) & curVeryEarlyAmpSlopeGroup>0 ...
            & curForceLateSlopeGroup<0;

        tracksG1fromG1_1=curTracksNApaxG1(validTrackID);

        %Potential G2 in G1
        ampSlopeG1 = arrayfun(@(x) x.ampSlope,curTracksNApaxG1);
        lifeTimeG1 = arrayfun(@(x) x.lifeTime,curTracksNApaxG1);
        ampEndingG1 = arrayfun(@(x) x.ampTotal(x.endingFrameExtra),curTracksNApaxG1);
        ampStartingG1 = arrayfun(@(x) x.ampTotal(x.startingFrameExtra),curTracksNApaxG1);
        idxIncreasingAmpG2inG1 = ampSlopeG1>0 & ampEndingG1>ampStartingG1;
        % idxLowInitForceG2= initForceG2<500;
        % maturing NA should have at least 3 min of lifetime
        thresLT_G2_frames_inG1 = 3*60/tInterval;
        idxLongLifeTimeG2inG1=lifeTimeG1>thresLT_G2_frames_inG1;
        idGroup2fromG1 = idxIncreasingAmpG2inG1 & idxLongLifeTimeG2inG1; %& idxLowInitForceG2 
        tracksG2fromG1_1 = curTracksNApaxG1(idGroup2fromG1);

        realTypesFromG1 = zeros(size(tracksG1fromG1_1));
        for kk=1:numel(tracksG1fromG1_1)
            tracksG1fromG1_1(kk)=readIntensityFromTracks(tracksG1fromG1_1(kk),imgMap,1,'extraLength',30,'movieData',curMD);
            tracksG1fromG1_1(kk)=readIntensityFromTracks(tracksG1fromG1_1(kk),tMap,2,'movieData',curMD);
            [firstIncreseTimeIntAgainstForce,forceTransmittingAll...
            ,firstIncreseTimeIntAll,firstIncreseTimeForce,bkgMaxIntAll,bkgMaxForce] =...
                calculateFirstIncreaseTimeTracks(tracksG1fromG1_1(kk),splineParamInit,preDetecFactor,tInterval);

            tracksG1fromG1_1(kk).forceTransmitting=forceTransmittingAll;
            tracksG1fromG1_1(kk).firstIncreseTimeInt=firstIncreseTimeIntAll;
            tracksG1fromG1_1(kk).firstIncreseTimeForce=firstIncreseTimeForce;
            tracksG1fromG1_1(kk).firstIncreseTimeIntAgainstForce=firstIncreseTimeIntAgainstForce;
            tracksG1fromG1_1(kk).bkgMaxInt=bkgMaxIntAll;
            tracksG1fromG1_1(kk).bkgMaxForce=bkgMaxForce;

            hSum=showSingleAdhesionTrackSummary(curMD,tracksG1fromG1_1(kk),imgMap,tMap,kk);
            hBcc=figure; plot(crossVariance(tracksG1fromG1_1(kk).ampTotal,tracksG1fromG1_1(kk).forceMag,13))
            [ cur_tLagBccPosG1,cur_tLagBccNegG1,cur_BccInLaterG1,cur_BccInPriorG1 ] ...
                = calculateBccPeakLag( tracksG1fromG1_1(kk).ampTotal, tracksG1fromG1_1(kk).forceMag , tracksG1fromG1_1(kk).startingFrameExtra,curMD.timeInterval_);        
            disp(['tLagPos = ' num2str(cur_tLagBccPosG1) ', BccLater = ' num2str(cur_BccInLaterG1) ' tLagNeg=' num2str(cur_tLagBccNegG1) ' BccPrior:' num2str(cur_BccInPriorG1) ])

            s=(inputdlg('What class is this NA? (1-9): '));
            try
                realTypesFromG1(kk)=str2num(s{1});
            catch
                realTypesFromG1(kk)=3;
            end
            close(hSum); close(hBcc);
            disp(['G1: ' num2str(kk) ' out of ' num2str(numel(tracksG1fromG1_1)) ' done.'])
        end
        tracksG2fromG1_2 = tracksG1fromG1_1(realTypesFromG1==2);
        tracksG1fromG1_2 = tracksG1fromG1_1(realTypesFromG1==1);

        % Now curTracksNApaxG2
        ampSlopeG2 = arrayfun(@(x) x.ampSlope,curTracksNApaxG2);
        initForceG2 = arrayfun(@(x) x.forceMag(x.startingFrame),curTracksNApaxG2);
        lifeTimeG2 = arrayfun(@(x) x.lifeTime,curTracksNApaxG2);
        ampEndingG2 = arrayfun(@(x) x.ampTotal(x.endingFrameExtra),curTracksNApaxG2);
        ampStartingG2 = arrayfun(@(x) x.ampTotal(x.startingFrameExtra),curTracksNApaxG2);
        idxIncreasingAmpG2 = ampSlopeG2>0 & ampEndingG2>ampStartingG2;
        % idxLowInitForceG2= initForceG2<500;
        % maturing NA should have at least 3 min of lifetime
        thresLT_G2_frames = 2*60/tInterval;
        idxLongLifeTimeG2=lifeTimeG2>thresLT_G2_frames;
        realG2pax = curForceTransmittingG2 & curStartingAmpG2<thresAmp & stateFA;
        idGroup2valid = idxIncreasingAmpG2 & idxLongLifeTimeG2 & realG2pax; %& idxLowInitForceG2 

        curTracksNApaxG2f = curTracksNApaxG2(idGroup2valid);

        realTypesFromG2 = zeros(numel(curTracksNApaxG2f),1);
    %     realTypesFromG2 = false(size(realG2pax));
        for kk=1:numel(curTracksNApaxG2f)
            curTracksNApaxG2f(kk)=readIntensityFromTracks(curTracksNApaxG2f(kk),imgMap,1,'extraLength',30,'movieData',curMD);
            curTracksNApaxG2f(kk)=readIntensityFromTracks(curTracksNApaxG2f(kk),tMap,2,'movieData',curMD);

            [firstIncreseTimeIntAgainstForce,forceTransmittingAll...
            ,firstIncreseTimeIntAll,firstIncreseTimeForce,bkgMaxIntAll,bkgMaxForce] =...
                calculateFirstIncreaseTimeTracks(curTracksNApaxG2f(kk),splineParamInit,preDetecFactor,tInterval);
            curTracksNApaxG2f(kk).forceTransmitting=forceTransmittingAll;
            curTracksNApaxG2f(kk).firstIncreseTimeInt=firstIncreseTimeIntAll;
            curTracksNApaxG2f(kk).firstIncreseTimeForce=firstIncreseTimeForce;
            curTracksNApaxG2f(kk).firstIncreseTimeIntAgainstForce=firstIncreseTimeIntAgainstForce;
            curTracksNApaxG2f(kk).bkgMaxInt=bkgMaxIntAll;
            curTracksNApaxG2f(kk).bkgMaxForce=bkgMaxForce;

            hSum=showSingleAdhesionTrackSummary(curMD,curTracksNApaxG2f(kk),imgMap,tMap,kk);
            hBcc=figure; plot(crossVariance(curTracksNApaxG2f(kk).ampTotal,curTracksNApaxG2f(kk).forceMag,13))
            [ cur_tLagBccPosG2,cur_tLagBccNegG2,cur_BccInLaterG2,cur_BccInPriorG2 ] ...
                = calculateBccPeakLag( curTracksNApaxG2f(kk).ampTotal, curTracksNApaxG2f(kk).forceMag , curTracksNApaxG2f(kk).startingFrameExtra,curMD.timeInterval_);        
            disp(['tLagPosG2:' num2str(cur_tLagBccPosG2) ', tLagBccNeg: ' num2str(cur_tLagBccNegG2) ', cur_BccInLaterG2: ' num2str(cur_BccInLaterG2) ' cur_BccInPriorG2:' num2str(cur_BccInPriorG2) ])
            s=(inputdlg('What class is this NA? (1-9): '));
            realTypesFromG2(kk)=str2num(s{1});
            close(hSum); close(hBcc);
        end
        tracksG2fromG2 = curTracksNApaxG2f(realTypesFromG2==2);
        tracksG1fromG2 = curTracksNApaxG2f(realTypesFromG2==1);

        curTracksNApaxG2f = [tracksG2fromG1_1; tracksG2fromG1_2; tracksG2fromG2];%[tracksG2fromG1; tracksG2fromG2];
        curTracksNApaxG1 = [tracksG1fromG1_2; tracksG1fromG2];%[tracksG1fromG1; tracksG1fromG2];
        tracksG2=curTracksNApaxG2f; tracksG1=curTracksNApaxG1;
        save([paxFolder{ii} filesep 'data' filesep 'tracksG2real.mat'],'tracksG2');
        save([paxFolder{ii} filesep 'data' filesep 'tracksG1real.mat'],'tracksG1');
    end
    % Get the time measurement
    tInterval = curMD.timeInterval_;
    [curTimeDelayG2] = calculateFirstIncreaseTimeTracks(curTracksNApaxG2f,splineParamInit,preDetecFactor,tInterval);
    if isempty(curTimeDelayG2)
        disp('No curTimeDelayG2')
    else
        disp(nanmedian(curTimeDelayG2))
    end
    [curTimeDelayG1] = calculateFirstIncreaseTimeTracks(curTracksNApaxG1,splineParamInit,preDetecFactor,tInterval);
    disp(nanmedian(curTimeDelayG1))

    InitTimeLagMaturePax{ii} = curTimeDelayG2;
    disp([paxFolder{ii} ' has been analyzed.'])
    
    tlagInitPax{ii} = curTimeDelayG1; %cur_tlagInitPax;

    cur_tlagPeakPax = load([paxFolder{ii} '/data/timePeaksG1.mat'],'peakTimeIntAgainstForceAll');
    cur_tlagPeakPax = cur_tlagPeakPax.peakTimeIntAgainstForceAll;
    tlagPeakPax{ii} = cur_tlagPeakPax;
    
    % Get the assembly rate and force growth rate
    [~,curEarlyAssmRateG1] = arrayfun(@(x) regression(1:3,x.ampTotal(x.startingFrameExtra:x.startingFrameExtra+2)),curTracksNApaxG1);
    [~,curEarlyAssmRateG2] = arrayfun(@(x) regression(1:5,x.ampTotal(x.startingFrameExtra:x.startingFrameExtra+4)),curTracksNApaxG2f);
    earlyAssmRateG1All{ii} = curEarlyAssmRateG1; %cur_tlagInitPax;
    earlyAssmRateG2All{ii} = curEarlyAssmRateG2; %cur_tlagInitPax;
    
    [~,earlyForceRateG1] = arrayfun(@(x) regression(1:3,x.forceMag(x.startingFrameExtra:x.startingFrameExtra+2)),curTracksNApaxG1);
    [~,earlyForceRateG2] = arrayfun(@(x) regression(1:3,x.forceMag(x.startingFrameExtra:x.startingFrameExtra+2)),curTracksNApaxG2f);
    earlyForceRateG1All{ii} = earlyForceRateG1; %cur_tlagInitPax;
    earlyForceRateG2All{ii} = earlyForceRateG2; %cur_tlagInitPax;

    
    % Bcc measurement - G1
    tLagBccNeg = NaN(numel(curTracksNApaxG1),1);
    BccInPrior = NaN(numel(curTracksNApaxG1),1);
    tLagBccPos = NaN(numel(curTracksNApaxG1),1);
    BccInLater = NaN(numel(curTracksNApaxG1),1);
    for pp=1:numel(curTracksNApaxG1)
        curTrack = curTracksNApaxG1(pp);
        if length(curTrack.forceMag)<length(curTrack.ampTotal)
            tempForce=curTrack.forceMag;
            curTrack.forceMag=curTrack.ampTotal;
            curTrack.forceMag(1:length(tempForce))=tempForce;
        end
        curBcc = crossVariance(curTrack.ampTotal,curTrack.forceMag,11);
        %Will find the peak against curTrack.startingFrameExtra to the
        %front and back
        maxBcc=nanmax(curBcc);
        indCands=locmax1d(curBcc,5);
        % Exclude the locmaxes at the ends
        firstInd = find(~isnan(curBcc),1);
        lastInd = find(~isnan(curBcc),1,'last');
        indCands(ismember(indCands,[firstInd,lastInd]))=[];
        % if the loc max is close to curTrack.startingFrameExtra and the
        % actual value is above 0.8*maxBcc, we'll admit the time as
        % BccShift time
        
        % Sort based on proximity to curTrack.startingFrameExtra
        SFE = curTrack.startingFrameExtra;
        frameDiff = indCands-SFE;
        
        if any(frameDiff<0)
            indNegDiff= indCands(frameDiff<0);
            for curIndCand = indNegDiff'
                if curBcc(curIndCand)>0.5*maxBcc
                    tLagBccNeg(pp) = frameDiff(indCands==curIndCand)*curMD.timeInterval_;
                    BccInPrior(pp) = curBcc(curIndCand);
                    break
                end
            end
        end
        indPosDiff= indCands(frameDiff>=0);

        for curIndCand = indPosDiff'
            if curBcc(curIndCand)>0.7*maxBcc
                tLagBccPos(pp) = frameDiff(indCands==curIndCand)*curMD.timeInterval_;
                BccInLater(pp) = curBcc(curIndCand);
                break
            end
        end
    end
    tLagBccPosG1All{ii} = tLagBccPos; 
    tLagBccNegG1All{ii} = tLagBccNeg; 
    BccInLaterG1All{ii} = BccInLater; 
    BccInPriorG1All{ii} = BccInPrior; 

    % Bcc measurement - G2
    tLagBccNegG2 = NaN(numel(curTracksNApaxG2f),1);
    BccInPriorG2 = NaN(numel(curTracksNApaxG2f),1);
    tLagBccPosG2 = NaN(numel(curTracksNApaxG2f),1);
    BccInLaterG2 = NaN(numel(curTracksNApaxG2f),1);
    for pp=1:numel(curTracksNApaxG2f)
        curTrack = curTracksNApaxG2f(pp);
        if length(curTrack.forceMag)<length(curTrack.ampTotal)
            tempForce=curTrack.forceMag;
            curTrack.forceMag=curTrack.ampTotal;
            curTrack.forceMag(1:length(tempForce))=tempForce;
        end
        curBcc = crossVariance(curTrack.ampTotal,curTrack.forceMag,11);
        %Will find the peak against curTrack.startingFrameExtra to the
        %front and back
        maxBcc=nanmax(curBcc);
        indCands=locmax1d(curBcc,5);
        % Exclude the locmaxes at the ends
        firstInd = find(~isnan(curBcc),1);
        lastInd = find(~isnan(curBcc),1,'last');
        indCands(ismember(indCands,[firstInd,lastInd]))=[];
        % if the loc max is close to curTrack.startingFrameExtra and the
        % actual value is above 0.8*maxBcc, we'll admit the time as
        % BccShift time
        
        % Sort based on proximity to curTrack.startingFrameExtra
        SFE = curTrack.startingFrameExtra;
        frameDiff = indCands-SFE;
        
        if any(frameDiff<0)
            indNegDiff= indCands(frameDiff<0);
            for curIndCand = indNegDiff'
                if curBcc(curIndCand)>0.4*maxBcc
                    tLagBccNegG2(pp) = frameDiff(indCands==curIndCand)*curMD.timeInterval_;
                    BccInPriorG2(pp) = curBcc(curIndCand);
                    break
                end
            end
        end
        indPosDiff= indCands(frameDiff>=0);

        for curIndCand = indPosDiff'
            if curBcc(curIndCand)>0.5*maxBcc
                tLagBccPosG2(pp) = frameDiff(indCands==curIndCand)*curMD.timeInterval_;
                BccInLaterG2(pp) = curBcc(curIndCand);
                break
            end
        end
    end
    tLagBccPosG2All{ii} = tLagBccPosG2; 
    tLagBccNegG2All{ii} = tLagBccNegG2; 
    BccInLaterG2All{ii} = BccInLaterG2; 
    BccInPriorG2All{ii} = BccInPriorG2; 
    
    clear curTracksNApaxG2f curTracksNApaxG1 imgMap tMap
end
totalInitTLagPax=[tlagInitPax{1}; tlagInitPax{2}; tlagInitPax{3}; tlagInitPax{4}];
totalPeakTLagPax=[tlagPeakPax{1}; tlagPeakPax{2}; tlagPeakPax{3}; tlagPeakPax{4}];
totalInitTimeLagMaturePax =[InitTimeLagMaturePax{1}; InitTimeLagMaturePax{2}; 
    InitTimeLagMaturePax{3}; InitTimeLagMaturePax{4}];

tLagBccPosG1Pax =cell2mat(tLagBccPosG1All(1:numPaxMovies));
tLagBccNegG1Pax =cell2mat(tLagBccNegG1All(1:numPaxMovies));
tLagBccPosG2Pax =cell2mat(tLagBccPosG2All(1:numPaxMovies));
tLagBccNegG2Pax =cell2mat(tLagBccNegG2All(1:numPaxMovies));
BccInPriorG1Pax =cell2mat(BccInPriorG1All(1:numPaxMovies));
BccInLaterG1Pax =cell2mat(BccInLaterG1All(1:numPaxMovies));
BccInPriorG2Pax =cell2mat(BccInPriorG2All(1:numPaxMovies));
BccInLaterG2Pax =cell2mat(BccInLaterG2All(1:numPaxMovies));

totalEarlyAssmRateG1Pax=[earlyAssmRateG1All{1}; earlyAssmRateG1All{2}; earlyAssmRateG1All{3}; earlyAssmRateG1All{4}];
totalEarlyAssmRateG2Pax=[earlyAssmRateG2All{1}; earlyAssmRateG2All{2}; earlyAssmRateG2All{3}; earlyAssmRateG2All{4}];
totalEarlyForceRateG1Pax=[earlyForceRateG1All{1}; earlyForceRateG1All{2}; earlyForceRateG1All{3}; earlyForceRateG1All{4}];
totalEarlyForceRateG2Pax=[earlyForceRateG2All{1}; earlyForceRateG2All{2}; earlyForceRateG2All{3}; earlyForceRateG2All{4}];
clear earlyAssmRateG2All earlyAssmRateG1All earlyForceRateG1All earlyForceRateG2All
%% Box Plot!
%% filtering with high and low limit (-60 and 60)
lowTLag=-80; hiTLag=60;
totalInitTLagTal2 = totalInitTLagTal(totalInitTLagTal>lowTLag & totalInitTLagTal<hiTLag);
totalInitTLagVin2 = totalInitTLagVin(totalInitTLagVin>lowTLag & totalInitTLagVin<hiTLag);
totalInitTLagPax2 = totalInitTLagPax(totalInitTLagPax>lowTLag & totalInitTLagPax<hiTLag);
%% Actual Fig 3
% making it to matrix
[lengthLongest,iLongest]=max([length(totalInitTLagPax2),length(totalInitTLagVin2),length(totalInitTLagTal2)]);
% it was vinculin that was longest
matrixPaxTalVin = NaN(lengthLongest,3);
matrixPaxTalVin(1:length(totalInitTLagPax2),1) = totalInitTLagPax2;
matrixPaxTalVin(1:length(totalInitTLagVin2),2) = totalInitTLagVin2;
matrixPaxTalVin(1:length(totalInitTLagTal2),3) = totalInitTLagTal2;
boxWidth=0.5;

axes('Position',[40/600,125/600,140/600,65/600])
boxplot(matrixPaxTalVin,'orientation','horizontal','whisker',0.5,'notch','on',...
    'labels',{['Pax (N=' num2str(length(totalInitTLagPax2)) ')'],...
    ['Vin (N=' num2str(length(totalInitTLagVin2)) ')'],...
    ['Tal (N=' num2str(length(totalInitTLagTal2)) ')']},...
    'symbol','','widths',boxWidth,'jitter',1,'colors','k')
txt = findobj(gca,'Type','text');
set(txt(3:end),'VerticalAlignment', 'Middle');
set(findobj(gca,'LineStyle','--'),'LineStyle','-')
set(findobj(gca,'tag','Median'),'LineWidth',2)
xlim([-55 12])
xlabel('Time lag in initial rise (s)')
title('non-maturing adhesions (G1)')
set(gca,'FontSize',7)
%% statistical testing for inital lag
% normality test
[hNormalPax,pNormalPax] = kstest(totalInitTLagPax2);
% hNormalPax = 1
% pNormalPax = 7.8863e-242
% So it's not normal! 
[pPaxVsVin,hPaxVsVin] = ranksum(totalInitTLagPax2,totalInitTLagVin2);
[pPaxVsTal,hPaxVsTal] = ranksum(totalInitTLagPax2,totalInitTLagTal2);
[pVinVsTal,hVinVsTal] = ranksum(totalInitTLagVin2,totalInitTLagTal2);

%% p-values
% findobj(gca,'tag','Median')
% findobj(gca,'tag','Outlier')
% get(findobj(gca,'Group','Pax'),'OutlierValue')
% [x y]=ginput(1);
medPaxInit=nanmedian(totalInitTLagPax);
medVinInit=nanmedian(totalInitTLagVin);
medTalInit=nanmedian(totalInitTLagTal);

disp([cellfun(@(x) nanmedian(x),tlagInitPax); cellfun(@(x) length(x),tlagInitPax)])
disp([cellfun(@(x) nanmedian(x),tlagInitVin); cellfun(@(x) length(x),tlagInitVin)])
disp([cellfun(@(x) nanmedian(x),tlagInitTal); cellfun(@(x) length(x),tlagInitTal)])

text(14,1.5,['p=' num2str(pPaxVsVin,2)]) %pax vs vin
text(14,2.5,['p=' num2str(pVinVsTal,2)]) %tal vs vin
text(14,2,['p=' num2str(pPaxVsTal,2)]) %pax vs tal

% median values
text(medPaxInit,1.5,num2str(medPaxInit,'%10.0f'),'HorizontalAlignment','center') %pax vs vin
text(medVinInit,2.5,num2str(medVinInit,'%10.0f'),'HorizontalAlignment','center') %tal vs vin
text(medTalInit,3.5,num2str(medTalInit,'%10.0f'),'HorizontalAlignment','center') %pax vs tal

set(findobj(gca,'Type','Text'),'FontSize',7)
set(gca,'FontSize',7)
%% boxplot for Bcc
axes('Position', [280/600 125/600 110/600 65/600]); % 
% axes('Position', [510/600 325/600 80/600 65/600]); % 
figure
boxPlotCellArray({tLagBccPosG1Talin,tLagBccPosG1Vin,tLagBccPosG1Pax},{'Tal','Vin','Pax'})
figure
boxPlotCellArray({tLagBccPosG2Talin,tLagBccPosG2Vin,tLagBccPosG2Pax},{'Tal','Vin','Pax'})
figure
boxPlotCellArray({tLagBccNegG1Talin,tLagBccNegG1Vin,tLagBccNegG1Pax},{'Tal','Vin','Pax'})
figure
boxPlotCellArray({tLagBccNegG2Talin,tLagBccNegG2Vin,tLagBccNegG2Pax},{'Tal','Vin','Pax'})

%% saving
% title('Time lag in initial rises')%% saving
print('-depsc', '-loose', [f1path filesep 'Fig1Middle_spatial_corr2.eps'])
savefig([f1path filesep 'Fig1Middle_spatial_corr2.fig'])
print('-dtiff', '-loose', '-r300', [f1path filesep 'Fig1Middle_spatial_corr2.tif'])
% close gcf
save([f1path filesep 'fig1.mat'])
%% saving
mkdir(fS2path)
print('-depsc', '-loose', [fS2path filesep 'FigS2AllGroupTimeSeries.eps'])
savefig([fS2path filesep 'FigS2AllGroupTimeSeries.fig'])
print('-dtiff', '-loose', '-r300', [fS2path filesep 'FigS2AllGroupTimeSeries.tif'])
close 
%% Fig1 later: Feature difference
%% Loading all tracks!
for ii=1:4
    tic;
    curtracksNApax = load([paxFolder{ii} filesep 'data' filesep 'tracksNA.mat'],'tracksNA');
    tracksNApaxAll{ii} = curtracksNApax.tracksNA; toc

    idGroupsPaxStruct = load([paxFolder{ii} filesep 'data' filesep 'idsClassified.mat'],'idGroup1','idGroup2','idGroup3',...
    'idGroup4','idGroup5','idGroup6','idGroup7','idGroup8','idGroup9');
    idGroupsPaxAll{ii}{1} = idGroupsPaxStruct.idGroup1;
    idGroupsPaxAll{ii}{2} = idGroupsPaxStruct.idGroup2;
    idGroupsPaxAll{ii}{3} = idGroupsPaxStruct.idGroup3;
    idGroupsPaxAll{ii}{4} = idGroupsPaxStruct.idGroup4;
    idGroupsPaxAll{ii}{5} = idGroupsPaxStruct.idGroup5;
    idGroupsPaxAll{ii}{6} = idGroupsPaxStruct.idGroup6;
    idGroupsPaxAll{ii}{7} = idGroupsPaxStruct.idGroup7;
    idGroupsPaxAll{ii}{8} = idGroupsPaxStruct.idGroup8;
    idGroupsPaxAll{ii}{9} = idGroupsPaxStruct.idGroup9;
end

for ii=1:6
    tic;
    curtracksNAvin = load([vinFolder{ii} filesep 'data' filesep 'tracksNA.mat'],'tracksNA');
    tracksNAvinAll{ii} = curtracksNAvin.tracksNA; toc

    idGroupsVinStruct = load([vinFolder{ii} filesep 'data' filesep 'idsClassified.mat'],'idGroup1','idGroup2','idGroup3',...
    'idGroup4','idGroup5','idGroup6','idGroup7','idGroup8','idGroup9');
    idGroupsVinAll{ii}{1} = idGroupsVinStruct.idGroup1;
    idGroupsVinAll{ii}{2} = idGroupsVinStruct.idGroup2;
    idGroupsVinAll{ii}{3} = idGroupsVinStruct.idGroup3;
    idGroupsVinAll{ii}{4} = idGroupsVinStruct.idGroup4;
    idGroupsVinAll{ii}{5} = idGroupsVinStruct.idGroup5;
    idGroupsVinAll{ii}{6} = idGroupsVinStruct.idGroup6;
    idGroupsVinAll{ii}{7} = idGroupsVinStruct.idGroup7;
    idGroupsVinAll{ii}{8} = idGroupsVinStruct.idGroup8;
    idGroupsVinAll{ii}{9} = idGroupsVinStruct.idGroup9;
end
for ii=1:1
    tic;
    curtracksNAtal = load([talFolder{ii} filesep 'data' filesep 'tracksNA.mat'],'tracksNA');
    tracksNAtalAll{ii} = curtracksNAtal.tracksNA; toc

    idGroupsTalStruct = load([talFolder{ii} filesep 'data' filesep 'idsClassified.mat'],'idGroup1','idGroup2','idGroup3',...
    'idGroup4','idGroup5','idGroup6','idGroup7','idGroup8','idGroup9');
    idGroupsTalAll{ii}{1} = idGroupsTalStruct.idGroup1;
    idGroupsTalAll{ii}{2} = idGroupsTalStruct.idGroup2;
    idGroupsTalAll{ii}{3} = idGroupsTalStruct.idGroup3;
    idGroupsTalAll{ii}{4} = idGroupsTalStruct.idGroup4;
    idGroupsTalAll{ii}{5} = idGroupsTalStruct.idGroup5;
    idGroupsTalAll{ii}{6} = idGroupsTalStruct.idGroup6;
    idGroupsTalAll{ii}{7} = idGroupsTalStruct.idGroup7;
    idGroupsTalAll{ii}{8} = idGroupsTalStruct.idGroup8;
    idGroupsTalAll{ii}{9} = idGroupsTalStruct.idGroup9;
end
%% Filtering each idGroups based on obvious criteria
% First criterion: edge & adhesion shouldn't retract too much in G1 and G2
% Second: edge & adhesion shouldn't retract in G3
for ii=1:4
    curEdgeAdvG1=arrayfun(@(x) mean(x.edgeAdvanceDist),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{1}));
    curAdhAdvG1=arrayfun(@(x) mean(x.advanceDist),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{1}));
    indG1RetractTooMuch=(curEdgeAdvG1<-2 & curEdgeAdvG1<curAdhAdvG1) | curEdgeAdvG1<-10;
    curIndG1=find(idGroupsPaxAll{ii}{1});
    idGroupsPaxAll{ii}{1}(curIndG1(indG1RetractTooMuch))=false;
    
    curEdgeAdvG2=arrayfun(@(x) mean(x.edgeAdvanceDist),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{2}));
    curAdhAdvG2=arrayfun(@(x) mean(x.advanceDist),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{2}));
    indG2RetractTooMuch=(curEdgeAdvG2<-2 & curEdgeAdvG2<curAdhAdvG2) | curEdgeAdvG2<-10;
    curIndG2=find(idGroupsPaxAll{ii}{2});
    idGroupsPaxAll{ii}{2}(curIndG2(indG2RetractTooMuch))=false;

    curEdgeAdvG3=arrayfun(@(x) mean(x.edgeAdvanceDist),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{3}));
    curAdhAdvG3=arrayfun(@(x) mean(x.advanceDist),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{3}));
    indG3RetractTooMuch=curEdgeAdvG3<0  | curAdhAdvG3<-1;
    curIndG3=find(idGroupsPaxAll{ii}{3});
    idGroupsPaxAll{ii}{3}(curIndG3(indG3RetractTooMuch))=false;
end

for ii=1:6
    curEdgeAdvG1=arrayfun(@(x) mean(x.edgeAdvanceDist),tracksNAvinAll{ii}(idGroupsVinAll{ii}{1}));
    curAdhAdvG1=arrayfun(@(x) mean(x.advanceDist),tracksNAvinAll{ii}(idGroupsVinAll{ii}{1}));
    indG1RetractTooMuch=(curEdgeAdvG1<-2 & curEdgeAdvG1<curAdhAdvG1) | curEdgeAdvG1<-10;
    curIndG1=find(idGroupsVinAll{ii}{1});
    idGroupsVinAll{ii}{1}(curIndG1(indG1RetractTooMuch))=false;
    
    curEdgeAdvG2=arrayfun(@(x) mean(x.edgeAdvanceDist),tracksNAvinAll{ii}(idGroupsVinAll{ii}{2}));
    curAdhAdvG2=arrayfun(@(x) mean(x.advanceDist),tracksNAvinAll{ii}(idGroupsVinAll{ii}{2}));
    indG2RetractTooMuch=(curEdgeAdvG2<-2 & curEdgeAdvG2<curAdhAdvG2) | curEdgeAdvG2<-10;
    curIndG2=find(idGroupsVinAll{ii}{2});
    idGroupsVinAll{ii}{2}(curIndG2(indG2RetractTooMuch))=false;

    curEdgeAdvG3=arrayfun(@(x) mean(x.edgeAdvanceDist),tracksNAvinAll{ii}(idGroupsVinAll{ii}{3}));
    curAdhAdvG3=arrayfun(@(x) mean(x.advanceDist),tracksNAvinAll{ii}(idGroupsVinAll{ii}{3}));
    indG3RetractTooMuch=curEdgeAdvG3<0  | curAdhAdvG3<-1;
    curIndG3=find(idGroupsVinAll{ii}{3});
    idGroupsVinAll{ii}{3}(curIndG3(indG3RetractTooMuch))=false;
end

for ii=1:6
    curEdgeAdvG1=arrayfun(@(x) mean(x.edgeAdvanceDist),tracksNAtal{ii}(idGroupsTalAll{ii}{1}));
    curAdhAdvG1=arrayfun(@(x) mean(x.advanceDist),tracksNAtal{ii}(idGroupsTalAll{ii}{1}));
    indG1RetractTooMuch=(curEdgeAdvG1<-2 & curEdgeAdvG1<curAdhAdvG1) | curEdgeAdvG1<-10;
    curIndG1=find(idGroupsTalAll{ii}{1});
    idGroupsTalAll{ii}{1}(curIndG1(indG1RetractTooMuch))=false;
    
    curEdgeAdvG2=arrayfun(@(x) mean(x.edgeAdvanceDist),tracksNAtal{ii}(idGroupsTalAll{ii}{2}));
    curAdhAdvG2=arrayfun(@(x) mean(x.advanceDist),tracksNAtal{ii}(idGroupsTalAll{ii}{2}));
    indG2RetractTooMuch=(curEdgeAdvG2<-2 & curEdgeAdvG2<curAdhAdvG2) | curEdgeAdvG2<-10;
    curIndG2=find(idGroupsTalAll{ii}{2});
    idGroupsTalAll{ii}{2}(curIndG2(indG2RetractTooMuch))=false;

    curEdgeAdvG3=arrayfun(@(x) mean(x.edgeAdvanceDist),tracksNAtal{ii}(idGroupsTalAll{ii}{3}));
    curAdhAdvG3=arrayfun(@(x) mean(x.advanceDist),tracksNAtal{ii}(idGroupsTalAll{ii}{3}));
    indG3RetractTooMuch=curEdgeAdvG3<0  | curAdhAdvG3<-1;
    curIndG3=find(idGroupsTalAll{ii}{3});
    idGroupsTalAll{ii}{3}(curIndG3(indG3RetractTooMuch))=false;
end
%% intensity and force profiles - vinculin
% axes('Position', [130/600 330/600 65/600 65/600]); % 
% % plotIntensityForce(tracksNAvin(idGroupsVin5{1} | idGroupsVin5{7}),[],false,false,'UseCurrentAxis',true,'Source',{'ampTotal'},'plotCohorts',true,'tInterval',MDvin.timeInterval_,'prePostFrames',10)
% % set(findobj(gca,'Type','Text'),'FontSize',7)
% axes('Position', [130/600 230/600 65/600 65/600]); % 
% plotIntensityForce(tracksNAvin(idGroupsVin5{1} | idGroupsVin5{7}),[],false,false,'UseCurrentAxis',true,'Source',{'forceMag'},'plotCohorts',true,'tInterval',MDvin.timeInterval_,'prePostFrames',10)
%% only F-T - vinculin
%% vinculin - G2
axes('Position', [230/600 330/600 65/600 65/600]); % 
plotIntensityForce(tracksNAvin((idGroupsVin5{2}) ),[],false,false,'UseCurrentAxis',true,...
    'Source',{'ampTotal'},'plotCohorts',true,'tInterval',MDvin.timeInterval_,'prePostFrames',10,'yNormalization',false)
xlim([-5 220])
ylim([230 680])
axes('Position', [230/600 230/600 65/600 65/600]); % 
plotIntensityForce(tracksNAvin((idGroupsVin5{2})),[],false,false,'UseCurrentAxis',true,...
    'Source',{'forceMag'},'plotCohorts',true,'tInterval',MDvin.timeInterval_,'prePostFrames',10,'yNormalization',false)
xlim([-5 220])
ylim([5 305])
%% intensity and force profiles - paxillin
% axes('Position', [130/600 130/600 65/600 65/600]); % 
% plotIntensityForce(tracksNApax(idGroupsPax{1} | idGroupsPax{7}),[],false,false,'UseCurrentAxis',true,'Source',{'ampTotal'},'plotCohorts',true,'tInterval',MDpax.timeInterval_,'prePostFrames',10)
% % set(findobj(gca,'Type','Text'),'FontSize',7)
% axes('Position', [130/600 30/600 65/600 65/600]); % 
% plotIntensityForce(tracksNApax(idGroupsPax{1} | idGroupsPax{7}),[],false,false,'UseCurrentAxis',true,'Source',{'forceMag'},'plotCohorts',true,'tInterval',MDpax.timeInterval_,'prePostFrames',10)
%% paxillin - G2 time series
axes('Position', [230/600 130/600 65/600 65/600]); % 
plotIntensityForce(tracksNApax((idGroupsPax{2})),[],false,false,'UseCurrentAxis',true,...
    'Source',{'ampTotal'},'plotCohorts',true,'tInterval',MDpax.timeInterval_,'prePostFrames',10,'yNormalization',false)
xlim([-5 140])
ylim([300 750])
axes('Position', [230/600 30/600 65/600 65/600]); % 
plotIntensityForce(tracksNApax((idGroupsPax{2})),[],false,false,'UseCurrentAxis',true,...
    'Source',{'forceMag'},'plotCohorts',true,'tInterval',MDpax.timeInterval_,'prePostFrames',10,'yNormalization',false)
xlim([-5 130])
ylim([10 150])



%% Looking at each cell
tlagBound = 300; 
filteredPeakTLagPax=totalPeakTLagPax(totalPeakTLagPax>-tlagBound & totalPeakTLagPax<tlagBound);
filteredPeakTLagVin=totalPeakTLagVin(totalPeakTLagVin>-tlagBound & totalPeakTLagVin<tlagBound);
filteredPeakTLagTal=totalPeakTLagTal(totalPeakTLagTal>-tlagBound & totalPeakTLagTal<tlagBound);

medPaxPeak=median(totalPeakTLagPax);
medVinPeak=median(totalPeakTLagVin);
medTalPeak=median(totalPeakTLagTal);
% median(filteredPeakTLagPax)
% median(filteredPeakTLagVin)
% median(filteredPeakTLagTal)

% mean(totalPeakTLagPax)
% mean(totalPeakTLagVin)
% mean(totalPeakTLagTal)

disp([cellfun(@(x) nanmedian(x),tlagPeakPax); cellfun(@(x) length(x),tlagPeakPax)])
disp([cellfun(@(x) nanmedian(x),tlagPeakVin); cellfun(@(x) length(x),tlagPeakVin)])
disp([cellfun(@(x) nanmedian(x),tlagPeakTal); cellfun(@(x) length(x),tlagPeakTal)])

%% statistical testing for peak lag
% normality test
[hNormalPaxPeak,pNormalPaxPeak] = kstest(filteredPeakTLagPax);
[hNormalVinPeak,pNormalVinPeak] = kstest(filteredPeakTLagVin);
[hNormalTalPeak,pNormalTalPeak] = kstest(filteredPeakTLagTal);
% hNormalPax = 1
% pNormalPax = 7.8863e-242
% So it's not normal! 
[pPaxVsVinPeak,hPaxVsVinPeak] = ranksum(filteredPeakTLagPax,filteredPeakTLagVin);
[pPaxVsTalPeak,hPaxVsTalPeak] = ranksum(filteredPeakTLagPax,filteredPeakTLagTal);
[pVinVsTalPeak,hVinVsTalPeak] = ranksum(filteredPeakTLagVin,filteredPeakTLagTal);

[pPaxVsVinPeak2,hPaxVsVinPeak2] = ranksum(totalPeakTLagPax,totalPeakTLagVin);
[pPaxVsTalPeak2,hPaxVsTalPeak2] = ranksum(totalPeakTLagPax,totalPeakTLagTal);
[pVinVsTalPeak2,hVinVsTalPeak2] = ranksum(totalPeakTLagVin,totalPeakTLagTal);



%% Time lag in G2
%% Maturing adhesions - t_init - loading
% splineParamInit=0.1;
% preDetecFactor = 0.01;
% for ii=1:4
%     try
% %         curInitTimeLagMature = load([paxFolder{ii} '/data/t_initForMaturingAdhesion.mat']);
% %         curInitTimeLagMaturePax = curInitTimeLagMature.firstIncreseTimeIntAgainstForceAllG2;
% %         InitTimeLagMaturePax{ii} = curInitTimeLagMaturePax;
% 
%         % try to filter G2 based on starting amplitude, peaks and slopes...
%         curTracksNApaxG2 = load([paxFolder{ii} filesep 'data' filesep 'tracksG2.mat'],'tracksG2');
%         curTracksNApaxG2 = curTracksNApaxG2.tracksG2; 
%         curForceTransmittingG2 = arrayfun(@(x) ~isempty(x.firstIncreseTimeIntAgainstForce), curTracksNApaxG2);
%         curStartingAmpG2 = arrayfun(@(x) x.ampTotal(x.startingFrameExtra), curTracksNApaxG2);
% 
%         curTracksNApaxG1 = load([paxFolder{ii} filesep 'data' filesep 'tracksG1.mat'],'tracksG1');
%         curTracksNApaxG1 = curTracksNApaxG1.tracksG1; 
%         curTimeDelay = arrayfun(@(x) x.firstIncreseTimeIntAgainstForce, curTracksNApaxG1,'UniformOutput',false);
%         curStartingAmpG1 = arrayfun(@(x) x.ampTotal(x.startingFrameExtra), curTracksNApaxG1);
%         thresAmp = mean(curStartingAmpG1)+3*std(curStartingAmpG1);
%         
%         % filtering for G2
%         stateFA = arrayfun(@(x) any(cellfun(@(y) strcmp(y,'FC') || strcmp(y,'FA'), x.state)), curTracksNApaxG2);
%         realG2pax = curForceTransmittingG2 & curStartingAmpG2<thresAmp & stateFA;
%         curTimeDelay = arrayfun(@(x) x.firstIncreseTimeIntAgainstForce, curTracksNApaxG2(realG2pax));
%         tInterval = MDpax.timeInterval_;
%         [~,curTimeDelayG2] = calculateFirstIncreaseTimeTracks(curTracksNAvinG2(realG2pax),splineParamInit,preDetecFactor,tInterval);
%         disp(nanmedian(curTimeDelayG2))
%         [~,curTimeDelayG1] = calculateFirstIncreaseTimeTracks(curTracksNAvinG1,splineParamInit,preDetecFactor,tInterval);
%         disp(nanmedian(curTimeDelayG1))
%         
%         InitTimeLagMaturePax{ii} = curTimeDelayG2;
%         clear curTracksNApaxG2 curTracksNApaxG1
%         disp([paxFolder{ii} ' has been analyzed for t_initForMaturingAdhesion.'])
%     catch
%         disp([paxFolder{ii} ' has not been analyzed for t_initForMaturingAdhesion.'])
%     end
% end
% totalInitTimeLagMaturePax =[InitTimeLagMaturePax{1}; InitTimeLagMaturePax{2}; 
%     InitTimeLagMaturePax{3}; InitTimeLagMaturePax{4}];

% for ii=1:6
%     try
% %         curInitTimeLagMature = load([vinFolder{ii} '/data/t_initForMaturingAdhesion.mat']);
% %         curInitTimeLagMatureVin = curInitTimeLagMature.firstIncreseTimeIntAgainstForceAllG2;
% %         InitTimeLagMatureVin{ii} = curInitTimeLagMatureVin;
% 
%         % try to filter G2 based on starting amplitude, peaks and slopes...
%         curTracksNAvinG2 = load([vinFolder{ii} filesep 'data' filesep 'tracksG2.mat'],'tracksG2');
%         curTracksNAvinG2 = curTracksNAvinG2.tracksG2; 
%         curForceTransmittingG2 = arrayfun(@(x) ~isempty(x.firstIncreseTimeIntAgainstForce), curTracksNAvinG2);
%         curStartingAmpG2 = arrayfun(@(x) x.ampTotal(x.startingFrameExtra), curTracksNAvinG2);
% 
%         curTracksNAvinG1 = load([vinFolder{ii} filesep 'data' filesep 'tracksG1.mat'],'tracksG1');
%         curTracksNAvinG1 = curTracksNAvinG1.tracksG1; 
%         curTimeDelay = arrayfun(@(x) x.firstIncreseTimeIntAgainstForce, curTracksNAvinG1,'UniformOutput',false);
%         curStartingAmpG1 = arrayfun(@(x) x.ampTotal(x.startingFrameExtra), curTracksNAvinG1);
%         thresAmp = mean(curStartingAmpG1)+3*std(curStartingAmpG1);
%         
%         % filtering for G2
%         stateFA = arrayfun(@(x) any(cellfun(@(y) strcmp(y,'FC') || strcmp(y,'FA'), x.state)), curTracksNAvinG2);
%         realG2vin = curForceTransmittingG2 & curStartingAmpG2<thresAmp & stateFA;
%         curTimeDelay = arrayfun(@(x) x.firstIncreseTimeIntAgainstForce, curTracksNAvinG2(realG2vin));
%         tInterval = MDvin.timeInterval_;
%         [~,curTimeDelayG2] = calculateFirstIncreaseTimeTracks(curTracksNAvinG2(realG2vin),splineParamInit,preDetecFactor,tInterval);
%         disp(nanmedian(curTimeDelayG2))
%         [~,curTimeDelayG1] = calculateFirstIncreaseTimeTracks(curTracksNAvinG1,splineParamInit,preDetecFactor,tInterval);
%         disp(nanmedian(curTimeDelayG1))
% %         curTimeDelay = arrayfun(@(x) x.firstIncreseTimeIntAgainstForce, curTracksNAvinG2(curForceTransmittingG2 & curStartingAmpG2<thresAmp));
%         InitTimeLagMatureVin{ii} = curTimeDelayG2;
%         clear curTracksNAvinG2 curTracksNAvinG1
%         disp([vinFolder{ii} ' has been analyzed for t_initForMaturingAdhesion.'])
%     catch
%         disp([vinFolder{ii} ' has not been analyzed for t_initForMaturingAdhesion.'])
%     end
% end
% totalInitTimeLagMatureVin =[InitTimeLagMatureVin{1}; InitTimeLagMatureVin{2}; 
%     InitTimeLagMatureVin{3}; InitTimeLagMatureVin{4}; InitTimeLagMatureVin{5}; InitTimeLagMatureVin{6}];

%% talin only...
% for ii=1:1
%     try
% %         curInitTimeLagMature = load([talFolder{ii} '/data/t_initForMaturingAdhesion.mat']);
% %         curInitTimeLagMatureTal = curInitTimeLagMature.firstIncreseTimeIntAgainstForceAllG2;
% %         InitTimeLagMatureTal{ii} = curInitTimeLagMatureTal;
%         
%         % try to filter G2 based on starting amplitude, peaks and slopes...
%         curTracksNAtalG2 = load([talFolder{ii} filesep 'data' filesep 'tracksG2.mat'],'tracksG2');
%         curTracksNAtalG2 = curTracksNAtalG2.tracksG2; 
%         curForceTransmittingG2 = arrayfun(@(x) ~isempty(x.firstIncreseTimeIntAgainstForce), curTracksNAtalG2);
%         curStartingAmpG2 = arrayfun(@(x) x.ampTotal(x.startingFrameExtra), curTracksNAtalG2);
% 
%         curTracksNAtalG1 = load([talFolder{ii} filesep 'data' filesep 'tracksG1.mat'],'tracksG1');
%         curTracksNAtalG1 = curTracksNAtalG1.tracksG1; 
%         curTimeDelay = arrayfun(@(x) x.firstIncreseTimeIntAgainstForce, curTracksNAtalG1,'UniformOutput',false);
%         curStartingAmpG1 = arrayfun(@(x) x.ampTotal(x.startingFrameExtra), curTracksNAtalG1);
%         thresAmp = mean(curStartingAmpG1)+3*std(curStartingAmpG1);
%         
%         % filtering for G2
%         stateFA = arrayfun(@(x) any(cellfun(@(y) strcmp(y,'FC') || strcmp(y,'FA'), x.state)), curTracksNAtalG2);
%         realG2tal = curForceTransmittingG2 & curStartingAmpG2<thresAmp & stateFA;
% %         curTimeDelay = arrayfun(@(x) x.firstIncreseTimeIntAgainstForce, curTracksNAtalG2(realG2tal));
%         tInterval = MDtal.timeInterval_;
%         [~,curTimeDelayG2] = calculateFirstIncreaseTimeTracks(curTracksNAtalG2(realG2tal),splineParamInit,preDetecFactor,tInterval);
%         disp(nanmedian(curTimeDelayG2))
%         [~,curTimeDelayG1] = calculateFirstIncreaseTimeTracks(curTracksNAtalG1,splineParamInit,preDetecFactor,tInterval);
%         disp(nanmedian(curTimeDelayG1))
% %         curTimeDelay = arrayfun(@(x) x.firstIncreseTimeIntAgainstForce, curTracksNAtalG2(realG2tal));
% %         curTimeDelay = arrayfun(@(x) x.firstIncreseTimeIntAgainstForce, curTracksNAtalG2(curForceTransmittingG2 & curStartingAmpG2<thresAmp));
%         InitTimeLagMatureTal{ii} = curTimeDelayG2;
%         tracksNAtalG2{ii}=curTracksNAtalG2(curForceTransmittingG2 & curStartingAmpG2<thresAmp);
%         clear curTracksNAtalG2
%         disp([talFolder{ii} ' has been analyzed for t_initForMaturingAdhesion.'])
%     catch
%         disp([talFolder{ii} ' has not been analyzed for t_initForMaturingAdhesion.'])
%     end
% end
% totalInitTimeLagMatureTal =[InitTimeLagMatureTal{1}];
% totalInitTimeLagMatureTal =[InitTimeLagMatureTal{1}; InitTimeLagMatureTal{2}*2; 
%     InitTimeLagMatureTal{3}*2; InitTimeLagMatureTal{4}*2; InitTimeLagMatureTal{5}*2; ...
%     InitTimeLagMatureTal{6}*2; InitTimeLagMatureTal{7}*2];
%% filtering with high and low limit (-60 and 60)
lowTLag=-80; hiTLag=60;
totalInitTimeLagMatureTal2 = totalInitTimeLagMatureTal(totalInitTimeLagMatureTal>lowTLag & totalInitTimeLagMatureTal<hiTLag);
totalInitTimeLagMatureVin2 = totalInitTimeLagMatureVin(totalInitTimeLagMatureVin>lowTLag & totalInitTimeLagMatureVin<hiTLag);
totalInitTimeLagMaturePax2 = totalInitTimeLagMaturePax(totalInitTimeLagMaturePax>lowTLag & totalInitTimeLagMaturePax<hiTLag);

disp([cellfun(@(x) nanmedian(x),InitTimeLagMaturePax); cellfun(@(x) length(x),InitTimeLagMaturePax)])
disp([cellfun(@(x) nanmedian(x),InitTimeLagMatureVin); cellfun(@(x) length(x),InitTimeLagMatureVin)])
disp([cellfun(@(x) nanmedian(x),InitTimeLagMatureTal); cellfun(@(x) length(x),InitTimeLagMatureTal)])
%% G2 representation - talin
CurrentFrameTal2=67;
iRepTal2=3;
% data loading
talForceStack = load([talFolder{iRepTal2} filesep 'fMap' filesep 'tMap.mat'],'tMap');
talForceStack = talForceStack.tMap;
talImgStack = load([talFolder{iRepTal2} filesep 'pax' filesep 'paxImgStack.mat'],'paxImgStack');
talImgStack = talImgStack.paxImgStack;
forceFieldTal = load([MDpathTal filesep 'TFMPackage' filesep 'forceField' filesep 'forceField.mat']);
forceFieldTal = forceFieldTal.forceField;
displFieldTal = load([MDpathTal filesep 'TFMPackage' filesep 'correctedDisplacementField' filesep 'displField.mat']);
displFieldTal = displFieldTal.displField;


figure('Position', pos, 'PaperPositionMode', 'auto','Color','w');
axes('Position', [0 5/6 1/6 1/6]); avgWidth=0; % avgWidth
imshow(imcomplement(mean(talImgStack(:,:,CurrentFrameTal2-avgWidth:CurrentFrameTal2+avgWidth),3)),[])

% Use the filtered tracks for G2 above
tracksTalStruct = load([talFolder{iRepTal2} filesep 'data' filesep 'tracksG2real.mat'],'tracksG2');
tracksNAtal2 = tracksTalStruct.tracksG2;
% See the iInit for these tracks
disp(InitTimeLagMatureTal{iRepTal2})
% See the general trend of 8th track
iNATal=9;
showSingleAdhesionTrackSummary(MDtal,tracksNAtal2(iNATal),talImgStack,talForceStack,iNATal);























%% Maturing adhesions - t_init - plotting
[lengthLongest,iLongest]=max([length(totalInitTimeLagMaturePax2),length(totalInitTimeLagMatureVin2),length(totalInitTimeLagMatureTal2)]);
% it was vinculin that was longest
matrixInitTimeLagG2 = NaN(lengthLongest,3);
matrixInitTimeLagG2(1:length(totalInitTimeLagMaturePax2),1) = totalInitTimeLagMaturePax2;
matrixInitTimeLagG2(1:length(totalInitTimeLagMatureVin2),2) = totalInitTimeLagMatureVin2;
matrixInitTimeLagG2(1:length(totalInitTimeLagMatureTal2),3) = totalInitTimeLagMatureTal2;
% axes('Position',[0.08,0.45,0.40,0.18])
axes('Position', [80/600 25/600 140/600 65/600]); % 
% axes('Position', [510/600 325/600 80/600 65/600]); % 
boxplot(matrixInitTimeLagG2,'orientation','horizontal','whisker',0.5,'notch','on',...
    'labels',{['Pax (N=' num2str(length(totalInitTimeLagMaturePax)) ')'],['Vin (N=' num2str(length(totalInitTimeLagMatureVin)) ')'],...
    ['Tal (N=' num2str(length(totalInitTimeLagMatureTal)) ')']},'symbol','','widths',boxWidth,'jitter',1,'colors','k')
txt = findobj(gca,'Type','text');
set(txt(3:end),'VerticalAlignment', 'Middle');
set(findobj(gca,'LineStyle','--'),'LineStyle','-')
set(findobj(gca,'tag','Median'),'LineWidth',2)
xlim([-35 30])
xlabel('Time lag in initial rise (s)')
[pPaxVsVin2,hPaxVsVin2] = ranksum(totalInitTimeLagMaturePax2,totalInitTimeLagMatureVin2);
[pPaxVsTal2,hPaxVsTal2] = ranksum(totalInitTimeLagMaturePax2,totalInitTimeLagMatureTal2);
[pVinVsTal2,hVinVsTal2] = ranksum(totalInitTimeLagMatureVin2,totalInitTimeLagMatureTal2);

text(14,1.5,['p=' num2str(pPaxVsVin2,2)]) %pax vs vin
text(14,2.5,['p=' num2str(pVinVsTal2,2)]) %tal vs vin
text(14,2,['p=' num2str(pPaxVsTal2,2)]) %pax vs tal

% median values
medPaxInit2=median(totalInitTimeLagMaturePax2);
medVinInit2=median(totalInitTimeLagMatureVin2);
medTalInit2=median(totalInitTimeLagMatureTal2);
text(medPaxInit2,1.5,num2str(medPaxInit2,'%10.0f'),'HorizontalAlignment','center') %pax vs vin
text(medVinInit2,2.5,num2str(medVinInit2,'%10.0f'),'HorizontalAlignment','center') %tal vs vin
text(medTalInit2,3.5,num2str(medTalInit2,'%10.0f'),'HorizontalAlignment','center') %pax vs tal

set(findobj(gca,'Type','Text'),'FontSize',7)
set(gca,'FontSize',7)
title('maturing adhesions (G2)')
% totalInitTimeLagMature{1} = totalInitTimeLagMaturePax2;
% totalInitTimeLagMature{2} = totalInitTimeLagMatureVin2;
% totalInitTimeLagMature{3} = totalInitTimeLagMatureTal2;
% figure
% boxPlotCellArray(totalInitTimeLagMature,{'pax','vin','tal'},1,0)
%% Molecular association rate
%% Re-calculate/normalize the rates. 3 time frames were creating almost random distribution
tR = 5; % time of regression
iT_vp = 0;
%% pax recal
for ii=1:numel(paxFolder)
    tracksG2=load([paxFolder{ii} filesep 'data' filesep 'tracksG2real.mat'],'tracksG2');
    curTracksNApaxG2 = tracksG2.tracksG2;
    tracksG1 = load([paxFolder{ii} filesep 'data' filesep 'tracksG1real.mat'],'tracksG1');
    curTracksNApaxG1=tracksG1.tracksG1;
    curMDpath = fileparts(fileparts(paxFolder{ii}));
    curMD = load([curMDpath filesep 'movieData.mat']);
    curMD = curMD.MD;
    tInterval = curMD.timeInterval_;
    
%     [curTracksNApaxG1,timeG1] = calculateFirstIncreaseTimeTracks(curTracksNApaxG1,splineParamInit,preDetecFactor,tInterval);
%     [curTracksNApaxG2,timeG2] = calculateFirstIncreaseTimeTracks(curTracksNApaxG2,splineParamInit,preDetecFactor,tInterval);
%     transG1=~isnan(timeG1);
%     transG2=~isnan(timeG2);
%     longEnoughG1 = arrayfun(@(x) (x.endingFrameExtra-x.firstIncreseTimeInt+1)>iT_vp+tR,curTracksNApaxG1(transG1));
%     longEnoughG2 = arrayfun(@(x) (x.endingFrameExtra-x.firstIncreseTimeInt+1)>iT_vp+tR,curTracksNApaxG2(transG2));
%     indTransG1=find(transG1);
%     indTransG2=find(transG2);
%     transG1(indTransG1(~longEnoughG1))=0;
%     transG2(indTransG2(~longEnoughG2))=0;
%     [~,curEarlyAssmRateG1] = arrayfun(@(x) regression(1:tR,...
%         log(x.ampTotal(round(x.firstIncreseTimeInt/tInterval)+iT_vp:round(x.firstIncreseTimeInt/tInterval)+tR-1+iT_vp)...
%         /x.ampTotal(round(x.firstIncreseTimeInt/tInterval)+iT_vp))),curTracksNApaxG1(transG1));
%     [~,curEarlyAssmRateG2] = arrayfun(@(x) regression(1:tR,...
%         log(x.ampTotal(round(x.firstIncreseTimeInt/tInterval)+iT_vp:round(x.firstIncreseTimeInt/tInterval)+tR-1+iT_vp)...
%         /x.ampTotal(round(x.firstIncreseTimeInt/tInterval)+iT_vp))),curTracksNApaxG2(transG2));
    % Get the assembly rate and force growth rate
    longEnoughG1 = arrayfun(@(x) (x.endingFrameExtra-x.startingFrameExtra+1)>iT_vp+tR,curTracksNApaxG1);
    longEnoughG2 = arrayfun(@(x) (x.endingFrameExtra-x.startingFrameExtra+1)>iT_vp+tR,curTracksNApaxG2);
    % Get the assembly rate and force growth rate
    [~,curEarlyAssmRateG1] = arrayfun(@(x) regression(1:tR,...
        log(x.ampTotal(x.startingFrameExtra+iT_vp:x.startingFrameExtra+tR-1+iT_vp)/x.ampTotal(x.startingFrameExtra+iT_vp))),curTracksNApaxG1(longEnoughG1));
    [~,curEarlyAssmRateG2] = arrayfun(@(x) regression(1:tR,...
        log(x.ampTotal(x.startingFrameExtra+iT_vp:x.startingFrameExtra+tR-1+iT_vp)/x.ampTotal(x.startingFrameExtra+iT_vp))),curTracksNApaxG2(longEnoughG2));

    earlyAssmRateG1All{ii} = curEarlyAssmRateG1; %cur_tlagInitPax;
    earlyAssmRateG2All{ii} = curEarlyAssmRateG2; %cur_tlagInitPax;
    
    [~,earlyForceRateG1] = arrayfun(@(x) regression(1:tR,x.forceMag(x.startingFrameExtra+iT_vp:x.startingFrameExtra+tR-1+iT_vp)),curTracksNApaxG1(longEnoughG1));
    [~,earlyForceRateG2] = arrayfun(@(x) regression(1:tR,x.forceMag(x.startingFrameExtra+iT_vp:x.startingFrameExtra+tR-1+iT_vp)),curTracksNApaxG2(longEnoughG2));
    earlyForceRateG1All{ii} = earlyForceRateG1; %cur_tlagInitPax;
    earlyForceRateG2All{ii} = earlyForceRateG2; %cur_tlagInitPax;
    disp([paxFolder{ii} ' done.'])
end
totalEarlyAssmRateG1Pax=[earlyAssmRateG1All{1}; earlyAssmRateG1All{2}; earlyAssmRateG1All{3}; earlyAssmRateG1All{4}];
totalEarlyAssmRateG2Pax=[earlyAssmRateG2All{1}; earlyAssmRateG2All{2}; earlyAssmRateG2All{3}; earlyAssmRateG2All{4}];
totalEarlyForceRateG1Pax=[earlyForceRateG1All{1}; earlyForceRateG1All{2}; earlyForceRateG1All{3}; earlyForceRateG1All{4}];
totalEarlyForceRateG2Pax=[earlyForceRateG2All{1}; earlyForceRateG2All{2}; earlyForceRateG2All{3}; earlyForceRateG2All{4}];

%% vinculin
% tR = 10; % time of regression
% iT_vp = iT;
tRf = tR; % time of regression
for ii=1:numel(vinFolder)
%     tracksG2=load([vinFolder{ii} filesep 'data' filesep 'tracksG2.mat'],'tracksG2');
    tracksG2=load([vinFolder{ii} filesep 'data' filesep 'tracksG2real.mat'],'tracksG2');
    curTracksNAvinG2 = tracksG2.tracksG2;
%     tracksG1 = load([vinFolder{ii} filesep 'data' filesep 'tracksG1.mat'],'tracksG1');
    tracksG1 = load([vinFolder{ii} filesep 'data' filesep 'tracksG1real.mat'],'tracksG1');
    curTracksNAvinG1=tracksG1.tracksG1;
    curMDpath = fileparts(fileparts(vinFolder{ii}));
    curMD = load([curMDpath filesep 'movieData.mat']);
    curMD = curMD.MD;
    tInterval = curMD.timeInterval_;

    longEnoughG2 = arrayfun(@(x) (x.endingFrameExtra-x.startingFrameExtra+1)>iT_vp+tR,curTracksNAvinG2);
    longEnoughG1 = arrayfun(@(x) (x.endingFrameExtra-x.startingFrameExtra+1)>iT_vp+tR,curTracksNAvinG1);
    % Get the assembly rate and force growth rate
    [~,curEarlyAssmRateG1] = arrayfun(@(x) regression(1:tR,...
        log(x.ampTotal(x.startingFrameExtra+iT_vp:x.startingFrameExtra+tR-1+iT_vp)/x.ampTotal(x.startingFrameExtra+iT_vp))),curTracksNAvinG1(longEnoughG1));
    [~,curEarlyAssmRateG2] = arrayfun(@(x) regression(1:tR,...
        log(x.ampTotal(x.startingFrameExtra+iT_vp:x.startingFrameExtra+tR-1+iT_vp)/x.ampTotal(x.startingFrameExtra+iT_vp))),curTracksNAvinG2(longEnoughG2));
%     [curTracksNAvinG1,timeG1] = calculateFirstIncreaseTimeTracks(curTracksNAvinG1,splineParamInit,preDetecFactor,tInterval);
%     [curTracksNAvinG2,timeG2] = calculateFirstIncreaseTimeTracks(curTracksNAvinG2,splineParamInit,preDetecFactor,tInterval);
%     transG1=~isnan(timeG1);
%     transG2=~isnan(timeG2);
%     longEnoughG1 = arrayfun(@(x) (x.endingFrameExtra-x.firstIncreseTimeInt+1)>iT_vp+tR,curTracksNAvinG1(transG1));
%     longEnoughG2 = arrayfun(@(x) (x.endingFrameExtra-x.firstIncreseTimeInt+1)>iT_vp+tR,curTracksNAvinG2(transG2));
%     indTransG1=find(transG1);
%     indTransG2=find(transG2);
%     transG1(indTransG1(~longEnoughG1))=0;
%     transG2(indTransG2(~longEnoughG2))=0;
%     % Get the assembly rate and force growth rate
%     [~,curEarlyAssmRateG1] = arrayfun(@(x) regression(1:tR,...
%         log(x.ampTotal(round(x.firstIncreseTimeInt/tInterval)+iT_vp:round(x.firstIncreseTimeInt/tInterval)+tR-1+iT_vp)...
%         /x.ampTotal(round(x.firstIncreseTimeInt/tInterval)+iT_vp))),curTracksNAvinG1(transG1));
%     [~,curEarlyAssmRateG2] = arrayfun(@(x) regression(1:tR,...
%         log(x.ampTotal(round(x.firstIncreseTimeInt/tInterval)+iT_vp:round(x.firstIncreseTimeInt/tInterval)+tR-1+iT_vp)...
%         /x.ampTotal(round(x.firstIncreseTimeInt/tInterval)+iT_vp))),curTracksNAvinG2(transG2));

    earlyAssmRateG1All{ii} = curEarlyAssmRateG1; %cur_tlagInitPax;
    earlyAssmRateG2All{ii} = curEarlyAssmRateG2; %cur_tlagInitPax;

    longEnoughG2 = arrayfun(@(x) (x.endingFrameExtra-x.startingFrameExtra+1)>iT_vp+tRf,curTracksNAvinG2);
    longEnoughG1 = arrayfun(@(x) (x.endingFrameExtra-x.startingFrameExtra+1)>iT_vp+tRf,curTracksNAvinG1);
    [~,earlyForceRateG1] = arrayfun(@(x) regression(1:tRf,x.forceMag(x.startingFrameExtra+iT_vp:x.startingFrameExtra+tRf-1+iT_vp)),curTracksNAvinG1(longEnoughG1));
    [~,earlyForceRateG2] = arrayfun(@(x) regression(1:tRf,x.forceMag(x.startingFrameExtra+iT_vp:x.startingFrameExtra+tRf-1+iT_vp)),curTracksNAvinG2(longEnoughG2));
    earlyForceRateG1All{ii} = earlyForceRateG1; %cur_tlagInitPax;
    earlyForceRateG2All{ii} = earlyForceRateG2; %cur_tlagInitPax;
    disp([vinFolder{ii} ' done.'])
end
totalEarlyAssmRateG1Vin=[earlyAssmRateG1All{1}; earlyAssmRateG1All{2}; earlyAssmRateG1All{3}; earlyAssmRateG1All{4}; earlyAssmRateG1All{5}; earlyAssmRateG1All{6}];
totalEarlyAssmRateG2Vin=[earlyAssmRateG2All{1}; earlyAssmRateG2All{2}; earlyAssmRateG2All{3}; earlyAssmRateG2All{4}; earlyAssmRateG2All{5}; earlyAssmRateG2All{6}];
totalEarlyForceRateG1Vin=[earlyForceRateG1All{1}; earlyForceRateG1All{2}; earlyForceRateG1All{3}; earlyForceRateG1All{4}; earlyForceRateG1All{5}; earlyForceRateG1All{6}];
totalEarlyForceRateG2Vin=[earlyForceRateG2All{1}; earlyForceRateG2All{2}; earlyForceRateG2All{3}; earlyForceRateG2All{4}; earlyForceRateG2All{5}; earlyForceRateG2All{6}];
%% talin- recal
% tR = 10; % time of regression
iT_t = iT_vp;
for ii=1:1
    tracksG2=load([talFolder{ii} filesep 'data' filesep 'tracksG2real.mat'],'tracksG2');
    curTracksNAtalG2 = tracksG2.tracksG2.tracksG2;
    tracksG1 = load([talFolder{ii} filesep 'data' filesep 'tracksG1real.mat'],'tracksG1');
    curTracksNAtalG1=tracksG1.tracksG1.tracksG1;
    longEnoughG2 = arrayfun(@(x) (x.endingFrameExtra-x.startingFrameExtra+1)>iT_t+tR,curTracksNAtalG2);
    longEnoughG1 = arrayfun(@(x) (x.endingFrameExtra-x.startingFrameExtra+1)>iT_t+tR,curTracksNAtalG1);
    
    % Get the assembly rate and force growth rate
    [~,curEarlyAssmRateG1] = arrayfun(@(x) regression(1:tR,...
        log(x.ampTotal(x.startingFrameExtra+iT_t:x.startingFrameExtra+tR-1+iT_t)/x.ampTotal(x.startingFrameExtra+iT_t))),curTracksNAtalG1(longEnoughG1));
    [~,curEarlyAssmRateG2] = arrayfun(@(x) regression(1:tR,...
        log(x.ampTotal(x.startingFrameExtra+iT_t:x.startingFrameExtra+tR-1+iT_t)/x.ampTotal(x.startingFrameExtra+iT_t))),curTracksNAtalG2(longEnoughG2));
%     [curTracksNAtalG1,timeG1] = calculateFirstIncreaseTimeTracks(curTracksNAtalG1,splineParamInit,preDetecFactor,tInterval);
%     [curTracksNAtalG2,timeG2] = calculateFirstIncreaseTimeTracks(curTracksNAtalG2,splineParamInit,preDetecFactor,tInterval);
%     transG1=~isnan(timeG1);
%     transG2=~isnan(timeG2);
%     longEnoughG1 = arrayfun(@(x) (x.endingFrameExtra-x.firstIncreseTimeInt+1)>iT_vp+tR,curTracksNAtalG1(transG1));
%     longEnoughG2 = arrayfun(@(x) (x.endingFrameExtra-x.firstIncreseTimeInt+1)>iT_vp+tR,curTracksNAtalG2(transG2));
%     indTransG1=find(transG1);
%     indTransG2=find(transG2);
%     transG1(indTransG1(~longEnoughG1))=0;
%     transG2(indTransG2(~longEnoughG2))=0;
    
    % Get the assembly rate and force growth rate
%     [~,curEarlyAssmRateG1] = arrayfun(@(x) regression(1:tR,...
%         log(x.ampTotal(round(x.firstIncreseTimeInt/tInterval)+iT_vp:round(x.firstIncreseTimeInt/tInterval)+tR-1+iT_vp)...
%         /x.ampTotal(round(x.firstIncreseTimeInt/tInterval)+iT_vp))),curTracksNAtalG1(longEnoughG1));
%     [~,curEarlyAssmRateG2] = arrayfun(@(x) regression(1:tR,...
%         log(x.ampTotal(round(x.firstIncreseTimeInt/tInterval)+iT_vp:round(x.firstIncreseTimeInt/tInterval)+tR-1+iT_vp)...
%         /x.ampTotal(round(x.firstIncreseTimeInt/tInterval)+iT_vp))),curTracksNAtalG2(longEnoughG2));
    
    earlyAssmRateG1All{ii} = curEarlyAssmRateG1; %cur_tlagInitPax;
    earlyAssmRateG2All{ii} = curEarlyAssmRateG2; %cur_tlagInitPax;

    [~,earlyForceRateG1] = arrayfun(@(x) regression(1:tR,x.forceMag(x.startingFrameExtra+iT_t:x.startingFrameExtra+tR-1+iT_t)),curTracksNAtalG1(longEnoughG1));
    [~,earlyForceRateG2] = arrayfun(@(x) regression(1:tR,x.forceMag(x.startingFrameExtra+iT_t:x.startingFrameExtra+tR-1+iT_t)),curTracksNAtalG2(longEnoughG2));
    earlyForceRateG1All{ii} = earlyForceRateG1; %cur_tlagInitPax;
    earlyForceRateG2All{ii} = earlyForceRateG2; %cur_tlagInitPax;
    disp([vinFolder{ii} ' done.'])
end
totalEarlyAssmRateG1Tal=[earlyAssmRateG1All{1}];
totalEarlyAssmRateG2Tal=[earlyAssmRateG2All{1}];
totalEarlyForceRateG1Tal=[earlyForceRateG1All{1}];
totalEarlyForceRateG2Tal=[earlyForceRateG2All{1}];
%% Molecular association rate
earlyAmpSlopeG1G2={totalEarlyAssmRateG1Tal,totalEarlyAssmRateG2Tal,...
    totalEarlyAssmRateG1Vin,totalEarlyAssmRateG2Vin,totalEarlyAssmRateG1Pax,totalEarlyAssmRateG2Pax};
nameListG1G2={'Tal_G1','Tal_G2','Vin_G1','Vin_G2','Pax_G1','Pax_G2'};
% axes('Position', [380/600 225/600 80/600 65/600]); % 
axes('Position', [350/600 (115-50)/600 50/600 (80+50)/600]);
% barPlotCellArray(earlyAmpSlopeG1G2,nameListG1G2,1)
boxPlotCellArray(earlyAmpSlopeG1G2,nameListG1G2,1,0,0,2)

% do some statistics:
[pTalG1G2,hTalG1G2] = ranksum(totalEarlyAssmRateG1Tal,totalEarlyAssmRateG2Tal);
[pVinG1G2,hAmpVinG1G2] = ranksum(totalEarlyAssmRateG1Vin,totalEarlyAssmRateG2Vin);
[pPaxG1G2,hPaxG1G2] = ranksum(totalEarlyAssmRateG1Pax,totalEarlyAssmRateG2Pax);
%put it in the plot
text(1,0.08,['p=' num2str(pTalG1G2,2)]) %pax vs vin
text(3,0.09,['p=' num2str(pVinG1G2,2)]) %tal vs vin
text(5,0.12,['p=' num2str(pPaxG1G2,2)]) %pax vs tal

ylabel('Assembly rate (a.u./sec)')
title({'Assembly rate'; ['for first ' num2str(tR) ' sec']})
set(findobj(gca,'Type','Text'),'FontSize',7)
set(gca,'FontSize',7)
xlim([0.4 6.6])
ylim([-0.07 0.13])
%% dF/dt on initial molecular binding

earlyForceSlopeG1G2={totalEarlyForceRateG1Tal, totalEarlyForceRateG2Tal,...
    totalEarlyForceRateG1Vin,totalEarlyForceRateG2Vin,totalEarlyForceRateG1Pax,totalEarlyForceRateG2Pax};
% axes('Position', [490/600 (115-50)/600 50/600 (80+50)/600]); % 
figure('Position',[100 100 150 200])
boxPlotCellArray(earlyForceSlopeG1G2,nameListG1G2,1,0,0,2)
% barPlotCellArray(earlyForceSlopeG1G2,nameListG1G2,1)
ylabel('Force growth rate (Pa/sec)')
title({'Force growth rate'; ['for ' num2str(2*tR) ' sec after ' num2str(2*iT_vp) ' sec']})
% do some statistics:
[pForceTalG1G2,hForceTalG1G2] = ranksum(totalEarlyForceRateG1Tal,totalEarlyForceRateG2Tal);
[pForceVinG1G2,hForceVinG1G2] = ranksum(totalEarlyForceRateG1Vin,totalEarlyForceRateG2Vin);
[pForcePaxG1G2,hForcePaxG1G2] = ranksum(totalEarlyForceRateG1Pax,totalEarlyForceRateG2Pax);
%put it in the plot
text(1,3,['p=' num2str(pForceTalG1G2,2)]) %pax vs vin
text(3,3.5,['p=' num2str(pForceVinG1G2,2)]) %pax vs vin
text(5,4,['p=' num2str(pForcePaxG1G2,2)]) %pax vs vin

set(findobj(gca,'Type','Text'),'FontSize',7)
set(gca,'FontSize',7)
xlim([0.4 6.6])
% ylim([-12 12])
curName = ['ForceGrowth' num2str(2*tR) 'sec after ' num2str(2*iT_vp) ' sec'];
print('-depsc', '-loose', [f3_Tmeasure_path filesep 'FigS3' curName '.eps'])
print('-dtiff', '-loose', '-r300', [f3_Tmeasure_path filesep 'FigS3' curName '.tif'])
savefig([f3_Tmeasure_path filesep 'FigS3' curName '.fig'])

%% Replotting time series - talin G1
tracksG1 = load([talFolder{1} filesep 'data' filesep 'tracksG1real.mat'],'tracksG1');
tracksNAtalG1=tracksG1.tracksG1.tracksG1;
tracksG2 = load([talFolder{1} filesep 'data' filesep 'tracksG2real.mat'],'tracksG2');
tracksNAtalG2=tracksG2.tracksG2.tracksG2;
axes('Position', [130/600 530/600 65/600 65/600]); % 
plotIntensityForce(tracksNAtalG1,[],false,false,'UseCurrentAxis',true,...
    'Source',{'ampTotal'},'plotCohorts',true,'tInterval',MDtal.timeInterval_,'prePostFrames',10,...
    'yNormalization',true,'numCohorts',1)
xlim([-10 100])
ylim([210 322])
axes('Position', [230/600 530/600 65/600 65/600]); % 
plotIntensityForce(tracksNAtalG2,[],false,false,'UseCurrentAxis',true,...
    'Source',{'ampTotal'},'plotCohorts',true,'tInterval',MDtal.timeInterval_,'prePostFrames',10,'yNormalization',true,'numCohorts',1)
xlim([-10 100])
ylim([210 322])

axes('Position', [130/600 430/600 65/600 65/600]); % only force-transmitting
plotIntensityForce(tracksNAtalG1,[],false,false,'UseCurrentAxis',true,...
    'Source',{'forceMag'},'plotCohorts',true,'tInterval',MDtal.timeInterval_,'prePostFrames',10,'yNormalization',true,'numCohorts',1)
xlim([-10 100])
ylim([15 120])
axes('Position', [230/600 430/600 65/600 65/600]); % only force-transmitting
plotIntensityForce(tracksNAtalG2,[],false,false,'UseCurrentAxis',true,...
    'Source',{'forceMag'},'plotCohorts',true,'tInterval',MDtal.timeInterval_,'prePostFrames',10,'yNormalization',true,'numCohorts',1)
xlim([-10 100])
ylim([15 120])
%% Replotting time series - vinculin G1
tracksG1 = load([vinFolder{4} filesep 'data' filesep 'tracksG1real.mat'],'tracksG1');
tracksNAvinG1=tracksG1.tracksG1;
tracksG2 = load([vinFolder{4} filesep 'data' filesep 'tracksG2real.mat'],'tracksG2');
tracksNAvinG2=tracksG2.tracksG2;
axes('Position', [130/600 330/600 65/600 65/600]); % 
plotIntensityForce(tracksNAvinG1,[],false,false,'UseCurrentAxis',true,...
    'Source',{'ampTotal'},'plotCohorts',true,'tInterval',MDvin.timeInterval_,'prePostFrames',10,...
    'yNormalization',true,'numCohorts',1)
xlim([-10 120])
ylim([8000 13000])
axes('Position', [230/600 330/600 65/600 65/600]); % 
plotIntensityForce(tracksNAvinG2,[],false,false,'UseCurrentAxis',true,...
    'Source',{'ampTotal'},'plotCohorts',true,'tInterval',MDvin.timeInterval_,'prePostFrames',10,'yNormalization',true,'numCohorts',1)
xlim([-10 120])
ylim([8000 13000])

axes('Position', [130/600 230/600 65/600 65/600]); % only force-transmitting
plotIntensityForce(tracksNAvinG1,[],false,false,'UseCurrentAxis',true,...
    'Source',{'forceMag'},'plotCohorts',true,'tInterval',MDvin.timeInterval_,'prePostFrames',10,'yNormalization',true,'numCohorts',1)
xlim([-10 120])
ylim([15 1000])
axes('Position', [230/600 230/600 65/600 65/600]); % only force-transmitting
plotIntensityForce(tracksNAvinG2,[],false,false,'UseCurrentAxis',true,...
    'Source',{'forceMag'},'plotCohorts',true,'tInterval',MDvin.timeInterval_,'prePostFrames',10,'yNormalization',true,'numCohorts',1)
xlim([-10 120])
ylim([15 1000])
%% Replotting time series - paxillin G1
tracksG1 = load([paxFolder{3} filesep 'data' filesep 'tracksG1real.mat'],'tracksG1');
tracksNApaxG1=tracksG1.tracksG1;
tracksG2 = load([paxFolder{3} filesep 'data' filesep 'tracksG2real.mat'],'tracksG2');
tracksNApaxG2=tracksG2.tracksG2;
axes('Position', [130/600 130/600 65/600 65/600]); % 
plotIntensityForce(tracksNApaxG1,[],false,false,'UseCurrentAxis',true,...
    'Source',{'ampTotal'},'plotCohorts',true,'tInterval',MDvin.timeInterval_,'prePostFrames',10,...
    'yNormalization',true,'numCohorts',1)
xlim([-10 120])
ylim([240 580])
axes('Position', [230/600 130/600 65/600 65/600]); % 
plotIntensityForce(tracksNApaxG2,[],false,false,'UseCurrentAxis',true,...
    'Source',{'ampTotal'},'plotCohorts',true,'tInterval',MDvin.timeInterval_,'prePostFrames',10,'yNormalization',true,'numCohorts',1)
xlim([-10 120])
ylim([240 580])

axes('Position', [130/600 30/600 65/600 65/600]); % only force-transmitting
plotIntensityForce(tracksNApaxG1,[],false,false,'UseCurrentAxis',true,...
    'Source',{'forceMag'},'plotCohorts',true,'tInterval',MDvin.timeInterval_,'prePostFrames',10,'yNormalization',true,'numCohorts',1)
xlim([-10 120])
ylim([0 35])
axes('Position', [230/600 30/600 65/600 65/600]); % only force-transmitting
plotIntensityForce(tracksNApaxG2,[],false,false,'UseCurrentAxis',true,...
    'Source',{'forceMag'},'plotCohorts',true,'tInterval',MDvin.timeInterval_,'prePostFrames',10,'yNormalization',true,'numCohorts',1)
xlim([-10 120])
ylim([0 35])
%% save
f3_Tmeasure_path = '/home2/shan/ownCloud/documents/NascentAdhesionForce/Manuscript/Figures/Fig3timeMeasurement';
if ~exist(f3_Tmeasure_path,'dir')
    mkdir(f3_Tmeasure_path)
end
print('-depsc', '-loose', [f3_Tmeasure_path filesep 'Fig3TimeLags.eps'])
print('-dtiff', '-loose', '-r300', [f3_Tmeasure_path filesep 'Fig3TimeLags.tif'])
savefig([f3_Tmeasure_path filesep 'Fig3TimeLags.fig'])
%% plot for individual ones
pos(3:4) = [400 200];
figure('Position', pos, 'PaperPositionMode', 'auto','Color','w');
axes('Position', [30/400 130/200 65/400 65/200]); 
plotIntensityForce(tracksNApaxG1,[],false,false,'UseCurrentAxis',true,'Source',{'ampTotal'},...
    'plotCohorts',true,'tInterval',MDpax.timeInterval_,'prePostFrames',10,'onlyFirstMode',true,'yNormalization',true,'numCohorts',1)
set(gca,'XLim',[-10 100],'YLim',[100 800])

axes('Position', [30/400 30/200 65/400 65/200]); 
plotIntensityForce(tracksNApaxG1,[],false,false,'UseCurrentAxis',true,'Source',{'forceMag'},...
    'plotCohorts',true,'tInterval',MDpax.timeInterval_,'prePostFrames',10,'onlyFirstMode',true,'yNormalization',true,'numCohorts',1)
set(gca,'XLim',[-10 100],'YLim',[-3 140])

axes('Position', [130/400 130/200 65/400 65/200]); 
plotIntensityForce(tracksNApaxG1,[],false,false,'UseCurrentAxis',true,'Source',{'ampTotal'},...
    'plotCohorts',true,'tInterval',MDpax.timeInterval_,'prePostFrames',10,'onlyFirstMode',true,'yNormalization',true,'numCohorts',1)
set(gca,'XLim',[-10 20],'YLim',[100 800])

axes('Position', [130/400 30/200 65/400 65/200]);
plotIntensityForce(tracksNApaxG1,[],false,false,'UseCurrentAxis',true,'Source',{'forceMag'},...
    'plotCohorts',true,'tInterval',MDpax.timeInterval_,'prePostFrames',10,'onlyFirstMode',true,'yNormalization',true,'numCohorts',1)
set(gca,'XLim',[-10 20],'YLim',[-2 20])
%% G2
axes('Position', [230/400 130/200 65/400 65/200]); 
plotIntensityForce(tracksNApaxG2,[],false,false,'UseCurrentAxis',true,'Source',{'ampTotal'},...
    'plotCohorts',true,'tInterval',MDpax.timeInterval_,'prePostFrames',10,'onlyFirstMode',true,'yNormalization',true,'numCohorts',1)
set(gca,'XLim',[-10 100],'YLim',[100 800])

axes('Position', [230/400 30/200 65/400 65/200]);
plotIntensityForce(tracksNApaxG2,[],false,false,'UseCurrentAxis',true,'Source',{'forceMag'},...
    'plotCohorts',true,'tInterval',MDpax.timeInterval_,'prePostFrames',10,'onlyFirstMode',true,'yNormalization',true,'numCohorts',1)
set(gca,'XLim',[-10 100],'YLim',[-3 140])

axes('Position', [330/400 130/200 65/400 65/200]); 
plotIntensityForce(tracksNApaxG2,[],false,false,'UseCurrentAxis',true,'Source',{'ampTotal'},...
    'plotCohorts',true,'tInterval',MDpax.timeInterval_,'prePostFrames',10,'onlyFirstMode',true,'yNormalization',true,'numCohorts',1)
set(gca,'XLim',[-10 20],'YLim',[100 800])

axes('Position', [330/400 30/200 65/400 65/200]); 
plotIntensityForce(tracksNApaxG2,[],false,false,'UseCurrentAxis',true,'Source',{'forceMag'},...
    'plotCohorts',true,'tInterval',MDpax.timeInterval_,'prePostFrames',10,'onlyFirstMode',true,'yNormalization',true,'numCohorts',1)
set(gca,'XLim',[-10 20],'YLim',[-2 20])
% Message: There are heterogeneity, but overage response for initial 10 sec
% shows the force increases are the same for both cases.

%% save
f3_Tmeasure_path = '/home2/shan/ownCloud/documents/NascentAdhesionForce/Manuscript/Figures/FigS3individualPlots';
if ~exist(f3_Tmeasure_path,'dir')
    mkdir(f3_Tmeasure_path)
end
print('-depsc', '-loose', [f3_Tmeasure_path filesep 'FigS3individualPlots.eps'])
print('-dtiff', '-loose', '-r300', [f3_Tmeasure_path filesep 'FigS3individualPlots.tif'])
savefig([f3_Tmeasure_path filesep 'FigS3individualPlots.fig'])
%% SI - Peak lag
[lengthLongest,iLongest]=max([length(totalPeakTLagPax),length(totalPeakTLagVin),length(totalPeakTLagTal)]);
% it was vinculin that was longest
matrixPeakPaxVinTal = NaN(lengthLongest,3);
matrixPeakPaxVinTal(1:length(totalPeakTLagPax),1) = totalPeakTLagPax;
matrixPeakPaxVinTal(1:length(totalPeakTLagVin),2) = totalPeakTLagVin;
matrixPeakPaxVinTal(1:length(totalPeakTLagTal),3) = totalPeakTLagTal;

axes('Position',[280/600,200/600,100/600,85/600])

boxplot(matrixPeakPaxVinTal,'orientation','horizontal','whisker',0.5,'notch','on',...
    'labels',{['Pax (N=' num2str(length(totalPeakTLagPax)) ')'],['Vin (N=' num2str(length(totalPeakTLagVin)) ')'],['Tal (N=' num2str(length(totalPeakTLagTal)) ')']},...
    'symbol','','widths',boxWidth,'jitter',1,'colors','k')
set(findobj(gca,'LineStyle','--'),'LineStyle','-')
set(findobj(gca,'tag','Median'),'LineWidth',2)
xlim([-35 30])
xlabel('Time lag in peak (s)')
text(21,1.5,['p=' num2str(pPaxVsVinPeak2,2)]) %pax vs vin
text(21,2.5,['p=' num2str(pVinVsTalPeak2,2)]) %tal vs vin
text(21,2,['p=' num2str(pPaxVsTalPeak2,2)]) %pax vs tal

text(medPaxPeak,1.5,num2str(medPaxPeak,'%10.0f'),'HorizontalAlignment','center') %pax vs vin
text(medVinPeak,2.5,num2str(medVinPeak,'%10.0f'),'HorizontalAlignment','center') %tal vs vin
text(medTalPeak,3.5,num2str(medTalPeak,'%10.0f'),'HorizontalAlignment','center') %pax vs tal
set(findobj(gca,'Type','Text'),'FontSize',7)
set(gca,'FontSize',7)
%% Time to peak - loading
for ii=[1 2 3 4]
    cur_timeToPeakPax = load([paxFolder{ii} '/data/timeToPeaks.mat']);
    cur_timeToPeakPaxInt = cur_timeToPeakPax.timeToPeakForceFT;
    cur_timeToPeakPaxForce = cur_timeToPeakPax.timeToPeakIntFT;
    timeToPeakPaxInt{ii} = cur_timeToPeakPaxInt;
    timeToPeakPaxForce{ii} = cur_timeToPeakPaxForce;
end
totalTimeToPeakPaxInt=[timeToPeakPaxInt{1}; timeToPeakPaxInt{2}; timeToPeakPaxInt{3}; timeToPeakPaxInt{4}];
totalTimeToPeakPaxForce=[timeToPeakPaxForce{1}; timeToPeakPaxForce{2}; timeToPeakPaxForce{3}; timeToPeakPaxForce{4}];

for ii=1:6
    cur_timeToPeakVin = load([vinFolder{ii} '/data/timeToPeaks.mat']);
    cur_timeToPeakVinInt = cur_timeToPeakVin.timeToPeakForceFT;
    cur_timeToPeakVinForce = cur_timeToPeakVin.timeToPeakIntFT;
    timeToPeakVinInt{ii} = cur_timeToPeakVinInt;
    timeToPeakVinForce{ii} = cur_timeToPeakVinForce;
end
totalTimeToPeakVinInt=[timeToPeakVinInt{1}; timeToPeakVinInt{2}; timeToPeakVinInt{3}; 
    timeToPeakVinInt{4}; timeToPeakVinInt{5}; timeToPeakVinInt{6}];
totalTimeToPeakVinForce=[timeToPeakVinForce{1}; timeToPeakVinForce{2}; timeToPeakVinForce{3}; 
    timeToPeakVinForce{4}; timeToPeakVinForce{5}; timeToPeakVinForce{6}];

for ii=1:numel(talFolder)
    cur_timeToPeakTal = load([talFolder{ii} '/data/timeToPeaks.mat']);
    cur_timeToPeakTalInt = cur_timeToPeakTal.timeToPeakForceFT;
    cur_timeToPeakTalForce = cur_timeToPeakTal.timeToPeakIntFT;
    timeToPeakTalInt{ii} = cur_timeToPeakTalInt;
    timeToPeakTalForce{ii} = cur_timeToPeakTalForce;
end
totalTimeToPeakTalInt=[timeToPeakTalInt{1}; timeToPeakTalInt{2}; timeToPeakTalInt{3}; timeToPeakTalInt{4}; timeToPeakTalInt{5}; timeToPeakTalInt{6}];
totalTimeToPeakTalForce=[timeToPeakTalForce{1}; timeToPeakTalForce{2}; timeToPeakTalForce{3}; timeToPeakTalForce{4}; timeToPeakTalForce{5}; timeToPeakTalForce{6}];
%% Stat test for Time to peak
[pPaxVsVinTTPInt,hPaxVsVinTTPInt] = ranksum(totalTimeToPeakPaxInt,totalTimeToPeakVinInt);
[pPaxVsTalTTPInt,hPaxVsTalTTPInt] = ranksum(totalTimeToPeakPaxInt,totalTimeToPeakTalInt);
[pVinVsTalTTPInt,hVinVsTalTTPInt] = ranksum(totalTimeToPeakVinInt,totalTimeToPeakTalInt);

[pPaxVsVinTTPForce,hPaxVsVinTTPForce] = ranksum(totalTimeToPeakPaxForce,totalTimeToPeakVinForce);
[pPaxVsTalTTPForce,hPaxVsTalTTPForce] = ranksum(totalTimeToPeakPaxForce,totalTimeToPeakTalForce);
[pVinVsTalTTPForce,hVinVsTalTTPForce] = ranksum(totalTimeToPeakVinForce,totalTimeToPeakTalForce);
%% Time to peak
[lengthLongest,iLongest]=max([length(totalTimeToPeakPaxInt),length(totalTimeToPeakVinInt),length(totalTimeToPeakTalInt)]);
% it was vinculin that was longest
matrixPaxTalVinTTPeak = NaN(lengthLongest,3);
matrixPaxTalVinTTPeak(1:length(totalTimeToPeakPaxInt),1) = totalTimeToPeakPaxInt;
matrixPaxTalVinTTPeak(1:length(totalTimeToPeakVinInt),2) = totalTimeToPeakVinInt;
matrixPaxTalVinTTPeak(1:length(totalTimeToPeakTalInt),3) = totalTimeToPeakTalInt;

% axes('Position',[0.82,0.44,0.16,0.10])
axes('Position',[280/600,110/600,100/600,60/600])
boxplot(matrixPaxTalVinTTPeak,'orientation','horizontal','whisker',0.5,'notch','on',...
    'labels',{['Pax (N=' num2str(length(totalTimeToPeakPaxInt)) ')'],['Vin (N=' num2str(length(totalTimeToPeakVinInt)) ')'],...
    ['Tal (N=' num2str(length(totalTimeToPeakTalInt)) ')']},'symbol','','widths',boxWidth,'jitter',1,'colors','k')
xlim([-5 120])
text(95,2.5,['p=' num2str(pVinVsTalTTPInt,2)]) %pax vs vin
text(95,1.5,['p=' num2str(pPaxVsVinTTPInt,2)]) %tal vs vin
text(95,2,['p=' num2str(pPaxVsTalTTPInt,2)]) %pax vs tal
txt = findobj(gca,'Type','text');
set(txt(3:end),'VerticalAlignment', 'Middle');
set(findobj(gca,'LineStyle','--'),'LineStyle','-')
set(findobj(gca,'tag','Median'),'LineWidth',2)
xlabel('Time to peak in F.I. (s)')
set(gca,'FontSize',7)
set(txt,'FontSize',7)
% % median values
% text(medPaxInit,1.5,num2str(medPaxInit,'%10.0f'),'HorizontalAlignment','center') %pax vs vin
% text(medVinInit,2.5,num2str(medVinInit,'%10.0f'),'HorizontalAlignment','center') %tal vs vin
% text(medTalInit,3.5,num2str(medTalInit,'%10.0f'),'HorizontalAlignment','center') %pax vs tal
% 
% set(findobj(gca,'Type','Text'),'FontSize',7)

[lengthLongest,iLongest]=max([length(totalTimeToPeakPaxForce),length(totalTimeToPeakVinForce),length(totalTimeToPeakTalForce)]);
% it was vinculin that was longest
matrixPaxTalVinTTPeakForce = NaN(lengthLongest,3);
matrixPaxTalVinTTPeakForce(1:length(totalTimeToPeakPaxForce),1) = totalTimeToPeakPaxForce;
matrixPaxTalVinTTPeakForce(1:length(totalTimeToPeakVinForce),2) = totalTimeToPeakVinForce;
matrixPaxTalVinTTPeakForce(1:length(totalTimeToPeakTalForce),3) = totalTimeToPeakTalForce;

% axes('Position',[0.82,0.28,0.16,0.10])
axes('Position',[280/600,20/600,100/600,60/600])
boxplot(matrixPaxTalVinTTPeakForce,'orientation','horizontal','whisker',0.5,'notch','on',...
    'labels',{['Pax (N=' num2str(length(totalTimeToPeakPaxForce)) ')'],['Vin (N=' num2str(length(totalTimeToPeakVinForce)) ')'],...
    ['Tal (N=' num2str(length(totalTimeToPeakTalForce)) ')']},'symbol','','widths',boxWidth,'jitter',1,'colors','k')
xlim([-5 120])
text(95,2.5,['p=' num2str(pVinVsTalTTPForce,2)]) %pax vs vin
text(95,1.5,['p=' num2str(pPaxVsVinTTPForce,2)]) %tal vs vin
text(95,2,['p=' num2str(pPaxVsTalTTPForce,2)]) %pax vs tal
txt = findobj(gca,'Type','text');
set(txt(3:end),'VerticalAlignment', 'Middle');
set(findobj(gca,'LineStyle','--'),'LineStyle','-')
set(findobj(gca,'tag','Median'),'LineWidth',2)
xlabel('Time to peak in force (s)')
set(gca,'FontSize',7)
set(txt,'FontSize',7)
%% Assembly rate - loading
for ii=[1 2 3 4]
    curAssemblyRate = load([paxFolder{ii} '/data/nonTransmittingVsForceTransmitting.mat']);
    curAmpSlope = curAssemblyRate.avgEarlyAmpSlope_forceTransmitting;
    curForceSlope = curAssemblyRate.avgForceSlope_forceTransmitting;
    AmpSlopePax{ii} = curAmpSlope;
    ForceSlopePax{ii} = curForceSlope;
end
totalAmpSlopePax=[AmpSlopePax{1}; AmpSlopePax{2}; AmpSlopePax{3}; AmpSlopePax{4}];
totalForceSlopePax=[ForceSlopePax{1}; ForceSlopePax{2}; ForceSlopePax{3}; ForceSlopePax{4}];

for ii=1:6
    curAssemblyRateVin = load([vinFolder{ii} '/data/nonTransmittingVsForceTransmitting.mat']);
    curAmpSlope = curAssemblyRateVin.avgEarlyAmpSlope_forceTransmitting;
    curForceSlope = curAssemblyRateVin.avgForceSlope_forceTransmitting;
    AmpSlopeVin{ii} = curAmpSlope;
    ForceSlopeVin{ii} = curForceSlope;
end
totalAmpSlopeVin=[AmpSlopeVin{1}; AmpSlopeVin{2}; AmpSlopeVin{3}; 
    AmpSlopeVin{4}; AmpSlopeVin{5}; AmpSlopeVin{6}];
totalForceSlopeVin=[ForceSlopeVin{1}; ForceSlopeVin{2}; ForceSlopeVin{3}; 
    ForceSlopeVin{4}; ForceSlopeVin{5}; ForceSlopeVin{6}];

for ii=1:numel(talFolder)
    curAssemblyRateTal = load([talFolder{ii} '/data/nonTransmittingVsForceTransmitting.mat']);
    curAmpSlope = curAssemblyRateTal.avgEarlyAmpSlope_forceTransmitting;
    curForceSlope = curAssemblyRateTal.avgForceSlope_forceTransmitting;
    AmpSlopeTal{ii} = curAmpSlope;
    ForceSlopeTal{ii} = curForceSlope;
end
% totalAmpSlopeTal=[AmpSlopeTal{1}];
% totalForceSlopeTal=[ForceSlopeTal{1}];
totalAmpSlopeTal=[AmpSlopeTal{1}; AmpSlopeTal{2}; AmpSlopeTal{3}; 
    AmpSlopeTal{4}; AmpSlopeTal{5}; AmpSlopeTal{6}];
totalForceSlopeTal=[ForceSlopeTal{1}; ForceSlopeTal{2}; ForceSlopeTal{3}; 
    ForceSlopeTal{4}; ForceSlopeTal{5}; ForceSlopeTal{6}];

%% Statistical testing for slopes
[hNormalASlopePax,pNormalASlopePax] = kstest(totalAmpSlopePax);
[hNormalASlopeVin,pNormalASlopeVin] = kstest(totalAmpSlopeVin);
[hNormalASlopeTal,pNormalASlopeTal] = kstest(totalAmpSlopeTal);
% hNormalASlopePax = 1 -> It's not normal.
[pPaxVsVinASlope,hPaxVsVinASlope] = ranksum(totalAmpSlopePax,totalAmpSlopeVin);
[pPaxVsTalASlope,hPaxVsTalASlope] = ranksum(totalAmpSlopePax,totalAmpSlopeTal);
[pVinVsTalASlope,hVinVsTalASlope] = ranksum(totalAmpSlopeVin,totalAmpSlopeTal);

[hNormalFSlopeTal,pNormalFSlopeTal] = kstest(totalForceSlopeTal);
% hNormalFSlopeTal = 1 -> It's not normal.
[pPaxVsVinFSlope,hPaxVsVinFSlope] = ranksum(totalForceSlopePax,totalForceSlopeVin);
[pPaxVsTalFSlope,hPaxVsTalFSlope] = ranksum(totalForceSlopePax,totalForceSlopeTal);
[pVinVsTalFSlope,hVinVsTalFSlope] = ranksum(totalForceSlopeVin,totalForceSlopeTal);
%% Assembly rate - plotting
[lengthLongest,iLongest]=max([length(totalAmpSlopePax),length(totalAmpSlopeVin),length(totalAmpSlopeTal)]);
% it was vinculin that was longest
matrixAmpSlope = NaN(lengthLongest,3);
matrixAmpSlope(1:length(totalAmpSlopePax),3) = totalAmpSlopePax;
matrixAmpSlope(1:length(totalAmpSlopeVin),2) = totalAmpSlopeVin;
matrixAmpSlope(1:length(totalAmpSlopeTal),1) = totalAmpSlopeTal;

axes('Position',[400/600,330/600,100/600,70/600])
boxplot(matrixAmpSlope,'orientation','vertical','whisker',0.5,'notch','on',...
    'labels',{'Tal','Vin','Pax'},'symbol','','widths',boxWidth,'jitter',1,'colors','k')
% bar([mean(totalAmpSlopePax),mean(totalAmpSlopeVin),mean(totalAmpSlopeTal)]); hold on
% errorbar([mean(totalAmpSlopePax),mean(totalAmpSlopeVin),mean(totalAmpSlopeTal)],...
%     [std(totalAmpSlopePax)/sqrt(length(totalAmpSlopePax)),...
%     std(totalAmpSlopeVin)/sqrt(length(totalAmpSlopeVin)),...
%     std(totalAmpSlopeTal)/sqrt(length(totalAmpSlopeTal))],'LineStyle','none','Color','k');
% bar([median(totalAmpSlopePax),median(totalAmpSlopeVin),median(totalAmpSlopeTal)]); hold on
% errorbar([median(totalAmpSlopePax),median(totalAmpSlopeVin),median(totalAmpSlopeTal)],...
%     [std(totalAmpSlopePax)/sqrt(length(totalAmpSlopePax)),...
%     std(totalAmpSlopeVin)/sqrt(length(totalAmpSlopeVin)),...
%     std(totalAmpSlopeTal)/sqrt(length(totalAmpSlopeTal))],'LineStyle','none','Color','k');
% set(gca,'xticklabel',{'Tal';'Vin';'Pax'})
text(1,2300,['p=' num2str(pVinVsTalASlope,2)]) %pax vs vin
text(3,1500,['p=' num2str(pPaxVsVinASlope,2)]) %tal vs vin
text(2,2600,['p=' num2str(pPaxVsTalASlope,2)]) %pax vs tal
txt = findobj(gca,'Type','text');
set(txt(3:end),'VerticalAlignment', 'Middle');
set(findobj(gca,'LineStyle','--'),'LineStyle','-')
set(findobj(gca,'tag','Median'),'LineWidth',2)
% ylim([-20 180])
ylim([-1000 2800])
% xlabel('Time to peak (s)')
ylabel('F.I. growth (a.u./s)')
set(gca,'FontSize',7)
set(txt,'FontSize',7)

% Force
% totalForceSlopePax =arrayfun(@(x) (x.forceSlope),tracksNApax(idGroupsPax{1} | idGroupsPax{7}));
% totalForceSlopeVin =arrayfun(@(x) (x.forceSlope),tracksNAvin(idGroupsVin5{1} | idGroupsVin5{7}));
% totalForceSlopeTal =arrayfun(@(x) (x.forceSlope),tracksNAtal(idGroupsTal{1} | idGroupsTal{7}));

[lengthLongest,iLongest]=max([length(totalForceSlopePax),length(totalForceSlopeVin),length(totalForceSlopeTal)]);
% it was vinculin that was longest
matrixForceSlope = NaN(lengthLongest,3);
matrixForceSlope(1:length(totalForceSlopePax),3) = totalForceSlopePax;
matrixForceSlope(1:length(totalForceSlopeVin),2) = totalForceSlopeVin;
matrixForceSlope(1:length(totalForceSlopeTal),1) = totalForceSlopeTal;

axes('Position',[430/600,230/600,70/600,70/600])
% axes('Position',[0.91,0.06,0.08,0.13])
boxplot(matrixForceSlope,'orientation','vertical','whisker',0.5,'notch','on',...
    'labels',{'Tal','Vin','Pax'},'symbol','','widths',boxWidth,'jitter',1,'colors','k')
% bar([median(totalForceSlopePax),median(totalForceSlopeVin),median(totalForceSlopeTal)])
% bar([mean(totalForceSlopePax),mean(totalForceSlopeVin),mean(totalForceSlopeTal)]); hold on
% errorbar([mean(totalForceSlopePax),mean(totalForceSlopeVin),mean(totalForceSlopeTal)],...
%     [std(totalForceSlopePax)/sqrt(length(totalForceSlopePax)),...
%     std(totalForceSlopeVin)/sqrt(length(totalForceSlopeVin)),...
%     std(totalForceSlopeTal)/sqrt(length(totalForceSlopeTal))],'LineStyle','none','Color','k');
% set(gca,'xticklabel',{'Tal';'Vin';'Pax'})
% boxplot(matrixForceSlope,'orientation','vertical','whisker',1,'notch','on',...
%     'labels',{'Pax','Vin','Tal'},'symbol','','widths',boxWidth,'jitter',1,'colors','k')
set(findobj(gca,'LineStyle','--'),'LineStyle','-')
set(findobj(gca,'tag','Median'),'LineWidth',2)
ylim([-33 70])
ylabel('Traction growth (Pa/s)')
set(gca,'FontSize',7)
text(1,43,['p=' num2str(pVinVsTalFSlope,2)]) %pax vs vin
text(3,48,['p=' num2str(pPaxVsVinFSlope,2)]) %tal vs vin
text(2,60,['p=' num2str(pPaxVsTalFSlope,2)]) %pax vs tal
txt = findobj(gca,'Type','text');
set(txt(3:end),'VerticalAlignment', 'Middle');
set(txt,'FontSize',7)

% export_fig([f2path filesep 'Fig2TimeLags.eps']);% histogramPeakLagVinVsTal -transparent
%% save
print('-depsc', '-loose', [f3_Tmeasure_path filesep 'Fig2TimeLags.eps'])
print('-dtiff', '-loose', '-r300', [f3_Tmeasure_path filesep 'Fig2TimeLags.tif'])
savefig([f3_Tmeasure_path filesep 'Fig2TimeLags.fig'])

% title('Time lag in peaks')
%% Cross-cell validation figure with vin4 vs. vin5 - I did this in different m file
% % vin4
% v4i1=load([vinFolder{4} filesep 'data' filesep 'selectedGroups.mat']);
% idGroupSelected4i1={v4i1.idGroup1Selected,v4i1.idGroup2Selected,v4i1.idGroup3Selected,v4i1.idGroup4Selected,v4i1.idGroup5Selected,....
%                                     v4i1.idGroup6Selected, v4i1.idGroup7Selected,v4i1.idGroup8Selected,v4i1.idGroup9Selected};
% [Ti1,allData_i1]=extractFeatureNA(tracksNAvin{4},idGroupSelected4i1);
% % add vin4 with int4sec
% vinFolder4i4='/storage/disk3/NA_recruitment/analysis/Sangyoon/NA_recruitment/Alexia/2015-04-29/cell2/Int4sec/Colocalization/analysis2';
% v4i4=load([vinFolder4i4 filesep 'data' filesep 'selectedGroups.mat']);
% idGroupSelected4i4={v4i4.idGroup1Selected,v4i4.idGroup2Selected,v4i4.idGroup3Selected,v4i4.idGroup4Selected,v4i4.idGroup5Selected,....
%                                     v4i4.idGroup6Selected, v4i4.idGroup7Selected,v4i4.idGroup8Selected,v4i4.idGroup9Selected};
% curtracksNAvin = load([vinFolder4i4 filesep 'data' filesep 'tracksNA.mat'],'tracksNA');
% tracksNAvin4i4 = curtracksNAvin.tracksNA; 
% [Ti4,allData_i4]=extractFeatureNA(tracksNAvin4i4,idGroupSelected4i4);
% % add i1 analysis1
% vinFolder4i1a1='/storage/disk3/NA_recruitment/analysis/Sangyoon/NA_recruitment/Alexia/2015-04-29/cell2/Int1sec/Colocalization/analysis1';
% v4i1a1=load([vinFolder4i1a1 filesep 'data' filesep 'selectedGroups.mat']);
% idGroupSelected4i1a1={v4i1a1.idGroup1Selected,v4i1a1.idGroup2Selected,v4i1a1.idGroup3Selected,v4i1a1.idGroup4Selected,v4i1a1.idGroup5Selected,....
%                                     v4i1a1.idGroup6Selected, v4i1a1.idGroup7Selected,v4i1a1.idGroup8Selected,v4i1a1.idGroup9Selected};
% curtracksNAvin = load([vinFolder4i1a1 filesep 'data' filesep 'tracksNA.mat'],'tracksNA');
% tracksNAvin4i1a1 = curtracksNAvin.tracksNA; 
% [Ti1a1,allData_i1a1]=extractFeatureNA(tracksNAvin4i1a1,idGroupSelected4i1a1);
% 
% Tv4=[Ti4; Ti1];%Ti1a1]; % Merge
% [trainedClassifier_v4, validationAccuracy_v4, Cv4, orderV4] = trainClassifierNA(Tv4);
% disp(['Validation accuracy is ' num2str(validationAccuracy_v4) '.'])
% % figured out that I need label specifically 1sec v4 video
% [idGroup1,idGroup2,idGroup3,idGroup4,idGroup5,idGroup6,idGroup7,idGroup8,idGroup9]= classifyNascentAdhesionTracks(vinFolder{4},'tracksNA',tracksNAvin{4});
% % Now for vin5
% [idGroup1,idGroup2,idGroup3,idGroup4,idGroup5,idGroup6,idGroup7,idGroup8,idGroup9]= classifyNascentAdhesionTracks(vinFolder{5},'tracksNA',tracksNAvin{5});
% %% Confusion matrix for v4 to v4
%     % normalize confusion matrix
%     for ii=1:size(C,1)
%         C(ii,:) = C(ii,:)/sum(C(ii,:));
%     end
%     response = T.Group;
%     % Get the unique resonses
%     totalGroups = unique(response);
% 
%     figure; confAxis=axes; imagesc(C); title('Confusion Matrix')
%     set(confAxis,'xticklabel',totalGroups')
%     set(confAxis,'yticklabel',totalGroups')
%     c = colorbar;
%     c.Label.String = 'normalized prediction';
%     print('-depsc2', '-r300', [pathForColocalization filesep 'eps' filesep 'confusionMatrix.eps']);

%% Fig 3 Relationship with Non-force-transmitting adhesions and maturing adhesions
% close all
f3path = '/home2/shan/ownCloud/documents/NascentAdhesionForce/Manuscript/Figures/Fig3NonTransmittingAndMaturing';
if ~exist(f3path,'dir')
    mkdir(f3path)
end
%% Assembly and FI - non-force-transmitting - loading
for ii=1:4
    try
        curNTvsFT = load([paxFolder{ii} '/data/nonTransmittingVsForceTransmitting.mat']);
        curAvgEarlyAmpSlopeFTPax = curNTvsFT.avgEarlyAmpSlope_forceTransmitting;
        curAvgEarlyAmpSlopeNTPax = curNTvsFT.avgEarlyAmpSlope_nonTransmitting;
        AvgEarlyAmpSlopeFTPax{ii} = curAvgEarlyAmpSlopeFTPax;
        AvgEarlyAmpSlopeNTPax{ii} = curAvgEarlyAmpSlopeNTPax;
        
        curAvgPreFI_FTPax = curNTvsFT.avgPreDetecFluoInten_forceTransmitting;
        curAvgPreFI_NTPax = curNTvsFT.avgPreDetecFluoInten_nonTransmitting;
        AvgPreFI_FTPax{ii} = curAvgPreFI_FTPax;
        AvgPreFI_NTPax{ii} = curAvgPreFI_NTPax;
    catch
        disp([paxFolder{ii} ' has not been analyzed for nonTransmitting.'])
    end
end
totalAvgEarlyAmpSlopeFTPax=[AvgEarlyAmpSlopeFTPax{1}; AvgEarlyAmpSlopeFTPax{2}; 
    AvgEarlyAmpSlopeFTPax{3}; AvgEarlyAmpSlopeFTPax{4}];
totalAvgEarlyAmpSlopeNTPax=[AvgEarlyAmpSlopeNTPax{1}; AvgEarlyAmpSlopeNTPax{2}; 
    AvgEarlyAmpSlopeNTPax{3}; AvgEarlyAmpSlopeNTPax{4}];
totalAvgPreFI_FTPax=[AvgPreFI_FTPax{1}; AvgPreFI_FTPax{2}; 
    AvgPreFI_FTPax{3}; AvgPreFI_FTPax{4}];
totalAvgPreFI_NTPax=[AvgPreFI_NTPax{1}; AvgPreFI_NTPax{2}; 
    AvgPreFI_NTPax{3}; AvgPreFI_NTPax{4}];

for ii=1:6
    try
        curNTvsFT = load([vinFolder{ii} '/data/nonTransmittingVsForceTransmitting.mat']);
        curAvgEarlyAmpSlopeFTVin = curNTvsFT.avgEarlyAmpSlope_forceTransmitting;
        curAvgEarlyAmpSlopeNTVin = curNTvsFT.avgEarlyAmpSlope_nonTransmitting;
        AvgEarlyAmpSlopeFTVin{ii} = curAvgEarlyAmpSlopeFTVin;
        AvgEarlyAmpSlopeNTVin{ii} = curAvgEarlyAmpSlopeNTVin;
        
        curAvgPreFI_FTVin = curNTvsFT.avgPreDetecFluoInten_forceTransmitting;
        curAvgPreFI_NTVin = curNTvsFT.avgPreDetecFluoInten_nonTransmitting;
        AvgPreFI_FTVin{ii} = curAvgPreFI_FTVin;
        AvgPreFI_NTVin{ii} = curAvgPreFI_NTVin;
    catch
        disp([vinFolder{ii} ' has not been analyzed for nonTransmitting.'])
    end
end
totalAvgEarlyAmpSlopeFTVin=[AvgEarlyAmpSlopeFTVin{5}; AvgEarlyAmpSlopeFTVin{6}];
totalAvgEarlyAmpSlopeNTVin=[AvgEarlyAmpSlopeNTVin{5}; AvgEarlyAmpSlopeNTVin{6}];
totalAvgPreFI_FTVin=[AvgPreFI_FTVin{5}; AvgPreFI_FTVin{6}];
totalAvgPreFI_NTVin=[AvgPreFI_NTVin{5}; AvgPreFI_NTVin{6}];

for ii=1:7
    try
        curNTvsFT = load([talFolder{ii} '/data/nonTransmittingVsForceTransmitting.mat']);
        curAvgEarlyAmpSlopeFTTal = curNTvsFT.avgEarlyAmpSlope_forceTransmitting;
        curAvgEarlyAmpSlopeNTTal = curNTvsFT.avgEarlyAmpSlope_nonTransmitting;
        AvgEarlyAmpSlopeFTTal{ii} = curAvgEarlyAmpSlopeFTTal;
        AvgEarlyAmpSlopeNTTal{ii} = curAvgEarlyAmpSlopeNTTal;
        
        curAvgPreFI_FTTal = curNTvsFT.avgPreDetecFluoInten_forceTransmitting;
        curAvgPreFI_NTTal = curNTvsFT.avgPreDetecFluoInten_nonTransmitting;
        AvgPreFI_FTTal{ii} = curAvgPreFI_FTTal;
        AvgPreFI_NTTal{ii} = curAvgPreFI_NTTal;
    catch
        disp([talFolder{ii} ' has not been analyzed for nonTransmitting.'])
    end
end
totalAvgEarlyAmpSlopeFTTal=[AvgEarlyAmpSlopeFTTal{1}];
totalAvgEarlyAmpSlopeNTTal=[AvgEarlyAmpSlopeNTTal{1}];
totalAvgPreFI_FTTal=[AvgPreFI_FTTal{1}];
totalAvgPreFI_NTTal=[AvgPreFI_NTTal{1}];
%% Statistical testing for NT slopes
[hNormalPreFI_FTPax,pNormalPreFI_FTPax] = kstest(totalAvgPreFI_FTPax);
% hNormalPreFI_FTPax = 1 -> It's not normal.
[pPaxPreFI_FTvsNT,hPreFI_PaxFTvsNT] = ranksum(totalAvgPreFI_FTPax,totalAvgPreFI_NTPax);
[pVinPreFI_FTvsNT,hPreFI_VinFTvsNT] = ranksum(totalAvgPreFI_FTVin,totalAvgPreFI_NTVin);
[pTalPreFI_FTvsNT,hPreFI_TalFTvsNT] = ranksum(totalAvgPreFI_FTTal,totalAvgPreFI_NTTal);

[hNormalEarlyAmpSlopeFTPax,pNormalEarlyAmpSlopeFTPax] = kstest(totalAvgEarlyAmpSlopeFTPax);
% hNormalEarlyAmpSlopeFTPax = 1 -> It's not normal.
[pPaxAmpFTvsNT,hAmpPaxFTvsNT] = ranksum(totalAvgEarlyAmpSlopeFTPax,totalAvgEarlyAmpSlopeNTPax);
[pVinAmpFTvsNT,hAmpVinFTvsNT] = ranksum(totalAvgEarlyAmpSlopeFTVin,totalAvgEarlyAmpSlopeNTVin);
[pTalAmpFTvsNT,hAmpTalFTvsNT] = ranksum(totalAvgEarlyAmpSlopeFTTal,totalAvgEarlyAmpSlopeNTTal);

%% Assembly and FI - non-force-transmitting - loading
pos(3:4) = [400 600];
figure('Position', pos, 'PaperPositionMode', 'auto','Color','w');
axes('Position', [30/400 520/600 100/400 60/600]); % 
% axes('Position', [0.07 .7 .18 .22]); % 

[lengthLongest,iLongest]=max([length(totalAvgPreFI_FTPax),length(totalAvgPreFI_NTPax),...
    length(totalAvgPreFI_FTVin), length(totalAvgPreFI_NTVin),...
    length(totalAvgPreFI_FTTal), length(totalAvgPreFI_NTTal)]);
% it was vinculin that was longest
matrixAvgPreFI = NaN(lengthLongest,6);
matrixAvgPreFI(1:length(totalAvgPreFI_FTTal),1) = totalAvgPreFI_FTTal;
matrixAvgPreFI(1:length(totalAvgPreFI_NTTal),2) = totalAvgPreFI_NTTal;
matrixAvgPreFI(1:length(totalAvgPreFI_FTVin),3) = totalAvgPreFI_FTVin;
matrixAvgPreFI(1:length(totalAvgPreFI_NTVin),4) = totalAvgPreFI_NTVin;
matrixAvgPreFI(1:length(totalAvgPreFI_FTPax),5) = totalAvgPreFI_FTPax;
matrixAvgPreFI(1:length(totalAvgPreFI_NTPax),6) = totalAvgPreFI_NTPax;

boxplot(matrixAvgPreFI,'orientation','vertical','whisker',0.5,'notch','off',...
    'labels',{'Tal FT','Tal NT','Vin FT','Vin NT','Pax FT','Pax NT'},'symbol','','widths',boxWidth,'jitter',1,'colors','k',...
    'labelorientation','inline')
txt = findobj(gca,'Type','text');
set(txt(3:end),'VerticalAlignment', 'Middle');
set(findobj(gca,'LineStyle','--'),'LineStyle','-')
set(findobj(gca,'tag','Median'),'LineWidth',2)
ylim([0 1000])
ylabel('Initial F.I. (a.u.)')
set(gca,'FontSize',7)
set(txt,'FontSize',7)
title('non-transmitting NAs')
%% Assembly
axes('Position', [160/400 520/600 100/400 60/600]); % 
[lengthLongest,iLongest]=max([length(totalAvgEarlyAmpSlopeFTTal),length(totalAvgEarlyAmpSlopeNTTal),...
    length(totalAvgEarlyAmpSlopeFTVin), length(totalAvgEarlyAmpSlopeNTVin),...
    length(totalAvgEarlyAmpSlopeFTPax), length(totalAvgEarlyAmpSlopeNTPax)]);
% it was vinculin that was longest
matrixAvgEarlyAmpSlope = NaN(lengthLongest,6);
matrixAvgEarlyAmpSlope(1:length(totalAvgEarlyAmpSlopeFTTal),1) = totalAvgEarlyAmpSlopeFTTal;
matrixAvgEarlyAmpSlope(1:length(totalAvgEarlyAmpSlopeNTTal),2) = totalAvgEarlyAmpSlopeNTTal;
matrixAvgEarlyAmpSlope(1:length(totalAvgEarlyAmpSlopeFTVin),3) = totalAvgEarlyAmpSlopeFTVin;
matrixAvgEarlyAmpSlope(1:length(totalAvgEarlyAmpSlopeNTVin),4) = totalAvgEarlyAmpSlopeNTVin;
matrixAvgEarlyAmpSlope(1:length(totalAvgEarlyAmpSlopeFTPax),5) = totalAvgEarlyAmpSlopeFTPax;
matrixAvgEarlyAmpSlope(1:length(totalAvgEarlyAmpSlopeNTPax),6) = totalAvgEarlyAmpSlopeNTPax;

boxplot(matrixAvgEarlyAmpSlope,'orientation','vertical','whisker',0.5,'notch','off',...
    'labels',{'Tal FT','Tal NT','Vin FT','Vin NT','Pax FT','Pax NT'},'symbol','','widths',boxWidth,'jitter',1,'colors','k','labelorientation','inline')
txt = findobj(gca,'Type','text');
set(txt(3:end),'VerticalAlignment', 'Middle');
set(findobj(gca,'LineStyle','--'),'LineStyle','-')
set(findobj(gca,'tag','Median'),'LineWidth',2)
ylim([-1000 1500])
ylabel('Assembly rate (a.u./sec)')
set(gca,'FontSize',7)
set(txt,'FontSize',7)
title('non-transmitting NAs')
%% plotting the time series
% for jj=1:numel(tracksNAtalG2)
%     curTracks=tracksNAtalG2{jj};
%     if ~isempty(curTracks)
%         if jj==1
%             tIntTal=1.6;
%         else
%             tIntTal=2;
%         end
%         figure
%         subplot(2,1,1)
%         plotIntensityForce(curTracks,[],false,false,'UseCurrentAxis',true,...
%             'Source',{'ampTotal'},'plotCohorts',true,'tInterval',tIntTal,'prePostFrames',10,'yNormalization',true)
%         title(['Intensity ' num2str(jj)])
%         subplot(2,1,2)
%         plotIntensityForce(curTracks,[],false,false,'UseCurrentAxis',true,...
%             'Source',{'forceMag'},'plotCohorts',true,'tInterval',tIntTal,'prePostFrames',10,'yNormalization',true)
%         title(['Force ' num2str(jj)])
%     end
% end
%% Quickly check assembly rate in vinculin at the time of time zero for g1 and g2

%% time course - force registration
% % Get all tracks
% if ~exist('tracksNAvin5G1','var')
%     disp('loading tracksNAvin5G1 ...'); tic;
%     tracksNAvin5G1 = load([vinFolder{5} filesep 'data' filesep 'tracksG1.mat'],'tracksG1');
%     tracksNAvin5G1 = tracksNAvin5G1.tracksG1; 
%     tracksNAvin5G2 = load([vinFolder{5} filesep 'data' filesep 'tracksG2.mat'],'tracksG2');
%     tracksNAvin5G2 = tracksNAvin5G2.tracksG2; 
%     tracksNAvin5G7 = load([vinFolder{5} filesep 'data' filesep 'tracksG7.mat'],'tracksG7');
%     tracksNAvin5G7 = tracksNAvin5G7.tracksG7; 
%     toc
% end
% if ~exist('tracksNAtal1G1','var')
%     disp('loading tracksNAtal1G1 ...'); tic;
%     tracksNAtal1G1 = load([talFolder{1} filesep 'data' filesep 'tracksG1.mat'],'tracksG1');
%     tracksNAtal1G1 = tracksNAtal1G1.tracksG1; 
%     tracksNAtal1G2 = load([talFolder{1} filesep 'data' filesep 'tracksG2.mat'],'tracksG2');
%     tracksNAtal1G2 = tracksNAtal1G2.tracksG2; 
%     tracksNAtal1G7 = load([talFolder{1} filesep 'data' filesep 'tracksG7.mat'],'tracksG7');
%     tracksNAtal1G7 = tracksNAtal1G7.tracksG7; 
%     toc; 
% end
% if ~exist('tracksNApax4G1','var')
%     disp('loading tracksNApax4G1 ...'); tic;
%     tracksNApax4G1 = load([paxFolder{4} filesep 'data' filesep 'tracksG1.mat'],'tracksG1');
%     tracksNApax4G1 = tracksNApax4G1.tracksG1; 
%     tracksNApax4G2 = load([paxFolder{4} filesep 'data' filesep 'tracksG2.mat'],'tracksG2');
%     tracksNApax4G2 = tracksNApax4G2.tracksG2; 
%     tracksNApax4G7 = load([paxFolder{4} filesep 'data' filesep 'tracksG7.mat'],'tracksG7');
%     tracksNApax4G7 = tracksNApax4G7.tracksG7; 
%     toc; 
% end
% if ~exist('tracksNApax3G1','var')
%     disp('loading tracksNApax3G1 ...'); tic;
%     tracksNApax3G1 = load([paxFolder{3} filesep 'data' filesep 'tracksG1.mat'],'tracksG1');
%     tracksNApax3G1 = tracksNApax3G1.tracksG1; 
%     tracksNApax3G2 = load([paxFolder{3} filesep 'data' filesep 'tracksG2.mat'],'tracksG2');
%     tracksNApax3G2 = tracksNApax3G2.tracksG2; 
%     tracksNApax3G7 = load([paxFolder{3} filesep 'data' filesep 'tracksG7.mat'],'tracksG7');
%     tracksNApax3G7 = tracksNApax3G7.tracksG7; 
%     toc; 
% end
% if ~exist('tracksNApax2G1','var')
%     disp('loading tracksNApax2G1 ...'); tic;
%     tracksNApax2G1 = load([paxFolder{2} filesep 'data' filesep 'tracksG1.mat'],'tracksG1');
%     tracksNApax2G1 = tracksNApax2G1.tracksG1; 
%     tracksNApax2G2 = load([paxFolder{2} filesep 'data' filesep 'tracksG2.mat'],'tracksG2');
%     tracksNApax2G2 = tracksNApax2G2.tracksG2; 
%     tracksNApax2G7 = load([paxFolder{2} filesep 'data' filesep 'tracksG7.mat'],'tracksG7');
%     tracksNApax2G7 = tracksNApax2G7.tracksG7; 
%     toc; 
% end
% 
% %% 1. Get the late force slopes from G1, G2 and G7
% [tracksNA, lateForceSlopes] = calculateTrackLateForceSlope(tracksNA,tInterval);
% 
% %% 2. Re-separate G1 and G2 based on late force slopes
% 
% axes('Position', [50/400 400/600 140/400 60/600]); % 
% plotIntensityForce(tracksNAvin((idGroupsVin5{1} & forceTransmittingVin)),[],false,false,'UseCurrentAxis',true,...
%     'Source',{'forceMag'},'plotCohorts',true,'tInterval',MDpax.timeInterval_,'prePostFrames',10,'yNormalization',true)
% plotIntensityForce(tracksNAvin((idGroupsVin5{1} & ~forceTransmittingVin)),[],false,false,'UseCurrentAxis',true,...
%     'Source',{'forceMag'},'plotCohorts',true,'tInterval',MDpax.timeInterval_,'prePostFrames',10,'yNormalization',true)
% plotIntensityForce(tracksNAvin((idGroupsVin5{2})),[],false,false,'UseCurrentAxis',true,...
%     'Source',{'forceMag'},'plotCohorts',true,'tInterval',MDpax.timeInterval_,'prePostFrames',10,'yNormalization',true)
% title('Force')
% set(gca,'FontSize',7)
% axes('Position', [50/400 300/600 140/400 60/600]); % 
% title('Talin')
% set(gca,'FontSize',7)
% axes('Position', [50/400 200/600 140/400 60/600]); % 
% title('Vinculin')
% set(gca,'FontSize',7)
% axes('Position', [50/400 100/600 140/400 60/600]); % 
% title('Paxillin')
% set(gca,'FontSize',7)
% %% assembly rate, force growth rate - loading
% for ii=1:4
%     try
%         curAssemblyRate = load([paxFolder{ii} '/data/earlyAmpSlope.mat']);
%         curAmpSlope = curAssemblyRate.earlyAmpSlope;
%         AmpSlopePaxG2{ii} = curAmpSlope{2};
%         curForceGrowth = load([paxFolder{ii} '/data/forceSlope.mat']);
%         curForceSlope = curForceGrowth.forceSlope;
%         ForceSlopePaxG2{ii} = curForceSlope{2};
%     catch
%         disp([paxFolder{ii} ' has not been analyzed for earlyAmpSlope or forceSlope.'])
%     end
% end
% totalAmpSlopePaxG2=[AmpSlopePaxG2{1}; AmpSlopePaxG2{2}; AmpSlopePaxG2{3}; AmpSlopePaxG2{4}];
% totalForceSlopePaxG2=[ForceSlopePaxG2{1}; ForceSlopePaxG2{2}; ForceSlopePaxG2{3}; ForceSlopePaxG2{4}];
% 
% for ii=1:6
%     try
%         curAssemblyRate = load([vinFolder{ii} '/data/earlyAmpSlope.mat']);
%         curAmpSlope = curAssemblyRate.earlyAmpSlope;
%         AmpSlopeVinG2{ii} = curAmpSlope{2};
%         curForceGrowth = load([vinFolder{ii} '/data/forceSlope.mat']);
%         curForceSlope = curForceGrowth.forceSlope;
%         ForceSlopeVinG2{ii} = curForceSlope{2};
%     catch
%         disp([vinFolder{ii} ' has not been analyzed for earlyAmpSlope or forceSlope.'])
%     end
% end
% totalAmpSlopeVinG2=[AmpSlopeVinG2{1}; AmpSlopeVinG2{2}; AmpSlopeVinG2{3}; 
%     AmpSlopeVinG2{4}; AmpSlopeVinG2{5}; AmpSlopeVinG2{6}];
% totalForceSlopeVinG2=[ForceSlopeVinG2{1}; ForceSlopeVinG2{2}; ForceSlopeVinG2{3}; 
%     ForceSlopeVinG2{4}; ForceSlopeVinG2{5}; ForceSlopeVinG2{6}];
% 
% for ii=1:1
%     try
%         curAssemblyRate = load([talFolder{ii} '/data/earlyAmpSlope.mat']);
%         curAmpSlope = curAssemblyRate.earlyAmpSlope;
%         AmpSlopeTalG2{ii} = curAmpSlope{2};
%         curForceGrowth = load([talFolder{ii} '/data/forceSlope.mat']);
%         curForceSlope = curForceGrowth.forceSlope;
%         ForceSlopeTalG2{ii} = curForceSlope{2};
%     catch
%         disp([talFolder{ii} ' has not been analyzed for earlyAmpSlope or forceSlope.'])
%     end
% end
% totalAmpSlopeTalG2 = [AmpSlopeTalG2{1}];
% totalForceSlopeTalG2 = [ForceSlopeTalG2{1}];
% %% assembly rate - plotting
% % axes('Position',[0.5,0.08,0.40,0.29])
% 
% [lengthLongest,iLongest]=max([length(totalAmpSlopeTalG2),length(totalAvgEarlyAmpSlopeFTTal),...
%     length(totalAmpSlopeVinG2), length(totalAvgEarlyAmpSlopeFTVin),...
%     length(totalAmpSlopePaxG2), length(totalAvgEarlyAmpSlopeFTPax)]);
% % it was vinculin that was longest
% matrixAvgEarlyAmpSlopeG2 = NaN(lengthLongest,6);
% matrixAvgEarlyAmpSlopeG2(1:length(totalAmpSlopeTalG2),1) = totalAmpSlopeTalG2;
% matrixAvgEarlyAmpSlopeG2(1:length(totalAvgEarlyAmpSlopeFTTal),2) = totalAvgEarlyAmpSlopeFTTal;
% matrixAvgEarlyAmpSlopeG2(1:length(totalAmpSlopeVinG2),3) = totalAmpSlopeVinG2;
% matrixAvgEarlyAmpSlopeG2(1:length(totalAvgEarlyAmpSlopeFTVin),4) = totalAvgEarlyAmpSlopeFTVin;
% matrixAvgEarlyAmpSlopeG2(1:length(totalAmpSlopePaxG2),5) = totalAmpSlopePaxG2;
% matrixAvgEarlyAmpSlopeG2(1:length(totalAvgEarlyAmpSlopeFTPax),6) = totalAvgEarlyAmpSlopeFTPax;
% 
% % axes('Position',[0.08,0.45,0.40,0.18])
% axes('Position', [250/400 240/600 100/400 100/600]); % 
% boxplot(matrixAvgEarlyAmpSlopeG2,'orientation','vertical','whisker',0.5,'notch','on',...
%     'labels',{'Tal G2','Tal FT','Vin G2','Vin FT','Pax G2','Pax FT'},'symbol','','widths',boxWidth,'jitter',1,'colors','k', 'labelorientation','inline')
% txt = findobj(gca,'Type','text');
% set(txt(3:end),'VerticalAlignment', 'Middle');
% set(findobj(gca,'LineStyle','--'),'LineStyle','-')
% set(findobj(gca,'tag','Median'),'LineWidth',2)
% ylim([-500 2000])
% ylabel('Assembly rate (a.u./sec)')
% set(findobj(gca,'Type','Text'),'FontSize',7)
% set(gca,'FontSize',7)
%% save
print('-depsc', '-loose', [f3path filesep 'Fig3NTG1_G2.eps'])
print('-dtiff', '-loose', '-r300', [f3path filesep 'Fig3NTG1_G2.tif'])
savefig([f3path filesep 'Fig3NTG1_G2.fig'])
%% The plot above was obviously wrong - mainly because the slopes are calculated based on 'startingFrame'. 
% % The slopes should've been calculated based on 'startingFrameExtra'. I
% % will do the slope calculation again and plot them again. - sh 2/29/16 
% % Let's do it with existing ones: tracksNAvin(idGroupsVin5{1} & forceTransmittingVin)
% tracksNAvinG1FT=tracksNAvin(idGroupsVin5{1} & forceTransmittingVin);
% [~,earlyVinSlopeG1FT] = arrayfun(@(x) regression(1:3,x.ampTotal(x.startingFrameExtra:x.startingFrameExtra+2)),tracksNAvinG1FT);
% tracksNAvinG2FT=tracksNAvin(idGroupsVin5{2} & forceTransmittingVin);
% % We need to filter tracksNAvinG2FT out using 'state'
% % If there is no FC or FA in the state over the lifetime, they will be
% % filtered out.
% idxFAcontainingG2=arrayfun(@(x) any(strcmp(x.state,'FC') | strcmp(x.state,'FA')),tracksNAvinG2FT);
% tracksNAvinG2FT2=tracksNAvinG2FT(idxFAcontainingG2);
% % tracksNAvinG2FT=tracksNAvin5G2;
% [~,earlyVinSlopeG2FT] = arrayfun(@(x) regression(1:3,x.ampTotal(x.startingFrameExtra:x.startingFrameExtra+2)),tracksNAvinG2FT2);
% % lifeTimeVinG1FT=arrayfun(@(x) x.lifeTime,tracksNAvinG1FT);
% % lifeTimeVinG2FT=arrayfun(@(x) x.lifeTime,tracksNAvinG2FT2);
% % I have to do the same thing for talin
% tracksNAtal1G1FT=tracksNAtal(idGroupsTal{1} & forceTransmittingTal);
% [~,earlyTalSlopeG1FT] = arrayfun(@(x) regression(1:3,x.ampTotal(x.startingFrameExtra:x.startingFrameExtra+2)),tracksNAtal1G1FT);
% tracksNAtalG2FT=tracksNAtal(idGroupsTal{2} & forceTransmittingTal);
% idxFAcontainingG2tal=arrayfun(@(x) any(strcmp(x.state,'FC') | strcmp(x.state,'FA')),tracksNAtalG2FT);
% tracksNAtalG2FT2=tracksNAtalG2FT(idxFAcontainingG2tal);
% [~,earlyTalSlopeG2FT] = arrayfun(@(x) regression(1:3,x.ampTotal(x.startingFrameExtra:x.startingFrameExtra+2)),tracksNAtalG2FT2);
% % I have to do the same thing for paxillin
% tracksNApax1G1FT=tracksNApax(idGroupsPax{1} & forceTransmittingPax);
% [~,earlyPaxSlopeG1FT] = arrayfun(@(x) regression(1:3,x.ampTotal(x.startingFrameExtra:x.startingFrameExtra+2)),tracksNApax1G1FT);
% tracksNApaxG2FT=tracksNApax(idGroupsPax{2} & forceTransmittingPax);
% idxFAcontainingG2pax=arrayfun(@(x) any(strcmp(x.state,'FC') | strcmp(x.state,'FA')),tracksNApaxG2FT);
% tracksNApaxG2FT2=tracksNApaxG2FT(idxFAcontainingG2pax);
% [~,earlyPaxSlopeG2FT] = arrayfun(@(x) regression(1:3,x.ampTotal(x.startingFrameExtra:x.startingFrameExtra+2)),tracksNApaxG2FT2);
% % collect from one more cell
% tic
% tracksNAvin4 = load([vinFolder{4} filesep 'data' filesep 'tracksNA.mat'],'tracksNA');
% tracksNAvin4 = tracksNAvin4.tracksNA; toc
% idGroupsVin4struct = load([vinFolder{4} filesep 'data' filesep 'idsClassified.mat'],'idGroup1filtered','idGroup2filtered','idGroup3filtered',...
% 'idGroup4filtered','idGroup5filtered','idGroup6filtered','idGroup7filtered','idGroup8filtered','idGroup9filtered',...
% 'idGroup1','idGroup2','idGroup3',...
% 'idGroup4','idGroup5','idGroup6','idGroup7','idGroup8','idGroup9');
% idGroupsVin4{1} = idGroupsVin4struct.idGroup1filtered;
% idGroupsVin4{2} = idGroupsVin4struct.idGroup2filtered;
% idGroupsVin4{3} = idGroupsVin4struct.idGroup3filtered;
% idGroupsVin4{4} = idGroupsVin4struct.idGroup4filtered;
% idGroupsVin4{5} = idGroupsVin4struct.idGroup5filtered;
% idGroupsVin4{6} = idGroupsVin4struct.idGroup6filtered;
% idGroupsVin4{7} = idGroupsVin4struct.idGroup7filtered;
% idGroupsVin4{8} = idGroupsVin4struct.idGroup8filtered;
% idGroupsVin4{9} = idGroupsVin4struct.idGroup9filtered;
forceTransmittingVin4 = arrayfun(@(x) ~isempty(x.forceTransmitting) & x.forceTransmitting,tracksNAvin4,'UniformOutput',false);
idEmptyFTvin4=cellfun(@isempty,forceTransmittingVin4);
forceTransmittingVin4(idEmptyFTvin4)={false};
forceTransmittingVin4 = cell2mat(forceTransmittingVin4);

tracksNAvin4G1FT=tracksNAvin4(idGroupsVin4{1} & forceTransmittingVin4);
[~,earlyVin4SlopeG1FT] = arrayfun(@(x) regression(1:3,x.ampTotal(x.startingFrameExtra:x.startingFrameExtra+2)),tracksNAvin4G1FT);
tracksNAvin4G2FT=tracksNAvin4(idGroupsVin4{2} & forceTransmittingVin4);
idxFAcontainingVin4G2=arrayfun(@(x) any(strcmp(x.state,'FC') | strcmp(x.state,'FA')),tracksNAvin4G2FT);
tracksNAvin4G2FT2=tracksNAvin4G2FT(idxFAcontainingVin4G2);
[~,earlyVin4SlopeG2FT] = arrayfun(@(x) regression(1:5,x.ampTotal(x.startingFrameExtra:x.startingFrameExtra+4)),tracksNAvin4G2FT2);

%% distribution checking on vin slope
% [nPax,edges] = histcounts(earlyVinSlopeG1FT,60);
% figure, 
% h=histogram(earlyVinSlopeG1FT,edges,'Normalization','probability');
% hold on
% % [nVin] = histcounts(totalTLagVin,edges);
% histogram(earlyVinSlopeG2FT,edges,'Normalization','probability');
% %% distribution checking on tal slope
% [~,edges] = histcounts(earlyTalSlopeG1FT,60);
% figure, 
% histogram(earlyTalSlopeG1FT,edges,'Normalization','probability');
% hold on
% % [nVin] = histcounts(totalTLagVin,edges);
% histogram(earlyTalSlopeG2FT,edges,'Normalization','probability');
% %% plot about vin slope and tal slope
% earlyVinSlopeG1FTtot=[earlyVinSlopeG1FT; earlyVin4SlopeG1FT];
% earlyVinSlopeG2FTtot=[earlyVinSlopeG2FT; earlyVin4SlopeG2FT];
% [lengthLongest,iLongest]=max([length(earlyTalSlopeG2FT),length(earlyTalSlopeG1FT),...
%     length(earlyVinSlopeG2FTtot), length(earlyVinSlopeG1FTtot),...
%     length(earlyPaxSlopeG2FT), length(earlyPaxSlopeG1FT)]);
% % it was vinculin that was longest
% matrixEarlyAmpSlopeG1G2 = NaN(lengthLongest,6);
% matrixEarlyAmpSlopeG1G2(1:length(earlyTalSlopeG1FT),1) = earlyTalSlopeG1FT;
% matrixEarlyAmpSlopeG1G2(1:length(earlyTalSlopeG2FT),2) = earlyTalSlopeG2FT;
% matrixEarlyAmpSlopeG1G2(1:length(earlyVinSlopeG1FTtot),3) = earlyVinSlopeG1FTtot;
% matrixEarlyAmpSlopeG1G2(1:length(earlyVinSlopeG2FTtot),4) = earlyVinSlopeG2FTtot;
% matrixEarlyAmpSlopeG1G2(1:length(earlyPaxSlopeG1FT),5) = earlyPaxSlopeG1FT;
% matrixEarlyAmpSlopeG1G2(1:length(earlyPaxSlopeG2FT),6) = earlyPaxSlopeG2FT;
% 
% axes('Position', [250/400 240/600 100/400 100/600]); % 
% boxplot(matrixEarlyAmpSlopeG1G2,'orientation','vertical','whisker',0.5,'notch','on',...
%     'labels',{'Tal FT','Tal G2','Vin FT','Vin G2','Pax FT','Pax G2'},'symbol','','widths',boxWidth,'jitter',1,'colors','k', 'labelorientation','inline')
% txt = findobj(gca,'Type','text');
% set(txt(3:end),'VerticalAlignment', 'Middle');
% set(findobj(gca,'LineStyle','--'),'LineStyle','-')
% set(findobj(gca,'tag','Median'),'LineWidth',2)
% ylim([-50 1000])
% ylabel('Assembly rate (a.u./sec)')
% title({'Assembly rate'; 'for first 3 sec'})
% set(findobj(gca,'Type','Text'),'FontSize',7)
% set(gca,'FontSize',7)
% %% do some statistics:
% [pTalG1G2,hTalG1G2] = ranksum(earlyTalSlopeG1FT,earlyTalSlopeG2FT);
% [pVinG1G2,hAmpVinG1G2] = ranksum(earlyVinSlopeG1FT,earlyVinSlopeG2FT);
% [pPaxG1G2,hPaxG1G2] = ranksum(earlyPaxSlopeG1FT,earlyPaxSlopeG2FT);
% %put it in the plot
% text(1,200,['p=' num2str(pTalG1G2,2)]) %pax vs vin
% text(3,970,['p=' num2str(pVinG1G2,2)]) %tal vs vin
% text(5,300,['p=' num2str(pPaxG1G2,2)]) %pax vs tal
% set(findobj(gca,'Type','Text'),'FontSize',7)
%% Now dF/dt on initial vinculin binding
% % idxFAcontainingVinG1FT=arrayfun(@(x) any(strcmp(x.state(x.endingFrame-2:x.endingFrame),'FC') | strcmp(x.state(x.endingFrame-2:x.endingFrame),'FA')),tracksNAvinG1FT);
% % idxFAcontainingVinG1FT=arrayfun(@(x) any(strcmp(x.state,'FC') | strcmp(x.state,'FA')),tracksNAvinG1FT);
% % tracksNAvinG1FT2 =tracksNAvinG1FT(~idxFAcontainingVinG1FT);
% 
% [~,earlyForceSlopeG1FTvin] = arrayfun(@(x) regression(1:3,x.forceMag(x.startingFrameExtra:x.startingFrameExtra+2)),tracksNAvinG1FT);
% [~,earlyForceSlopeG1FTvin4] = arrayfun(@(x) regression(1:3,x.forceMag(x.startingFrameExtra:x.startingFrameExtra+2)),tracksNAvin4G1FT);
% [~,earlyForceSlopeG2FTvin] = arrayfun(@(x) regression(1:3,x.forceMag(x.startingFrameExtra:x.startingFrameExtra+2)),tracksNAvinG2FT2);
% [~,earlyForceSlopeG2FTvin4] = arrayfun(@(x) regression(1:3,x.forceMag(x.startingFrameExtra:x.startingFrameExtra+2)),tracksNAvin4G2FT2);
% %% plotting on dF/dt
% earlyForceSlopeG1FTvinTot=[earlyForceSlopeG1FTvin; earlyForceSlopeG1FTvin4];
% earlyForceSlopeG2FTvinTot=[earlyForceSlopeG2FTvin; earlyForceSlopeG2FTvin4];
% [lengthLongest,iLongest]=max([length(earlyForceSlopeG1FTvinTot),length(earlyForceSlopeG2FTvinTot)]);
% % it was vinculin that was longest
% matrixEarlyForceSlopeG1G2vin = NaN(lengthLongest,2);
% matrixEarlyForceSlopeG1G2vin(1:length(earlyForceSlopeG1FTvinTot),1) = earlyForceSlopeG1FTvinTot;
% matrixEarlyForceSlopeG1G2vin(1:length(earlyForceSlopeG2FTvinTot),2) = earlyForceSlopeG2FTvinTot;
% axes('Position', [250/400 50/600 100/400 100/600]); % 
% boxplot(matrixEarlyForceSlopeG1G2vin,'orientation','vertical','whisker',0.5,'notch','on',...
%     'labels',{'Vin FT','Vin G2'},'symbol','','widths',boxWidth,'jitter',1,'colors','k', 'labelorientation','inline')
% txt = findobj(gca,'Type','text');
% set(txt(3:end),'VerticalAlignment', 'Middle');
% set(findobj(gca,'LineStyle','--'),'LineStyle','-')
% set(findobj(gca,'tag','Median'),'LineWidth',2)
% ylim([-50 200])
% ylabel('Force growth rate (Pa/sec)')
% title({'Force growth rate during'; 'initial vinculin recruitment'})
% set(findobj(gca,'Type','Text'),'FontSize',7)
% set(gca,'FontSize',7)
% %% do some statistics:
% [pForceVinG1G2,hForceVinG1G2] = ranksum(earlyForceSlopeG1FTvinTot,earlyForceSlopeG2FTvinTot);
% %put it in the plot
% text(1,50,['p=' num2str(pForceVinG1G2,2)]) %pax vs vin
% set(findobj(gca,'Type','Text'),'FontSize',7)

%% lifetime checking
% figure, 
% histogram(lifeTimeVinG1FT,'Normalization','probability');
% hold on
% % [nVin] = histcounts(totalTLagVin,edges);
% histogram(lifeTimeVinG2FT,'Normalization','probability');
%% Supplement - individual examples
% talAxes=axes('Position', [1/4 0 1/4 1]);
% importfig(['/project/cellbiology/gdanuser/adhesion/Sangyoon/Alexia/31-01-2015/Vinc2/Colocalization/analysis1/eps/representativeMedian/track866early force peak.fig'],talAxes);
% vinAxes=axes('Position', [2/4 0 1/4 1]);
% importfig(['/project/cellbiology/gdanuser/adhesion/Sangyoon/Alexia/31-01-2015/Vinc2/Colocalization/analysis1/eps/representativeMedian/track866early force peak.fig'],vinAxes);
% paxAxes=axes('Position', [3/4 0 1/4 1]);
% importfig(['/project/cellbiology/gdanuser/adhesion/Sangyoon/Alexia/2015-07-24/Paxillin2/Colocalization/analysis2/eps/representativeMedian/track251veryGoodstartingFromZero.fig'],paxAxes);
%% save
print('-depsc', '-loose', [f3path filesep 'Fig3NTG1_G2.eps'])
print('-dtiff', '-loose', '-r300', [f3path filesep 'Fig3NTG1_G2.tif'])
savefig([f3path filesep 'Fig3NTG1_G2.fig'])

%% histogram - supplements
[nPax,edges] = histcounts(totalInitTLagPax2,60);
figure, 
h=histogram(totalInitTLagPax2,edges,'Normalization','probability');
hold on
% [nVin] = histcounts(totalTLagVin,edges);
histogram(totalInitTLagVin2,edges,'Normalization','probability');
set(gcf,'Position',[100 100 400 300])
legend('pax', 'vin')
xlabel('Time lag (sec)')
ylabel('Probability')
title(['Initial time lag for paxillin vs. vinculin'])
print('-depsc2', '-r150',[pathFigs filesep 'intensityDistribution' num2str(jj) 'thCellFrom' num2str(nStart) 'to' num2str(nEnd) 'thFrame.eps'])
% export_fig('eps',[pathFigs filesep 'intensityDistribution' num2str(jj) 'thCellFrom' num2str(nStart) 'to' num2str(nEnd) 'thFrame'])
% savefig([pathFigs filesep 'intensityDistribution' num2str(jj) 'thCellFrom' num2str(nStart) 'to' num2str(nEnd) 'thFrame.fig'])
%% histogram - supplements for peak
[nPax,edgesPeak] = histcounts(totalPeakTLagTal);
figure, 
h=histogram(totalPeakTLagTal,edgesPeak,'Normalization','probability');
hold on
% [nVin] = histcounts(totalTLagVin,edges);
histogram(totalPeakTLagVin,edgesPeak,'Normalization','probability');
set(gcf,'Position',[100 100 400 300])
legend('tal', 'vin')
xlabel('Time lag (sec)')
ylabel('Probability')
title(['Peak time lag for vinculin vs. talin'])
print('-depsc', '-loose', [f3_Tmeasure_path filesep 'histogramPeakLagVinVsTal.eps'])
print('-dtiff', '-r300', [f3_Tmeasure_path filesep 'histogramPeakLagVinVsTal'])
prevDir = pwd;
cd(f3_Tmeasure_path)
export_fig([f3_Tmeasure_path filesep 'histogramPeakLagVinVsTal.eps']);% histogramPeakLagVinVsTal -transparent
cd(prevDir)
