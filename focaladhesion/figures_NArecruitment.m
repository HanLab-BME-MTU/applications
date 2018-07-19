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
talFolder{2} = '/storage/disk2/Kevin/2017-02/Talin1_ChoK1_5kPa_1520_01/Colocalization/analysis2';
talFolder{3} = '/storage/disk2/Kevin/2017-02/Talin1_ChoK1_5kPa_1520_03/Colocalization/analysis2';
talFolder{4} = '/storage/disk2/Kevin/2017-02/Talin1_ChoK1_5kPa_1520_05_2017_02_09/Colocalization/analysis2';
talFolder{5} = '/storage/disk2/Kevin/2017-02/Talin1_ChoK1_5kPa_1520_06/Colocalization/analysis2';
talFolder{6} = '/storage/disk2/Kevin/2017-02/Talin1_ChoK1_5kPa_1520_06_2017_02_09/Colocalization/analysis2';
talFolder{7} = '/storage/disk2/Kevin/2017-02/Talin1_ChoK1_5kPa_1520_07/Colocalization/analysis2';
talFolder{8} = '/storage/disk2/Kevin/2017-02/Talin1_ChoK1_5kPa_1520_07_2017_02_09/Colocalization/analysis2';
talFolder{9} = '/storage/disk2/Kevin/2017-02/Talin1_ChoK1_5kPa_1520_09_2017_02_09/Colocalization/analysis2';
talFolder{10} = '/storage/disk2/Kevin/2017-02/Talin1_ChoK1_5kPa_1520_08_2017_02_09/Colocalization/analysis2';
% talFolder{2} = '/project/bioinformatics/Danuser_lab/P01adhesion/analysis/Sangyoon/NA_RecruitmentProject/Kevin/2015-07-22_16kPa/Cell2/Colocalization/analysis3';
% talFolder{3} = '/project/bioinformatics/Danuser_lab/P01adhesion/analysis/Sangyoon/NA_RecruitmentProject/Kevin/2015-07-22_16kPa/Cell4/Colocalization/analysis3';
% talFolder{4} = '/project/bioinformatics/Danuser_lab/P01adhesion/analysis/Sangyoon/NA_RecruitmentProject/Kevin/2015-07-22_16kPa/Cell5/Colocalization/analysis3';
% talFolder{5} = '/project/bioinformatics/Danuser_lab/P01adhesion/analysis/Sangyoon/NA_RecruitmentProject/Kevin/2015-07-22_16kPa/Cell8/Colocalization/analysis3';
% talFolder{6} = '/project/bioinformatics/Danuser_lab/P01adhesion/analysis/Sangyoon/NA_RecruitmentProject/Kevin/2015-07-22_16kPa/Cell10/Colocalization/analysis2';
% talFolder{7} = '/project/bioinformatics/Danuser_lab/P01adhesion/analysis/Sangyoon/NA_RecruitmentProject/Kevin/2015-07-22_16kPa/Cell6/Colocalization/analysis3';
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
%% Fig 1 High-throughput collection of NA tracks with classification
% Showing overall look over each protein
iRepTal=1; iRepVin=2; iRepPax=4;
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
%% ------------------------ Talin--------data loading------------------------------------------
talForceStack = load([talFolder{iRepTal} filesep 'fMap' filesep 'tMap.mat'],'tMap');
talForceStack = talForceStack.tMap;
talImgStack = load([talFolder{iRepTal} filesep 'pax' filesep 'paxImgStack.mat'],'paxImgStack');
talImgStack = talImgStack.paxImgStack;
forceFieldTal = load([MDpathTal filesep 'TFMPackage' filesep 'forceField' filesep 'forceField.mat']);
forceFieldTal = forceFieldTal.forceField;
displFieldTal = load([MDpathTal filesep 'TFMPackage' filesep 'correctedDisplacementField' filesep 'displField.mat']);
displFieldTal = displFieldTal.displField;
%% ------------------------ Talin---------showing-----------------------------------------
widthPixelTalin = unifiedWidthNanometer/MDtal.pixelSize_;
avgWidth=1;
CurrentFrameTal=300;
CurrentFrameTal = min(size(talImgStack,3)-avgWidth,CurrentFrameTal-avgWidth);
pos(3:4) = [600 600];
figure('Position', pos, 'PaperPositionMode', 'auto','Color','w');
% subplot(2,1,1)
% force
axes('Position', [30/pos(3) 420/pos(4) 160/pos(3) 160/pos(4)]); % talin
fMaxTal=300;
imshow(mean(talForceStack(:,:,CurrentFrameTal-avgWidth:CurrentFrameTal+avgWidth),3),[0 fMaxTal])
load('/project/bioinformatics/Danuser_lab/P01adhesion/analysis/Sangyoon/NA_RecruitmentProject/MyColormaps','mycmap')
colormap(gca,mycmap); 
hold on
overlayTractionVectors(forceFieldTal(CurrentFrameTal),displFieldTal(CurrentFrameTal),0.1,3,'Color',[0.7 0.7 0.7],'ShiftField',false);
hc=colorbar('east'); hc.Position(3)=0.01; hc.Ticks=[0 fMaxTal];
set(gca,'XLim',[30 30+widthPixelTalin],'YLim',[30 30+widthPixelTalin])
hold on
% line([30+20 30+20+round(5000/MDtal.pixelSize_)],[30+widthPixelTalin-20 30+widthPixelTalin-20],'LineWidth',2,'Color','w')
title('Talin')

axes('Position', [30/pos(3) 240/pos(4) 160/pos(3) 160/pos(4)]); % talin
% axes('Position', [0.03 .35 .3 .30]); % talin
% imshow(imcomplement(mean(talImgStack(:,:,195:200),3)),[-900 -100])
imshow(imcomplement(mean(talImgStack(:,:,CurrentFrameTal-avgWidth:CurrentFrameTal+avgWidth),3)),[-800 -100])
set(gca,'XLim',[30 30+widthPixelTalin],'YLim',[30 30+widthPixelTalin])
% We want to keep the same zoom scale for all three images
% set(gca,'XLim',[64 64+widthPixelTalin],'YLim',[120 120+widthPixelTalin])
hold on
line([30+20 30+20+round(5000/MDtal.pixelSize_)],[30+widthPixelTalin-20 30+widthPixelTalin-20],'LineWidth',2,'Color','k')
% % line([64+10 64+10+round(5000/MDtal.pixelSize_)],[120+widthPixelTalin-20 120+widthPixelTalin-20],'LineWidth',2,'Color','k')
text(30+20, 30+widthPixelTalin-20-30,'5 um','Color','k','Fontsize',7)
%% vinculin loading
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
CurrentFrame = 250;
% force
widthPixelVinculin = unifiedWidthNanometer/MDvin5.pixelSize_;
axes('Position', [(160+30+20)/pos(3) 420/pos(4) 160/pos(3) 160/pos(4)]); % talin
imshow(mean(vinForceStack(:,:,CurrentFrame-2:CurrentFrame+2),3),[0 1800])
colormap(gca,mycmap); 
hold on
overlayTractionVectors(forceFieldVin(CurrentFrame),displFieldVin(CurrentFrame),0.04,3,'Color',[0.7 0.7 0.7],'ShiftField',false);
hc=colorbar('east'); hc.Position(3)=0.01; hc.Ticks=[0 1800];
set(gca,'XLim',[21 21+widthPixelVinculin],'YLim',[32 32+widthPixelVinculin])
title('Vinculin')
roiXmin = 140;
roiXmax = 310;
roiYmin = 60;
roiYmax = 230;
minWidth = min([roiXmax-roiXmin,roiYmax-roiYmin]);
hold on
rectangle('Position',[roiXmin,roiYmin,minWidth,minWidth],'EdgeColor','r')

axes('Position', [(160+30+20)/pos(3) 240/pos(4) 160/pos(3) 160/pos(4)]); % talin
imshow(imcomplement(mean(vinImgStack(:,:,CurrentFrame-2:CurrentFrame+2),3)),[])
% imshow(imcomplement(mean(vinImgStack(:,:,CurrentFrame-2:CurrentFrame+2),3)),[-700 -100])
set(gca,'XLim',[21 21+widthPixelVinculin],'YLim',[32 32+widthPixelVinculin])
% set(gca,'XLim',[30 30+widthPixelTalin],'YLim',[30 30+widthPixelTalin])
hold on
% line([23+10 23+10+round(5000/MDvin.pixelSize_)],[38+15 38+15],'LineWidth',2,'Color','w')
% text(45+10, 16+15+15,'5 um','Color','w','Fontsize',7)
rectangle('Position',[roiXmin,roiYmin,minWidth,minWidth],'EdgeColor','r')
%% pax loading
paxForceStack = load([paxFolder{iRepPax} filesep 'fMap' filesep 'tMap.mat'],'tMap');
paxForceStack = paxForceStack.tMap;
paxImgStack = load([paxFolder{iRepPax} filesep 'pax' filesep 'paxImgStack.mat'],'paxImgStack');
paxImgStack = paxImgStack.paxImgStack;
forceFieldPax = load([MDpathPax filesep 'TFMPackage' filesep 'forceField' filesep 'forceField.mat']);
forceFieldPax = forceFieldPax.forceField;
displFieldPax = load([MDpathPax filesep 'TFMPackage' filesep 'correctedDisplacementField' filesep 'displField.mat']);
displFieldPax = displFieldPax.displField;
%% ------------------------ paxillin--------------------------------------------------
% axes('Position', [0.67 .65 .3 .35]); % talin
% subplot(2,1,1)
% force
CurrentFramePax=298;
widthPixelPaxilln = unifiedWidthNanometer/MDpax.pixelSize_;
axes('Position', [(160+30+20+160+20)/pos(3) 420/pos(4) 160/pos(3) 160/pos(4)]); % talin
imshow(mean(paxForceStack(:,:,295:300),3),[0 1000])
colormap(gca,mycmap); 
hold on
overlayTractionVectors(forceFieldPax(CurrentFramePax),displFieldPax(CurrentFramePax),0.04,3,'Color',[0.7 0.7 0.7],'ShiftField',false);

hc=colorbar('east'); hc.Position(3)=0.01; hc.Ticks=[0 1000];
set(gca,'XLim',[30 30+widthPixelPaxilln],'YLim',[30 30+widthPixelPaxilln])
title('Paxillin')

axes('Position', [(160+30+20+160+20)/pos(3) 240/pos(4) 160/pos(3) 160/pos(4)]); % talin
imshow(imcomplement(mean(paxImgStack(:,:,295:300),3)),[-1500 -100])
set(gca,'XLim',[30 30+widthPixelPaxilln],'YLim',[30 30+widthPixelPaxilln])
hold on
% line([23+10 23+10+round(5000/MDvin.pixelSize_)],[38+15 38+15],'LineWidth',2,'Color','w')
% text(45+10, 16+15+15,'5 um','Color','w','Fontsize',7)
%% Fig 1 d. vinculin + tracks
axes('Position', [30/pos(3) 40/pos(4) 160/pos(3) 160/pos(4)]); 
% axes('Position', [0.03 .4 .3 .25]); % talin
imshow(imcomplement(mean(vinImgStack(:,:,CurrentFrame-2:CurrentFrame+2),3)),[-29000 -3000])
set(gca,'XLim',[roiXmin roiXmin+minWidth],'YLim',[roiYmin roiYmin+minWidth])
hold on
% line([23+10 23+10+round(5000/MDvin.pixelSize_)],[38+15 38+15],'LineWidth',2,'Color','w')
% text(45+10, 16+15+15,'5 um','Color','w','Fontsize',7)
if ~exist('tracksNAvin','var')
    tic
    tracksNAvin = load([vinFolder{iRepVin} filesep 'data' filesep 'tracksNA.mat'],'tracksNA');
    tracksNAvin = tracksNAvin.tracksNA; toc
end
idCurrent = arrayfun(@(x) x.startingFrameExtra<=CurrentFrame & x.endingFrameExtra>=CurrentFrame,tracksNAvin);
% pstructVin = pointSourceDetection(mean(vinImgStack(:,:,195:200),3),1.8);
% plot(pstructVin.x,pstructVin.y,'bo')
markerSize=4;
plot(arrayfun(@(x) x.xCoord(CurrentFrame),tracksNAvin(idCurrent)),arrayfun(@(x) x.yCoord(CurrentFrame),tracksNAvin(idCurrent)),'ro','MarkerSize',markerSize)
xmat = cell2mat(arrayfun(@(x) x.xCoord(1:CurrentFrame),tracksNAvin(idCurrent),'UniformOutput',false));
ymat = cell2mat(arrayfun(@(x) x.yCoord(1:CurrentFrame),tracksNAvin(idCurrent),'UniformOutput',false));
plot(xmat',ymat','r','linewidth',0.25)

idAdhLogic = arrayfun(@(x) ~isempty(x.adhBoundary),tracksNAvin);
idAdhCur = arrayfun(@(x) ~isempty(x.adhBoundary{CurrentFrame}),tracksNAvin(idAdhLogic));
idAdh = find(idAdhLogic);
idAdhCur = idAdh(idAdhCur);
arrayfun(@(x) plot(x.adhBoundary{CurrentFrame}(:,1),x.adhBoundary{CurrentFrame}(:,2), 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5),tracksNAvin(idAdhCur))

%% Cell boundary
iSDCProc =MDvin.getProcessIndex('StageDriftCorrectionProcess',1,1);     
SDCProc=MDvin.processes_{iSDCProc};
iBeadChan=1;
s = load(SDCProc.outFilePaths_{3,iBeadChan},'T');    
T = s.T;
maxX = ceil(max(abs(T(:, 2))));
maxY = ceil(max(abs(T(:, 1))));
Tr = maketform('affine', [1 0 0; 0 1 0; fliplr(T(CurrentFrame, :)) 1]);
iChan=2;
iMask = MDvin.getProcessIndex('MaskRefinementProcess');
maskProc = MDvin.getProcess(iMask);
bwPI4 = maskProc.loadChannelOutput(iChan,CurrentFrame);
Ibw = padarray(bwPI4, [maxY, maxX]);
bwPI4 = imtransform(Ibw, Tr, 'XData',[1 size(Ibw, 2)],'YData', [1 size(Ibw, 1)]);
[B,~,nBD]  = bwboundaries(bwPI4,'noholes');
for kk=1:nBD
    boundary = B{kk};
    plot(boundary(:,2), boundary(:,1), 'b', 'LineWidth', 0.5) % cell boundary
end
title('Adhesion tracking')

line([roiXmin+5 roiXmin+5+round(2000/MDtal.pixelSize_)],[roiYmin+5 roiYmin+5],'LineWidth',2,'Color','k')
text(roiXmin+5, roiYmin+5+7,'2 um','Color','k','Fontsize',7)
%% Feature selection
% axes('Position', [0.35 .4 .3 .25]); % talin
% axis off
% title('Features')

%% Classification result
axes('Position', [(160+30+20+160+20)/pos(3) 40/pos(4) 160/pos(3) 160/pos(4)]); % talin
% axes('Position', [0.67 .4 .3 .25]); % talin
imshow(imcomplement(mean(vinImgStack(:,:,195:200),3)),[-29000 -3000])
set(gca,'XLim',[roiXmin roiXmin+minWidth],'YLim',[roiYmin roiYmin+minWidth])

% Drawing
if ~exist('allDataClass','var')
    try
        trainedClassifier =  load([vinFolder{iRepVin} filesep 'data'  filesep 'trainedClassifier_fromSelectedGroups.mat'],'trainedClassifier');
        trainedClassifier = trainedClassifier.trainedClassifier;
        % trainedClassifier = trainClassifierNA(T);
        [~,allData] = extractFeatureNA(tracksNAvin);
        allDataClass = predict(trainedClassifier,allData);
    catch
        try
            trainedClassifier =  load([vinFolder{iRepVin} filesep 'data'  filesep 'trainedClassifier.mat'],'trainedClassifier');
            trainedClassifier = trainedClassifier.trainedClassifier;
            % trainedClassifier = trainClassifierNA(T);
            [~,allData] = extractFeatureNA(tracksNAvin);
            allDataClass = predict(trainedClassifier,allData);
        catch
            try
                idG = load([vinFolder{iRepVin} filesep 'data' filesep 'idsClassified.mat']);
                idGroupSelected={idG.idGroup1filtered,idG.idGroup2filtered,idG.idGroup3filtered,idG.idGroup4filtered,idG.idGroup5filtered,idG.idGroup6filtered,....
                                                idG.idGroup7filtered,idG.idGroup8filtered,idG.idGroup9filtered};
                [T,allData]=extractFeatureNA(tracksNAvin,idGroupSelected,2,MDvin);
                [trainedClassifier, validationAccuracy, C, order] = trainClassifierNA(T);
                allDataClass = predict(trainedClassifier,allData);
                save([vinFolder{iRepVin} filesep 'data'  filesep 'trainedClassifierFrom_idsClassified.mat'],'trainedClassifier');
            catch
                idG = load([vinFolder{iRepVin} filesep 'data' filesep 'selectedGroups.mat']);
                idG2 = load([vinFolder{iRepVin} filesep 'data' filesep 'selectedGroupsNew.mat']);
                idG3 = load([vinFolder{iRepVin} filesep 'data' filesep 'selectedGroupsOld.mat']);
                idGroupSelected={[idG.idGroup1Selected idG2.idGroup1Selected idG3.idGroup1Selected],...
                    [idG.idGroup2Selected idG2.idGroup2Selected idG3.idGroup2Selected],...
                    [idG.idGroup3Selected idG2.idGroup3Selected idG3.idGroup3Selected],...
                    [idG.idGroup4Selected idG2.idGroup4Selected idG3.idGroup4Selected],...
                    [idG.idGroup5Selected idG2.idGroup5Selected idG3.idGroup5Selected],...
                    [idG.idGroup6Selected idG2.idGroup6Selected idG3.idGroup6Selected],...
                    [idG.idGroup7Selected idG2.idGroup7Selected idG3.idGroup7Selected],...
                    [idG.idGroup8Selected idG2.idGroup8Selected idG3.idGroup8Selected],...
                    [idG.idGroup9Selected idG2.idGroup9Selected idG3.idGroup9Selected]};
                [T,allData]=extractFeatureNA(tracksNAvin,idGroupSelected,2,MDvin);
                [trainedClassifier, validationAccuracy, C, order] = trainClassifierNA(T);
                allDataClass = predict(trainedClassifier,allData);
                save([vinFolder{iRepVin} filesep 'data'  filesep 'trainedClassifier_fromSelectedGroups.mat'],'trainedClassifier');
                idGroup1 = strcmp(allDataClass,'Group1');
                idGroup2 = strcmp(allDataClass,'Group2');
                idGroup3 = strcmp(allDataClass,'Group3');
                idGroup4 = strcmp(allDataClass,'Group4');
                idGroup5 = strcmp(allDataClass,'Group5');
                idGroup6 = strcmp(allDataClass,'Group6');
                idGroup7 = strcmp(allDataClass,'Group7');
                idGroup8 = strcmp(allDataClass,'Group8');
                idGroup9 = strcmp(allDataClass,'Group9');
                save([vinFolder{iRepVin} filesep 'data'  filesep 'idsClassified_FromOwnLabels.mat'],'idGroup1','idGroup2','idGroup3',...
                    'idGroup4','idGroup5','idGroup6','idGroup7','idGroup8','idGroup9','-v7.3')
            end
        end
    end
end
hold on
drawClassifiedTracks(allDataClass,tracksNAvin,CurrentFrame,gca,10);
title('Classified tracks')
%% saving
print('-depsc', '-loose', [f1path filesep 'Fig1Methods2.eps'])
savefig([f1path filesep 'Fig1Methods2.fig'])
print('-dtiff', '-loose', '-r300', [f1path filesep 'Fig1Methods2.tif'])
close gcf
%% Overall spatial correlation - talin
% talin image - accumulating
accumulTalIntensity=cell(numel(talFolder),1);
accumulTalForce=cell(numel(talFolder),1);
CurrentFrameTal=100;
iChan=2;
for ii=1:numel(talFolder)
    curMDpathTal = fileparts(fileparts(talFolder{ii}));
    curMD = MovieData.load([curMDpathTal filesep 'movieData.mat']);
    SDCProc=curMD.processes_{curMD.getProcessIndex('StageDriftCorrectionProcess',1,1)};
    s = load(SDCProc.outFilePaths_{3,iBeadChan},'T');    
    T = s.T;
    maxX = ceil(max(abs(T(:, 2))));
    maxY = ceil(max(abs(T(:, 1))));
    curNumFrames=curMD.nFrames_;

    try
        talForceStack = load([talFolder{ii} filesep 'fMap' filesep 'tMap.mat'],'tMap');
        talForceStack = talForceStack.tMap;
        talImgStack = load([talFolder{ii} filesep 'pax' filesep 'paxImgStack.mat'],'paxImgStack');
        talImgStack = talImgStack.paxImgStack;
    catch 
        continue;
    end
    
    for CurrentFrameTal=1:curNumFrames
        Tr = maketform('affine', [1 0 0; 0 1 0; fliplr(T(CurrentFrameTal, :)) 1]);
        
%         iMask = curMD.getProcessIndex('MaskRefinementProcess');
%         maskProc = curMD.getProcess(iMask);
%         bwPI4 = maskProc.loadChannelOutput(iChan,CurrentFrameTal);
        
        FASegPackage = curMD.getPackage(curMD.getPackageIndex('FocalAdhesionSegmentationPackage'));
        FASegProc = FASegPackage.processes_{6};
        bwPI4 = FASegProc.loadChannelOutput(iChan,CurrentFramePax);
        bwPI4 = bwPI4>0;


        Ibw = padarray(bwPI4, [maxY, maxX]);
        bwPI4 = imtransform(Ibw, Tr, 'XData',[1 size(Ibw, 2)],'YData', [1 size(Ibw, 1)]);
        
        try
            curTalImg=talImgStack(:,:,CurrentFrameTal);
            curTalForce=talForceStack(:,:,CurrentFrameTal);

            accumulTalIntensity{ii}{CurrentFrameTal} = curTalImg(bwPI4);
            accumulTalForce{ii}{CurrentFrameTal} = curTalForce(bwPI4);
            disp([num2str(ii) 'th cell, ' num2str(CurrentFrameTal) 'th frame done.'])
        catch
            continue
        end
    end
end
%% drawing
    % figure, plot(curTalImg(:),curTalForce(:),'.')
figure('Position', [100 100 900 200]), 
subplot(1,4,1);
xmin=50;xmax=450; ymax=4000;
iCurTal=1; iFrame=60;
curTalIntenAllFrames=cell2mat(accumulTalIntensity{iCurTal}');
curTalForceAllFrames=cell2mat(accumulTalForce{iCurTal}');
idxNZF=curTalForceAllFrames>0;
% densityplot(accumulVinIntensity{iCurVin},accumulVinForce{iCurVin},linspace(0, xmax, 100),linspace(0, ymax, 100),'DisplayFunction',@log); 
% densityplot(accumulTalIntensity{iCurTal}{iFrame},accumulTalForce{iCurTal}{iFrame},linspace(xmin, xmax, 100),linspace(0, ymax, 100),'DisplayFunction',@log); 
% densityplot(log(curTalIntenAllFrames),log(curTalForceAllFrames),linspace(xmin, xmax, 100),linspace(0, ymax, 100),'DisplayFunction',@log); 
densityplot(log(1+curTalIntenAllFrames(idxNZF)),log(1+curTalForceAllFrames(idxNZF)),'DisplayFunction',@log); 
% densityplot(accumulTalIntensity,accumulTalForce,linspace(0, xmax, 100),linspace(0, ymax, 100),'DisplayFunction',@log); 
% ylim([5 ymax-5]);
colormap(gca,mycmap); 
title(['Talin, r=' num2str(corr(log(1+curTalIntenAllFrames(idxNZF)),log(1+curTalForceAllFrames(idxNZF))))])

%% Overall spatial correlation - vinculin
accumulVinIntensity=cell(numel(vinFolder),1);
accumulVinForce=cell(numel(vinFolder),1);
CurrentFrameVin = 150;
for ii=1:numel(vinFolder)
    curMDpathVin = fileparts(fileparts(vinFolder{ii}));
    curMD = MovieData.load([curMDpathVin filesep 'movieData.mat']);
    SDCProc=curMD.processes_{curMD.getProcessIndex('StageDriftCorrectionProcess',1,1)};
    s = load(SDCProc.outFilePaths_{3,iBeadChan},'T');    
    T = s.T;
    maxX = ceil(max(abs(T(:, 2))));
    maxY = ceil(max(abs(T(:, 1))));
    curNumFrames=curMD.nFrames_;
    try
        vinForceStack = load([vinFolder{ii} filesep 'fMap' filesep 'tMap.mat'],'tMap');
        vinForceStack = vinForceStack.tMap;
        vinImgStack = load([vinFolder{ii} filesep 'pax' filesep 'paxImgStack.mat'],'paxImgStack');
        vinImgStack = vinImgStack.paxImgStack;
    catch
        continue
    end

    for CurrentFrameVin=1:curNumFrames
        Tr = maketform('affine', [1 0 0; 0 1 0; fliplr(T(CurrentFrameVin, :)) 1]);

%         iMask = curMD.getProcessIndex('MaskRefinementProcess');
%         maskProc = curMD.getProcess(iMask);
%         bwPI4 = maskProc.loadChannelOutput(iChan,CurrentFrameVin);

        FASegPackage = curMD.getPackage(curMD.getPackageIndex('FocalAdhesionSegmentationPackage'));
        FASegProc = FASegPackage.processes_{6};
        bwPI4 = FASegProc.loadChannelOutput(iChan,CurrentFramePax);
        bwPI4 = bwPI4>0;

        Ibw = padarray(bwPI4, [maxY, maxX]);
        bwPI4 = imtransform(Ibw, Tr, 'XData',[1 size(Ibw, 2)],'YData', [1 size(Ibw, 1)]);
        try
            curVinImg=vinImgStack(:,:,CurrentFrameVin);
            curVinForce=vinForceStack(:,:,CurrentFrameVin);

            accumulVinIntensity{ii}{CurrentFrameVin} = curVinImg(bwPI4);
            accumulVinForce{ii}{CurrentFrameVin} = curVinForce(bwPI4);
            disp([num2str(ii) 'th cell, ' num2str(CurrentFrameVin) 'th frame done.'])
        catch
            continue
        end
    end
end
%% drawing
    % figure, plot(curVinImg(:),curVinForce(:),'.')
% figure, plot(curVinImg(bwPI4),curVinForce(bwPI4),'.')
subplot(1,4,2)
xmax=60005; ymax=4000;
iCurVin=5;
curVinIntenAllFrames=cell2mat(accumulVinIntensity{iCurVin}');
curVinForceAllFrames=cell2mat(accumulVinForce{iCurVin}');
idxNZF=curVinForceAllFrames>0;
% densityplot(log(curVinIntenAllFrames),log(curVinForceAllFrames),linspace(xmin, xmax, 100),linspace(0, ymax, 100),'DisplayFunction',@log); 
densityplot(log(1+curVinIntenAllFrames(idxNZF)),log(1+curVinForceAllFrames(idxNZF)),'DisplayFunction',@log); 
% densityplot(accumulVinIntensity{iCurVin},accumulVinForce{iCurVin},linspace(0, xmax, 100),linspace(0, ymax, 100),'DisplayFunction',@log); 
% densityplot(accumulVinIntensity{iCurVin},accumulVinForce{iCurVin},'DisplayFunction',@log); 
% xlim([0 xmax-1000]);
% ylim([5 ymax-5]);
colormap(gca,mycmap); 
% densityplot(curVinImg(bwPI4),curVinForce(bwPI4),'DisplayFunction',@log); colormap(gca,mycmap); 
% ylim([0 1800]);
% title(['Vinculin, r=' num2str(corr(curVinImg(bwPI4),curVinForce(bwPI4)))])
title(['Vinculin, r=' num2str(corr(log(1+curVinIntenAllFrames(idxNZF)),log(1+curVinForceAllFrames(idxNZF))))])

%% Overall spatial correlation - paxillin
accumulPaxIntensity=cell(numel(paxFolder),1);
accumulPaxForce=cell(numel(paxFolder),1);
CurrentFramePax=100;
for ii=1:numel(paxFolder)
    curMDpathPax = fileparts(fileparts(paxFolder{ii}));
    curMD = MovieData.load([curMDpathPax filesep 'movieData.mat']);
    SDCProc=curMD.processes_{curMD.getProcessIndex('StageDriftCorrectionProcess',1,1)};
    s = load(SDCProc.outFilePaths_{3,iBeadChan},'T');    
    T = s.T;
    maxX = ceil(max(abs(T(:, 2))));
    maxY = ceil(max(abs(T(:, 1))));
    curNumFrames=curMD.nFrames_;
    try
        paxForceStack = load([paxFolder{ii} filesep 'fMap' filesep 'tMap.mat'],'tMap');
        paxForceStack = paxForceStack.tMap;
        paxImgStack = load([paxFolder{ii} filesep 'pax' filesep 'paxImgStack.mat'],'paxImgStack');
        paxImgStack = paxImgStack.paxImgStack;
    catch
        continue
    end

    for CurrentFramePax=1:curNumFrames
        Tr = maketform('affine', [1 0 0; 0 1 0; fliplr(T(CurrentFramePax, :)) 1]);
        
        FASegPackage = curMD.getPackage(curMD.getPackageIndex('FocalAdhesionSegmentationPackage'));
        FASegProc = FASegPackage.processes_{6};
        bwPI4 = FASegProc.loadChannelOutput(iChan,CurrentFramePax);
        bwPI4 = bwPI4>0;

%         iMask = curMD.getProcessIndex('MaskRefinementProcess');
%         maskProc = curMD.getProcess(iMask);
%         bwPI4 = maskProc.loadChannelOutput(iChan,CurrentFramePax);
        
        Ibw = padarray(bwPI4, [maxY, maxX]);
        bwPI4 = imtransform(Ibw, Tr, 'XData',[1 size(Ibw, 2)],'YData', [1 size(Ibw, 1)]);

        try
            curPaxImg=paxImgStack(:,:,CurrentFramePax);
            curPaxForce=paxForceStack(:,:,CurrentFramePax);

            accumulPaxIntensity{ii}{CurrentFramePax} = curPaxImg(bwPI4);
            accumulPaxForce{ii}{CurrentFramePax} = curPaxForce(bwPI4);
            disp([num2str(ii) 'th cell, ' num2str(CurrentFramePax) 'th frame done.'])
        catch
            continue
        end
    end
end
%%
    % figure, plot(curPaxImg(:),curPaxForce(:),'.')
subplot(1,4,3);
xmax=60000; ymax=5000;
iCurPax=2;
% densityplot(accumulVinIntensity{iCurVin},accumulVinForce{iCurVin},linspace(0, xmax, 100),linspace(0, ymax, 100),'DisplayFunction',@log); 
curPaxIntenAllFrames=cell2mat(accumulPaxIntensity{iCurPax}(199:201)');
curPaxForceAllFrames=cell2mat(accumulPaxForce{iCurPax}(199:201)');
idxNZF=curPaxForceAllFrames>0;
% densityplot(curPaxIntenAllFrames,curPaxForceAllFrames,linspace(xmin, xmax, 100),linspace(0, ymax, 100),'DisplayFunction',@log); 
densityplot(log(1+curPaxIntenAllFrames(idxNZF)),log(1+curPaxForceAllFrames(idxNZF)),'DisplayFunction',@log); 
% densityplot(accumulPaxIntensity{iCurPax},accumulPaxForce{iCurPax},linspace(0, xmax, 100),linspace(0, ymax, 100),'DisplayFunction',@log); 
colormap(gca,mycmap); 
% figure, plot(curPaxImg(:),curPaxForce(:),'.')
% subplot(1,3,3), densityplot(curPaxImg(bwPI4),curPaxForce(bwPI4),'DisplayFunction',@log); colormap(gca,mycmap); 
% ylim([0 1000]);
title(['Paxillin, r=' num2str(corr(log(1+curPaxIntenAllFrames(idxNZF)),log(1+curPaxForceAllFrames(idxNZF))))])
%% correlation plot
idValidTal=~cellfun(@isempty,accumulTalIntensity);
idValidVin=~cellfun(@isempty,accumulVinIntensity);
idValidPax=~cellfun(@isempty,accumulPaxIntensity);
p=0;
for ii=find(idValidTal)'
    curMDpathTal = fileparts(fileparts(talFolder{ii}));
    curMD = load([curMDpathTal filesep 'movieData.mat']);
    curMD = curMD.MD;
    nFrames=curMD.nFrames_;
    
    for jj=1:nFrames
        p=p+1;
        idxNZF=accumulTalForce{ii}{jj}>0;
        rho{1}(p)=corr(log(1+accumulTalIntensity{ii}{jj}(idxNZF)),log(1+accumulTalForce{ii}{jj}(idxNZF)));
    end
end
p=0;
for ii=find(idValidVin)'
    curMDpathVin = fileparts(fileparts(vinFolder{ii}));
    curMD = load([curMDpathVin filesep 'movieData.mat']);
    curMD = curMD.MD;
    nFrames=curMD.nFrames_;
    for jj=1:nFrames
        p=p+1;
        idxNZF=accumulVinForce{ii}{jj}>0;
        rho{2}(p)=corr(log(1+accumulVinIntensity{ii}{jj}(idxNZF)),log(1+accumulVinForce{ii}{jj}(idxNZF)));
    end
end
p=0;
for ii=find(idValidPax)'
    curMDpathPax = fileparts(fileparts(paxFolder{ii}));
    curMD = load([curMDpathPax filesep 'movieData.mat']);
    curMD = curMD.MD;
    nFrames=curMD.nFrames_;
    for jj=1:nFrames
        p=p+1;
        idxNZF=accumulPaxForce{ii}{jj}>0;
        rho{3}(p)=corr(log(1+accumulPaxIntensity{ii}{jj}(idxNZF)),log(1+accumulPaxForce{ii}{jj}(idxNZF)));
    end
end
subplot(1,4,4)
boxPlotCellArray(rho,{'Tal','Vin','Pax'})
%% saving
print('-depsc', '-loose', [f1path filesep 'Fig1Middle_spatial_corr2.eps'])
savefig([f1path filesep 'Fig1Middle_spatial_corr2.fig'])
print('-dtiff', '-loose', '-r300', [f1path filesep 'Fig1Middle_spatial_corr2.tif'])
close gcf

%% Example of intensity, force, edge position and adhesion movement in g1
if ~exist('idGroupsVin5struct','var')
    idGroupsVin5struct = load([vinFolder{iRepVin} filesep 'data' filesep 'idsClassified.mat'],'idGroup1filtered','idGroup2filtered','idGroup3filtered',...
    'idGroup4filtered','idGroup5filtered','idGroup6filtered','idGroup7filtered','idGroup8filtered','idGroup9filtered',...
    'idGroup1','idGroup2','idGroup3',...
    'idGroup4','idGroup5','idGroup6','idGroup7','idGroup8','idGroup9');
%     idGroupsVin5struct = load([vinFolder{iRepVin} filesep 'data' filesep 'idsClassified_FromOwnLabels.mat'],...
%     'idGroup1','idGroup2','idGroup3',...
%     'idGroup4','idGroup5','idGroup6','idGroup7','idGroup8','idGroup9');
%     try
%     idGroupsVin5{1} = idGroupsVin5struct.idGroup1filtered;
%     idGroupsVin5{2} = idGroupsVin5struct.idGroup2filtered;
%     idGroupsVin5{3} = idGroupsVin5struct.idGroup3filtered;
%     idGroupsVin5{4} = idGroupsVin5struct.idGroup4filtered;
%     idGroupsVin5{5} = idGroupsVin5struct.idGroup5filtered;
%     idGroupsVin5{6} = idGroupsVin5struct.idGroup6filtered;
%     idGroupsVin5{7} = idGroupsVin5struct.idGroup7filtered;
%     idGroupsVin5{8} = idGroupsVin5struct.idGroup8filtered;
%     idGroupsVin5{9} = idGroupsVin5struct.idGroup9filtered;
%     catch
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
%% Intensity
fS2path = '/home2/shan/ownCloud/documents/NascentAdhesionForce/Manuscript/Figures/FigS2NonMaturingVsMaturing';
pos(3:4) = [650 750];
close
figure('Position', pos, 'PaperPositionMode', 'auto','Color','w');
idGroup1 = find(idGroupsVin5{1});
curID=157;
% curID=curID+1;
% curID=157,170,351,370, 417, 424, 431; %good candidate
curTrack = tracksNAvin(idGroup1(curID));
axes('Position', [0.07 .85 .12 .11]);
plotIntensityForceSingle(curTrack,idGroup1(curID),1,0.01,0.3,MDvin.timeInterval_);
ylim([0 1000])
% force
axes('Position', [0.27 .85 .12 .11]);
plotIntensityForceSingle(curTrack,idGroup1(curID),2,0.01,0.3,MDvin.timeInterval_);
ylim([0 130])
% edge position
axes('Position', [0.47 .85 .12 .11]);
plotIntensityForceSingle(curTrack,idGroup1(curID),3,0.01,0.3,MDvin.timeInterval_);
% edge position
axes('Position', [0.67 .85 .12 .11]);
plotIntensityForceSingle(curTrack,idGroup1(curID),4,0.01,0.3,MDvin.timeInterval_);
% distToEdge
axes('Position', [0.87 .85 .12 .11]);
plotIntensityForceSingle(curTrack,idGroup1(curID),5,0.01,0.3,MDvin.timeInterval_);
ylim([0 50])

% G2 example
idGroup2 = find(idGroupsVin5{2});
curID2=14;
% curID2=curID2+1;
% curID2=4,5,7,14; %good candidate
curTrack = tracksNAvin(idGroup2(curID2));
axes('Position', [0.07 .65 .12 .11]);
plotIntensityForceSingle(curTrack,idGroup2(curID2),1,0.01,0.3,MDvin.timeInterval_);
ylim([0 1000])
% force
axes('Position', [0.27 .65 .12 .11]);
plotIntensityForceSingle(curTrack,idGroup2(curID2),2,0.01,0.3,MDvin.timeInterval_);
ylim([0 130])
% edge position
axes('Position', [0.47 .65 .12 .11]);
plotIntensityForceSingle(curTrack,idGroup2(curID2),3,0.01,0.3,MDvin.timeInterval_);
% edge position
axes('Position', [0.67 .65 .12 .11]);
plotIntensityForceSingle(curTrack,idGroup2(curID2),4,0.01,0.3,MDvin.timeInterval_);
% distToEdge
axes('Position', [0.87 .65 .12 .11]);
plotIntensityForceSingle(curTrack,idGroup1(curID2),5,0.01,0.3,MDvin.timeInterval_);
ylim([0 50])

% G3 example
idGroup3 = find(idGroupsVin5{3});
curID3=145;
% curID3=curID3+1;
% curID3=22,67,106,117,131, 145(best); %good candidate
curTrack = tracksNAvin(idGroup3(curID3));
axes('Position', [0.07 .45 .12 .11]);
plotIntensityForceSingle(curTrack,idGroup3(curID3),1,0.01,0.3,MDvin.timeInterval_);
ylim([0 1000])
% force
axes('Position', [0.27 .45 .12 .11]);
plotIntensityForceSingle(curTrack,idGroup3(curID3),2,0.01,0.3,MDvin.timeInterval_);
ylim([0 130])
% edge position
axes('Position', [0.47 .45 .12 .11]);
plotIntensityForceSingle(curTrack,idGroup3(curID3),3,0.01,0.3,MDvin.timeInterval_);
% edge position
axes('Position', [0.67 .45 .12 .11]);
plotIntensityForceSingle(curTrack,idGroup3(curID3),4,0.01,0.3,MDvin.timeInterval_);
% distToEdge
axes('Position', [0.87 .45 .12 .11]);
plotIntensityForceSingle(curTrack,idGroup3(curID3),5,0.01,0.3,MDvin.timeInterval_);
ylim([0 50])

% G7 example
idGroup7 = find(idGroupsVin5{7});
curID7=231;
% curID7=curID7+1;
% curID7=5,15,53, 63, 102, 114, 173, 178, 220 231; %good candidate
curTrack = tracksNAvin(idGroup7(curID7));
axes('Position', [0.07 .25 .12 .11]);
plotIntensityForceSingle(curTrack,idGroup7(curID7),1,0.01,0.3,MDvin.timeInterval_);
ylim([0 1000])
% force
axes('Position', [0.27 .25 .12 .11]);
plotIntensityForceSingle(curTrack,idGroup7(curID7),2,0.01,0.3,MDvin.timeInterval_);
ylim([0 130])
% ylim([0 50])
% edge position
axes('Position', [0.47 .25 .12 .11]);
plotIntensityForceSingle(curTrack,idGroup7(curID7),3,0.01,0.3,MDvin.timeInterval_);
% edge position
axes('Position', [0.67 .25 .12 .11]);
plotIntensityForceSingle(curTrack,idGroup7(curID7),4,0.01,0.3,MDvin.timeInterval_);
% distToEdge
axes('Position', [0.87 .25 .12 .11]);
plotIntensityForceSingle(curTrack,idGroup7(curID7),5,0.01,0.3,MDvin.timeInterval_);
ylim([0 50])

% G9 example
idGroup9 = find(idGroupsVin5{9});
curID9=558;
% curID9=curID9+1;
% curID9=113, 138,513,541; %good candidate
curTrack = tracksNAvin(idGroup9(curID9));
axes('Position', [0.07 .05 .12 .11]);
plotIntensityForceSingle(curTrack,idGroup9(curID9),1,0.01,0.3,MDvin.timeInterval_);
ylim([0 1000])
% force
axes('Position', [0.27 .05 .12 .11]);
plotIntensityForceSingle(curTrack,idGroup9(curID9),2,0.01,0.3,MDvin.timeInterval_);
ylim([0 350])
% ylim([0 50])
% edge position
axes('Position', [0.47 .05 .12 .11]);
plotIntensityForceSingle(curTrack,idGroup9(curID9),3,0.01,0.3,MDvin.timeInterval_);
% edge position
axes('Position', [0.67 .05 .12 .11]);
plotIntensityForceSingle(curTrack,idGroup9(curID9),4,0.01,0.3,MDvin.timeInterval_);
% distToEdge
axes('Position', [0.87 .05 .12 .11]);
plotIntensityForceSingle(curTrack,idGroup9(curID9),5,0.01,0.3,MDvin.timeInterval_);
ylim([0 50])

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
%% Edge advance for five NAs - loading
for ii=1:4
    edgeAdvancePaxG1{ii}=arrayfun(@(x) mean(x.edgeAdvanceDist),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{1}));
    edgeAdvancePaxG2{ii}=arrayfun(@(x) mean(x.edgeAdvanceDist),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{2}));
    edgeAdvancePaxG3{ii}=arrayfun(@(x) mean(x.edgeAdvanceDist),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{3}));
    edgeAdvancePaxG7{ii}=arrayfun(@(x) mean(x.edgeAdvanceDist),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{7}));
    edgeAdvancePaxG9{ii}=arrayfun(@(x) mean(x.edgeAdvanceDist),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{9}));
end

for ii=1:6
    edgeAdvanceVinG1{ii}=arrayfun(@(x) mean(x.edgeAdvanceDist),tracksNAvinAll{ii}(idGroupsVinAll{ii}{1}));
    edgeAdvanceVinG2{ii}=arrayfun(@(x) mean(x.edgeAdvanceDist),tracksNAvinAll{ii}(idGroupsVinAll{ii}{2}));
    edgeAdvanceVinG3{ii}=arrayfun(@(x) mean(x.edgeAdvanceDist),tracksNAvinAll{ii}(idGroupsVinAll{ii}{3}));
    edgeAdvanceVinG7{ii}=arrayfun(@(x) mean(x.edgeAdvanceDist),tracksNAvinAll{ii}(idGroupsVinAll{ii}{7}));
    edgeAdvanceVinG9{ii}=arrayfun(@(x) mean(x.edgeAdvanceDist),tracksNAvinAll{ii}(idGroupsVinAll{ii}{9}));
end

for ii=1:6
    edgeAdvanceTalG1{ii}=arrayfun(@(x) mean(x.edgeAdvanceDist),tracksNAtal{ii}(idGroupsTalAll{ii}{1}));
    edgeAdvanceTalG2{ii}=arrayfun(@(x) mean(x.edgeAdvanceDist),tracksNAtal{ii}(idGroupsTalAll{ii}{2}));
    edgeAdvanceTalG3{ii}=arrayfun(@(x) mean(x.edgeAdvanceDist),tracksNAtal{ii}(idGroupsTalAll{ii}{3}));
    edgeAdvanceTalG7{ii}=arrayfun(@(x) mean(x.edgeAdvanceDist),tracksNAtal{ii}(idGroupsTalAll{ii}{7}));
    edgeAdvanceTalG9{ii}=arrayfun(@(x) mean(x.edgeAdvanceDist),tracksNAtal{ii}(idGroupsTalAll{ii}{9}));
end
edgeAdvancePaxG1t=cell2mat(edgeAdvancePaxG1');
edgeAdvancePaxG2t=cell2mat(edgeAdvancePaxG2');
edgeAdvancePaxG3t=cell2mat(edgeAdvancePaxG3');
edgeAdvancePaxG7t=cell2mat(edgeAdvancePaxG7');
edgeAdvancePaxG9t=cell2mat(edgeAdvancePaxG9');

edgeAdvanceVinG1t=cell2mat(edgeAdvanceVinG1');
edgeAdvanceVinG2t=cell2mat(edgeAdvanceVinG2');
edgeAdvanceVinG3t=cell2mat(edgeAdvanceVinG3');
edgeAdvanceVinG7t=cell2mat(edgeAdvanceVinG7');
edgeAdvanceVinG9t=cell2mat(edgeAdvanceVinG9');

edgeAdvanceTalG1t=cell2mat(edgeAdvanceTalG1');
edgeAdvanceTalG2t=cell2mat(edgeAdvanceTalG2');
edgeAdvanceTalG3t=cell2mat(edgeAdvanceTalG3');
edgeAdvanceTalG7t=cell2mat(edgeAdvanceTalG7');
edgeAdvanceTalG9t=cell2mat(edgeAdvanceTalG9');

%% Potential figure 2 or S1 - Characteristics of five different NAs
% Edge advance
fS1path = '/home2/shan/ownCloud/documents/NascentAdhesionForce/Manuscript/Figures/FigS1fiveNAs';
if ~exist(fS1path,'dir')
    mkdir(fS1path)
end
%% Edge advance for five NAs - plotting for vinculin
pos(3:4) = [200 280];
figure('Position', pos, 'PaperPositionMode', 'auto','Color','w');
% axes('Position', [30/400 520/600 100/400 60/600]); % 
% plot(1:5, [mean(edgeAdvanceTalG1), mean(edgeAdvanceTalG2), mean(edgeAdvanceTalG3),...
%     mean(edgeAdvanceTalG7), mean(edgeAdvanceTalG9)],'ko')
hold on
nameList={['G1,N=' num2str(length(edgeAdvanceVinG1t))] ['G2,N=' num2str(length(edgeAdvanceVinG2t))] ...
    ['G3,N=' num2str(length(edgeAdvanceVinG3t))] ['G7,N=' num2str(length(edgeAdvanceVinG7t))] ['G9,N=' num2str(length(edgeAdvanceVinG9t))]};
edgeAdvanceVin={edgeAdvanceVinG1t,edgeAdvanceVinG2t,edgeAdvanceVinG3t,edgeAdvanceVinG7t,edgeAdvanceVinG9t};
errorBarPlotCellArray(edgeAdvanceVin,nameList)
% boxPlotCellArray(edgeAdvanceVin,nameList)
xlim([0.4 5.6])
% ylim([-0.5 2])
title({'Edge protrusion'; ['in vinculin, M=' num2str(numel(vinFolder))]})
ylabel('Edge protrusion (um)')
% set(gca,'Unit','pixel'); set(gca,'FontSize',7)
print('-depsc', '-loose', [fS1path filesep 'edgeMovementGroupsVin.eps'])
print('-dtiff', '-loose', '-r300', [fS1path filesep 'edgeMovementGroupsVin.tif'])
savefig([fS1path filesep 'edgeMovementGroupsVin.fig'])
print('-depsc', '-loose', [f1path filesep 'edgeMovementGroupsVin.eps'])
print('-dtiff', '-loose', '-r300', [f1path filesep 'edgeMovementGroupsVin.tif'])
savefig([f1path filesep 'edgeMovementGroupsVin.fig'])
%% Edge advance for five NAs - plotting for pax
pos(3:4) = [200 280];
figure('Position', pos, 'PaperPositionMode', 'auto','Color','w');
hold on
nameList={['G1,N=' num2str(length(edgeAdvancePaxG1t))] ['G2,N=' num2str(length(edgeAdvancePaxG2t))] ...
    ['G3,N=' num2str(length(edgeAdvancePaxG3t))] ['G7,N=' num2str(length(edgeAdvancePaxG7t))] ['G9,N=' num2str(length(edgeAdvancePaxG9t))]};
edgeAdvancePax={edgeAdvancePaxG1t,edgeAdvancePaxG2t,edgeAdvancePaxG3t,edgeAdvancePaxG7t,edgeAdvancePaxG9t};
errorBarPlotCellArray(edgeAdvancePax,nameList)
% boxPlotCellArray(edgeAdvancePax,nameList)
xlim([0.4 5.6])
% ylim([-0.5 2])
title({'Edge protrusion'; ['in paxillin, M=' num2str(numel(paxFolder))]})
ylabel('Edge protrusion (um)')
% set(gca,'Unit','pixel'); set(gca,'FontSize',7)
print('-depsc', '-loose', [fS1path filesep 'edgeMovementGroupsPax.eps'])
print('-dtiff', '-loose', '-r300', [fS1path filesep 'edgeMovementGroupsPax.tif'])
savefig([fS1path filesep 'edgeMovementGroupsPax.fig'])
%% Edge advance for five NAs - plotting for talin
pos(3:4) = [200 280];
figure('Position', pos, 'PaperPositionMode', 'auto','Color','w');
hold on
nameList={['G1,N=' num2str(length(edgeAdvanceTalG1t))] ['G2,N=' num2str(length(edgeAdvanceTalG2t))] ...
    ['G3,N=' num2str(length(edgeAdvanceTalG3t))] ['G7,N=' num2str(length(edgeAdvanceTalG7t))] ['G9,N=' num2str(length(edgeAdvanceTalG9t))]};
edgeAdvanceTal={edgeAdvanceTalG1t,edgeAdvanceTalG2t,edgeAdvanceTalG3t,edgeAdvanceTalG7t,edgeAdvanceTalG9t};
errorBarPlotCellArray(edgeAdvanceTal,nameList)
% boxPlotCellArray(edgeAdvanceTal,nameList)
% xlim([0.4 5.6])
% ylim([-0.5 2])
title({'Edge protrusion'; ['in talin, M=' num2str(numel(talFolder))]})
ylabel('Edge protrusion (um)')
% set(gca,'Unit','pixel'); set(gca,'FontSize',7)
print('-depsc', '-loose', [fS1path filesep 'edgeMovementGroupsTal.eps'])
print('-dtiff', '-loose', '-r300', [fS1path filesep 'edgeMovementGroupsTal.tif'])
savefig([fS1path filesep 'edgeMovementGroupsTal.fig'])
%% adhesion advance for five NAs - loading
for ii=1:4
    advancePaxG1{ii}=arrayfun(@(x) mean(x.advanceDist),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{1}));
    advancePaxG2{ii}=arrayfun(@(x) mean(x.advanceDist),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{2}));
    advancePaxG3{ii}=arrayfun(@(x) mean(x.advanceDist),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{3}));
    advancePaxG7{ii}=arrayfun(@(x) mean(x.advanceDist),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{7}));
    advancePaxG9{ii}=arrayfun(@(x) mean(x.advanceDist),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{9}));
end

for ii=1:6
    advanceVinG1{ii}=arrayfun(@(x) mean(x.advanceDist),tracksNAvinAll{ii}(idGroupsVinAll{ii}{1}));
    advanceVinG2{ii}=arrayfun(@(x) mean(x.advanceDist),tracksNAvinAll{ii}(idGroupsVinAll{ii}{2}));
    advanceVinG3{ii}=arrayfun(@(x) mean(x.advanceDist),tracksNAvinAll{ii}(idGroupsVinAll{ii}{3}));
    advanceVinG7{ii}=arrayfun(@(x) mean(x.advanceDist),tracksNAvinAll{ii}(idGroupsVinAll{ii}{7}));
    advanceVinG9{ii}=arrayfun(@(x) mean(x.advanceDist),tracksNAvinAll{ii}(idGroupsVinAll{ii}{9}));
end

for ii=1
    advanceTalG1{ii}=arrayfun(@(x) mean(x.advanceDist),tracksNAtal{ii}(idGroupsTalAll{ii}{1}));
    advanceTalG2{ii}=arrayfun(@(x) mean(x.advanceDist),tracksNAtal{ii}(idGroupsTalAll{ii}{2}));
    advanceTalG3{ii}=arrayfun(@(x) mean(x.advanceDist),tracksNAtal{ii}(idGroupsTalAll{ii}{3}));
    advanceTalG7{ii}=arrayfun(@(x) mean(x.advanceDist),tracksNAtal{ii}(idGroupsTalAll{ii}{7}));
    advanceTalG9{ii}=arrayfun(@(x) mean(x.advanceDist),tracksNAtal{ii}(idGroupsTalAll{ii}{9}));
end
advancePaxG1t=cell2mat(advancePaxG1');
advancePaxG2t=cell2mat(advancePaxG2');
advancePaxG3t=cell2mat(advancePaxG3');
advancePaxG7t=cell2mat(advancePaxG7');
advancePaxG9t=cell2mat(advancePaxG9');

advanceVinG1t=cell2mat(advanceVinG1');
advanceVinG2t=cell2mat(advanceVinG2');
advanceVinG3t=cell2mat(advanceVinG3');
advanceVinG7t=cell2mat(advanceVinG7');
advanceVinG9t=cell2mat(advanceVinG9');

advanceTalG1t=cell2mat(advanceTalG1');
advanceTalG2t=cell2mat(advanceTalG2');
advanceTalG3t=cell2mat(advanceTalG3');
advanceTalG7t=cell2mat(advanceTalG7');
advanceTalG9t=cell2mat(advanceTalG9');

%% adhesion movement for five NAs - plotting for vin
pos(3:4) = [200 280];
figure('Position', pos, 'PaperPositionMode', 'auto','Color','w');
nameList={['G1,N=' num2str(length(advanceVinG1t))] ['G2,N=' num2str(length(advanceVinG2t))] ...
    ['G3,N=' num2str(length(advanceVinG3t))] ['G7,N=' num2str(length(advanceVinG7t))] ['G9,N=' num2str(length(advanceVinG9t))]};
advanceVin={advanceVinG1t,advanceVinG2t,advanceVinG3t,advanceVinG7t,advanceVinG9t};
errorBarPlotCellArray(advanceVin,nameList)
% boxPlotCellArray(advanceVin,nameList)
xlim([0.4 5.6])
% ylim([-0.5 2])
title({'Adhesion movement'; ['in vinculin, M=' num2str(numel(vinFolder))]})
ylabel({'Adhesion movement'; 'toward edge (um)'})
% set(gca,'Unit','pixel'); set(gca,'FontSize',7)
print('-depsc', '-loose', [f1path filesep 'adhesionMovementGroupsVin.eps'])
print('-dtiff', '-loose', '-r300', [f1path filesep 'adhesionMovementGroupsVin.tif'])
savefig([f1path filesep 'adhesionMovementGroupsVin.fig'])
print('-depsc', '-loose', [fS1path filesep 'adhesionMovementGroupsVin.eps'])
print('-dtiff', '-loose', '-r300', [fS1path filesep 'adhesionMovementGroupsVin.tif'])
savefig([fS1path filesep 'adhesionMovementGroupsVin.fig'])
%% adhesion movement for five NAs - plotting for pax
pos(3:4) = [200 280];
figure('Position', pos, 'PaperPositionMode', 'auto','Color','w');
nameList={['G1,N=' num2str(length(advancePaxG1t))] ['G2,N=' num2str(length(advancePaxG2t))] ...
    ['G3,N=' num2str(length(advancePaxG3t))] ['G7,N=' num2str(length(advancePaxG7t))] ['G9,N=' num2str(length(advancePaxG9t))]};
advancePax={advancePaxG1t,advancePaxG2t,advancePaxG3t,advancePaxG7t,advancePaxG9t};
errorBarPlotCellArray(advancePax,nameList)
% boxPlotCellArray(advancePax,nameList)
xlim([0.4 5.6])
% ylim([-0.5 2])
title({'Adhesion movement'; ['in paxillin, M=' num2str(numel(paxFolder))]})
ylabel({'Adhesion movement'; 'toward edge (um)'})
% set(gca,'Unit','pixel'); set(gca,'FontSize',7)
print('-depsc', '-loose', [fS1path filesep 'adhesionMovementGroupsPax.eps'])
print('-dtiff', '-loose', '-r300', [fS1path filesep 'adhesionMovementGroupsPax.tif'])
savefig([fS1path filesep 'adhesionMovementGroupsPax.fig'])
%% adhesion movement for five NAs - plotting for tal
figure('Position', pos, 'PaperPositionMode', 'auto','Color','w');
nameList={['G1,N=' num2str(length(advanceTalG1t))] ['G2,N=' num2str(length(advanceTalG2t))] ...
    ['G3,N=' num2str(length(advanceTalG3t))] ['G7,N=' num2str(length(advanceTalG7t))] ['G9,N=' num2str(length(advanceTalG9t))]};
advanceTal={advanceTalG1t,advanceTalG2t,advanceTalG3t,advanceTalG7t,advanceTalG9t};
errorBarPlotCellArray(advanceTal,nameList)
% boxPlotCellArray(advanceTal,nameList)
xlim([0.4 5.6])
title({'Adhesion movement'; ['in talin, M=' num2str(numel(talFolder))]})
ylabel({'Adhesion movement'; 'toward edge (um)'})
print('-depsc', '-loose', [fS1path filesep 'adhesionMovementGroupsTal.eps'])
print('-dtiff', '-loose', '-r300', [fS1path filesep 'adhesionMovementGroupsTal.tif'])
savefig([fS1path filesep 'adhesionMovementGroupsTal.fig'])
%% distanceToEdge for five NAs - loading
for ii=1:4
    distToEdgePaxG1{ii}=arrayfun(@(x) mean(x.distToEdge),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{1}));
    distToEdgePaxG2{ii}=arrayfun(@(x) mean(x.distToEdge),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{2}));
    distToEdgePaxG3{ii}=arrayfun(@(x) mean(x.distToEdge),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{3}));
    distToEdgePaxG7{ii}=arrayfun(@(x) mean(x.distToEdge),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{7}));
    distToEdgePaxG9{ii}=arrayfun(@(x) mean(x.distToEdge),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{9}));
end

for ii=1:6
    distToEdgeVinG1{ii}=arrayfun(@(x) mean(x.distToEdge),tracksNAvinAll{ii}(idGroupsVinAll{ii}{1}));
    distToEdgeVinG2{ii}=arrayfun(@(x) mean(x.distToEdge),tracksNAvinAll{ii}(idGroupsVinAll{ii}{2}));
    distToEdgeVinG3{ii}=arrayfun(@(x) mean(x.distToEdge),tracksNAvinAll{ii}(idGroupsVinAll{ii}{3}));
    distToEdgeVinG7{ii}=arrayfun(@(x) mean(x.distToEdge),tracksNAvinAll{ii}(idGroupsVinAll{ii}{7}));
    distToEdgeVinG9{ii}=arrayfun(@(x) mean(x.distToEdge),tracksNAvinAll{ii}(idGroupsVinAll{ii}{9}));
end

for ii=1
    distToEdgeTalG1{ii}=arrayfun(@(x) mean(x.distToEdge),tracksNAtal{ii}(idGroupsTalAll{ii}{1}));
    distToEdgeTalG2{ii}=arrayfun(@(x) mean(x.distToEdge),tracksNAtal{ii}(idGroupsTalAll{ii}{2}));
    distToEdgeTalG3{ii}=arrayfun(@(x) mean(x.distToEdge),tracksNAtal{ii}(idGroupsTalAll{ii}{3}));
    distToEdgeTalG7{ii}=arrayfun(@(x) mean(x.distToEdge),tracksNAtal{ii}(idGroupsTalAll{ii}{7}));
    distToEdgeTalG9{ii}=arrayfun(@(x) mean(x.distToEdge),tracksNAtal{ii}(idGroupsTalAll{ii}{9}));
end
distToEdgePaxG1t=cell2mat(distToEdgePaxG1');
distToEdgePaxG2t=cell2mat(distToEdgePaxG2');
distToEdgePaxG3t=cell2mat(distToEdgePaxG3');
distToEdgePaxG7t=cell2mat(distToEdgePaxG7');
distToEdgePaxG9t=cell2mat(distToEdgePaxG9');

distToEdgeVinG1t=cell2mat(distToEdgeVinG1');
distToEdgeVinG2t=cell2mat(distToEdgeVinG2');
distToEdgeVinG3t=cell2mat(distToEdgeVinG3');
distToEdgeVinG7t=cell2mat(distToEdgeVinG7');
distToEdgeVinG9t=cell2mat(distToEdgeVinG9');

distToEdgeTalG1t=cell2mat(distToEdgeTalG1');
distToEdgeTalG2t=cell2mat(distToEdgeTalG2');
distToEdgeTalG3t=cell2mat(distToEdgeTalG3');
distToEdgeTalG7t=cell2mat(distToEdgeTalG7');
distToEdgeTalG9t=cell2mat(distToEdgeTalG9');

%% distToEdgeGroups for five NAs - plotting
figure('Position', pos, 'PaperPositionMode', 'auto','Color','w');
nameList={['G1,N=' num2str(length(distToEdgeVinG1t))] ['G2,N=' num2str(length(distToEdgeVinG2t))] ...
    ['G3,N=' num2str(length(distToEdgeVinG3t))] ['G7,N=' num2str(length(distToEdgeVinG7t))] ['G9,N=' num2str(length(distToEdgeVinG9t))]};
distToEdgeVin={distToEdgeVinG1t,distToEdgeVinG2t,distToEdgeVinG3t,distToEdgeVinG7t,distToEdgeVinG9t};
errorBarPlotCellArray(distToEdgeVin,nameList)
% boxPlotCellArray(distToEdgeVin,nameList)
xlim([0.4 5.6])
% ylim([-0.5 2])
title({'Distance to edge'; ['in vinculin, M=' num2str(numel(vinFolder))]})
ylabel({'Distance to nearest'; ' edge (um)'})
% set(gca,'Unit','pixel'); set(gca,'FontSize',7)
print('-depsc', '-loose', [f1path filesep 'distToEdgeGroupsVin.eps'])
print('-dtiff', '-loose', '-r300', [f1path filesep 'distToEdgeGroupsVin.tif'])
savefig([f1path filesep 'distToEdgeGroupsVin.fig'])
print('-depsc', '-loose', [fS1path filesep 'distToEdgeGroupsVin.eps'])
print('-dtiff', '-loose', '-r300', [fS1path filesep 'distToEdgeGroupsVin.tif'])
savefig([fS1path filesep 'distToEdgeGroupsVin.fig'])
%% distToEdgeGroups for five NAs - plotting for pax
figure('Position', pos, 'PaperPositionMode', 'auto','Color','w');
nameList={['G1,N=' num2str(length(distToEdgePaxG1t))] ['G2,N=' num2str(length(distToEdgePaxG2t))] ...
    ['G3,N=' num2str(length(distToEdgePaxG3t))] ['G7,N=' num2str(length(distToEdgePaxG7t))] ['G9,N=' num2str(length(distToEdgePaxG9t))]};
distToEdgePax={distToEdgePaxG1t,distToEdgePaxG2t,distToEdgePaxG3t,distToEdgePaxG7t,distToEdgePaxG9t};
errorBarPlotCellArray(distToEdgePax,nameList)
% boxPlotCellArray(distToEdgePax,nameList)
xlim([0.4 5.6])
% ylim([-0.5 2])
title({'Distance to edge'; ['in paxillin, M=' num2str(numel(paxFolder))]})
ylabel({'Distance to nearest'; ' edge (um)'})
% set(gca,'Unit','pixel'); set(gca,'FontSize',7)
print('-depsc', '-loose', [fS1path filesep 'distToEdgeGroupsPax.eps'])
print('-dtiff', '-loose', '-r300', [fS1path filesep 'distToEdgeGroupsPax.tif'])
savefig([fS1path filesep 'distToEdgeGroupsPax.fig'])
%% distToEdgeGroups for five NAs - plotting for tal
figure('Position', pos, 'PaperPositionMode', 'auto','Color','w');
nameList={['G1,N=' num2str(length(distToEdgeTalG1t))] ['G2,N=' num2str(length(distToEdgeTalG2t))] ...
    ['G3,N=' num2str(length(distToEdgeTalG3t))] ['G7,N=' num2str(length(distToEdgeTalG7t))] ['G9,N=' num2str(length(distToEdgeTalG9t))]};
distToEdgeTal={distToEdgeTalG1t,distToEdgeTalG2t,distToEdgeTalG3t,distToEdgeTalG7t,distToEdgeTalG9t};
errorBarPlotCellArray(distToEdgeTal,nameList)
% boxPlotCellArray(distToEdgeTal,nameList)
xlim([0.4 5.6])
% ylim([-0.5 2])
title({'Distance to edge'; ['in talin, M=' num2str(numel(talFolder))]})
ylabel({'Distance to nearest'; ' edge (um)'})
% set(gca,'Unit','pixel'); set(gca,'FontSize',7)
print('-depsc', '-loose', [fS1path filesep 'distToEdgeGroupsTal.eps'])
print('-dtiff', '-loose', '-r300', [fS1path filesep 'distToEdgeGroupsTal.tif'])
savefig([fS1path filesep 'distToEdgeGroupsTal.fig'])
%% force increase in G3, G7-protruding and G7-stalling
for ii=1:6
    curForceInc = load([vinFolder{ii} filesep 'data' filesep 'forceIncreaseG3G7.mat']);
    if ii==1
        nameG3G7=curForceInc.nameG3G7;
    end
    forceIncVin{1}{ii} = curForceInc.forceIncreaseG3G7Cell{1};
    forceIncVin{2}{ii} = curForceInc.forceIncreaseG3G7Cell{2};
    forceIncVin{3}{ii} = curForceInc.forceIncreaseG3G7Cell{3};
    edgeAdhInc = load([vinFolder{ii} filesep 'data' filesep 'edgeAdhIncreaseG3G7.mat']);
    if ii==1
        nameG3G7edgeAdh=edgeAdhInc.nameG3G7edgeAdh;
    end
    edgeAdhIncVin{1}{ii} = edgeAdhInc.edgeAdhIncreaseG3G7Cell{1};
    edgeAdhIncVin{2}{ii} = edgeAdhInc.edgeAdhIncreaseG3G7Cell{2};
    edgeAdhIncVin{3}{ii} = edgeAdhInc.edgeAdhIncreaseG3G7Cell{3};
    edgeAdhIncVin{4}{ii} = edgeAdhInc.edgeAdhIncreaseG3G7Cell{4};
    edgeAdhIncVin{5}{ii} = edgeAdhInc.edgeAdhIncreaseG3G7Cell{5};
    edgeAdhIncVin{6}{ii} = edgeAdhInc.edgeAdhIncreaseG3G7Cell{6};
end
forceIncVinGroup=cellfun(@(x) cell2mat(x),forceIncVin,'uniformoutput',false);
edgeAdhIncVinGroup=cellfun(@(x) cell2mat(x),edgeAdhIncVin,'uniformoutput',false);
%% plotting force increase and edge increase in G3 and G7
figure, subplot(1,2,1)
boxPlotCellArray(forceIncVinGroup,nameG3G7,1,false,true);
title('Increase in force in G3 and two separate phases in G7')
ylabel('Change in force (Pa/min)')
subplot(1,2,2)
boxPlotCellArray(edgeAdhIncVinGroup,nameG3G7edgeAdh,1,false,true);
title('Advance in edge and adhesion in G3 and two separate phases in G7')
ylabel('Change in edge (um/min)')

hgsave([fS2path,'/forceIncreaseG3G7'],'-v7.3')
print('-depsc','-loose',[fS2path filesep 'forceIncreaseG3G7.eps']);% histogramPeakLagVinVsTal -transparent

%% force for five NAs - loading
for ii=1:4
    forceMagPaxG1{ii}=arrayfun(@(x) nanmean(x.forceMag),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{1}));
    forceMagPaxG2{ii}=arrayfun(@(x) nanmean(x.forceMag),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{2}));
    forceMagPaxG3{ii}=arrayfun(@(x) nanmean(x.forceMag),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{3}));
    forceMagPaxG7{ii}=arrayfun(@(x) nanmean(x.forceMag),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{7}));
    forceMagPaxG9{ii}=arrayfun(@(x) nanmean(x.forceMag),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{9}));
end

for ii=1:6
    forceMagVinG1{ii}=arrayfun(@(x) nanmean(x.forceMag),tracksNAvinAll{ii}(idGroupsVinAll{ii}{1}));
    forceMagVinG2{ii}=arrayfun(@(x) nanmean(x.forceMag),tracksNAvinAll{ii}(idGroupsVinAll{ii}{2}));
    forceMagVinG3{ii}=arrayfun(@(x) nanmean(x.forceMag),tracksNAvinAll{ii}(idGroupsVinAll{ii}{3}));
    forceMagVinG7{ii}=arrayfun(@(x) nanmean(x.forceMag),tracksNAvinAll{ii}(idGroupsVinAll{ii}{7}));
    forceMagVinG9{ii}=arrayfun(@(x) nanmean(x.forceMag),tracksNAvinAll{ii}(idGroupsVinAll{ii}{9}));
end

for ii=1
    forceMagTalG1{ii}=arrayfun(@(x) nanmean(x.forceMag),tracksNAtal{ii}(idGroupsTalAll{ii}{1}));
    forceMagTalG2{ii}=arrayfun(@(x) nanmean(x.forceMag),tracksNAtal{ii}(idGroupsTalAll{ii}{2}));
    forceMagTalG3{ii}=arrayfun(@(x) nanmean(x.forceMag),tracksNAtal{ii}(idGroupsTalAll{ii}{3}));
    forceMagTalG7{ii}=arrayfun(@(x) nanmean(x.forceMag),tracksNAtal{ii}(idGroupsTalAll{ii}{7}));
    forceMagTalG9{ii}=arrayfun(@(x) nanmean(x.forceMag),tracksNAtal{ii}(idGroupsTalAll{ii}{9}));
end
forceMagPaxG1t=cell2mat(forceMagPaxG1');
forceMagPaxG2t=cell2mat(forceMagPaxG2');
forceMagPaxG3t=cell2mat(forceMagPaxG3');
forceMagPaxG7t=cell2mat(forceMagPaxG7');
forceMagPaxG9t=cell2mat(forceMagPaxG9');

vinInterest=[2 3 6];
forceMagVinG1t=cell2mat(forceMagVinG1(vinInterest)');
forceMagVinG2t=cell2mat(forceMagVinG2(vinInterest)');
forceMagVinG3t=cell2mat(forceMagVinG3(vinInterest)');
forceMagVinG7t=cell2mat(forceMagVinG7(vinInterest)');
forceMagVinG9t=cell2mat(forceMagVinG9(vinInterest)');

forceMagTalG1t=cell2mat(forceMagTalG1');
forceMagTalG2t=cell2mat(forceMagTalG2');
forceMagTalG3t=cell2mat(forceMagTalG3');
forceMagTalG7t=cell2mat(forceMagTalG7');
forceMagTalG9t=cell2mat(forceMagTalG9');

%% forceMagGroups for five NAs - plotting
figure('Position', pos, 'PaperPositionMode', 'auto','Color','w');
nameList={['G1,N=' num2str(length(forceMagVinG1t))] ['G2,N=' num2str(length(forceMagVinG2t))] ...
    ['G3,N=' num2str(length(forceMagVinG3t))] ['G7,N=' num2str(length(forceMagVinG7t))] ['G9,N=' num2str(length(forceMagVinG9t))]};
forceMagVin={forceMagVinG1t,forceMagVinG2t,forceMagVinG3t,forceMagVinG7t,forceMagVinG9t};
errorBarPlotCellArray(forceMagVin,nameList)
xlim([0.4 5.6])
% ylim([-0.5 2])
title({'Traction magnitude'; ['in vinculin experiment,M=' num2str(numel(vinFolder))]})
ylabel({'Traction magnitude (Pa)'})
% set(gca,'Unit','pixel'); set(gca,'FontSize',7)
print('-depsc', '-loose', [f1path filesep 'forceMagGroupsVin.eps'])
print('-dtiff', '-loose', '-r300', [f1path filesep 'forceMagGroupsVin.tif'])
savefig([f1path filesep 'forceMagGroupsVin.fig'])
print('-depsc', '-loose', [fS1path filesep 'forceMagGroupsVin.eps'])
print('-dtiff', '-loose', '-r300', [fS1path filesep 'forceMagGroupsVin.tif'])
savefig([fS1path filesep 'forceMagGroupsVin.fig'])
%% forceMagGroups for five NAs - plotting for talin
figure('Position', pos, 'PaperPositionMode', 'auto','Color','w');
nameList={['G1,N=' num2str(length(forceMagTalG1t))] ['G2,N=' num2str(length(forceMagTalG2t))] ...
    ['G3,N=' num2str(length(forceMagTalG3t))] ['G7,N=' num2str(length(forceMagTalG7t))] ['G9,N=' num2str(length(forceMagTalG9t))]};
forceMagTal={forceMagTalG1t,forceMagTalG2t,forceMagTalG3t,forceMagTalG7t,forceMagTalG9t};
errorBarPlotCellArray(forceMagTal,nameList)
xlim([0.4 5.6])
% ylim([-0.5 2])
title({'Traction magnitude'; ['in talin experiment,M=' num2str(numel(talFolder))]})
ylabel({'Traction magnitude (Pa)'})
print('-depsc', '-loose', [fS1path filesep 'forceMagGroupsTal.eps'])
print('-dtiff', '-loose', '-r300', [fS1path filesep 'forceMagGroupsTal.tif'])
savefig([fS1path filesep 'forceMagGroupsTal.fig'])
%% forceMagGroups for five NAs - plotting for paxillin
figure('Position', pos, 'PaperPositionMode', 'auto','Color','w');
nameList={['G1,N=' num2str(length(forceMagPaxG1t))] ['G2,N=' num2str(length(forceMagPaxG2t))] ...
    ['G3,N=' num2str(length(forceMagPaxG3t))] ['G7,N=' num2str(length(forceMagPaxG7t))] ['G9,N=' num2str(length(forceMagPaxG9t))]};
forceMagPax={forceMagPaxG1t,forceMagPaxG2t,forceMagPaxG3t,forceMagPaxG7t,forceMagPaxG9t};
errorBarPlotCellArray(forceMagPax,nameList)
xlim([0.4 5.6])
% ylim([-0.5 2])
title({'Traction magnitude'; ['in paxillin experiment,M=' num2str(numel(paxFolder))]})
ylabel({'Traction magnitude (Pa)'})
print('-depsc', '-loose', [fS1path filesep 'forceMagGroupsPax.eps'])
print('-dtiff', '-loose', '-r300', [fS1path filesep 'forceMagGroupsPax.tif'])
savefig([fS1path filesep 'forceMagGroupsPax.fig'])
%% Force slope for five NAs
for ii=1:6
    curForceSlope = load([vinFolder{ii} filesep 'data' filesep 'forceSlopeAllGroups.mat']);
    forceSlopeVin{1}{ii} = curForceSlope.forceSlope{1};
    forceSlopeVin{2}{ii} = curForceSlope.forceSlope{2};
    forceSlopeVin{3}{ii} = curForceSlope.forceSlope{3};
    forceSlopeVin{4}{ii} = curForceSlope.forceSlope{7};
    forceSlopeVin{5}{ii} = curForceSlope.forceSlope{9};
end
forceSlopeVinGroup=cellfun(@(x) cell2mat(x'),forceSlopeVin,'uniformoutput',false);
%% plotting force increase and edge increase in G3 and G7
figure, 
% barPlotCellArray(forceSlopeVinGroup,nameList,1,false,true);
errorBarPlotCellArray(forceSlopeVinGroup,nameList);
title({'Force growth rate'; ['M=' num2str(numel(vinFolder))]})
ylabel({'Force growth rate (Pa/min)'})

hgsave([fS2path,'/forceGrowthRate5NAs'],'-v7.3')
print('-depsc','-loose',[fS2path filesep 'forceGrowthRate5NAs.eps']);% histogramPeakLagVinVsTal -transparent


%% ampTotal for five NAs - loading
for ii=1:4
    ampTotalPaxG1{ii}=arrayfun(@(x) nanmean(x.ampTotal),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{1}));
    ampTotalPaxG2{ii}=arrayfun(@(x) nanmean(x.ampTotal),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{2}));
    ampTotalPaxG3{ii}=arrayfun(@(x) nanmean(x.ampTotal),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{3}));
    ampTotalPaxG7{ii}=arrayfun(@(x) nanmean(x.ampTotal),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{7}));
    ampTotalPaxG9{ii}=arrayfun(@(x) nanmean(x.ampTotal),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{9}));
end

for ii=1:6
    ampTotalVinG1{ii}=arrayfun(@(x) nanmean(x.ampTotal),tracksNAvinAll{ii}(idGroupsVinAll{ii}{1}));
    ampTotalVinG2{ii}=arrayfun(@(x) nanmean(x.ampTotal),tracksNAvinAll{ii}(idGroupsVinAll{ii}{2}));
    ampTotalVinG3{ii}=arrayfun(@(x) nanmean(x.ampTotal),tracksNAvinAll{ii}(idGroupsVinAll{ii}{3}));
    ampTotalVinG7{ii}=arrayfun(@(x) nanmean(x.ampTotal),tracksNAvinAll{ii}(idGroupsVinAll{ii}{7}));
    ampTotalVinG9{ii}=arrayfun(@(x) nanmean(x.ampTotal),tracksNAvinAll{ii}(idGroupsVinAll{ii}{9}));
end

for ii=1
    ampTotalTalG1{ii}=arrayfun(@(x) nanmean(x.ampTotal),tracksNAtal{ii}(idGroupsTalAll{ii}{1}));
    ampTotalTalG2{ii}=arrayfun(@(x) nanmean(x.ampTotal),tracksNAtal{ii}(idGroupsTalAll{ii}{2}));
    ampTotalTalG3{ii}=arrayfun(@(x) nanmean(x.ampTotal),tracksNAtal{ii}(idGroupsTalAll{ii}{3}));
    ampTotalTalG7{ii}=arrayfun(@(x) nanmean(x.ampTotal),tracksNAtal{ii}(idGroupsTalAll{ii}{7}));
    ampTotalTalG9{ii}=arrayfun(@(x) nanmean(x.ampTotal),tracksNAtal{ii}(idGroupsTalAll{ii}{9}));
end
ampTotalPaxG1t=cell2mat(ampTotalPaxG1');
ampTotalPaxG2t=cell2mat(ampTotalPaxG2');
ampTotalPaxG3t=cell2mat(ampTotalPaxG3');
ampTotalPaxG7t=cell2mat(ampTotalPaxG7');
ampTotalPaxG9t=cell2mat(ampTotalPaxG9');

ampTotalVinG1t=cell2mat(ampTotalVinG1');
ampTotalVinG2t=cell2mat(ampTotalVinG2');
ampTotalVinG3t=cell2mat(ampTotalVinG3');
ampTotalVinG7t=cell2mat(ampTotalVinG7');
ampTotalVinG9t=cell2mat(ampTotalVinG9');

ampTotalTalG1t=cell2mat(ampTotalTalG1');
ampTotalTalG2t=cell2mat(ampTotalTalG2');
ampTotalTalG3t=cell2mat(ampTotalTalG3');
ampTotalTalG7t=cell2mat(ampTotalTalG7');
ampTotalTalG9t=cell2mat(ampTotalTalG9');

%% ampTotalGroups for five NAs - plotting - vinculin
pos(3:4) = [200 250];
figure('Position', pos, 'PaperPositionMode', 'auto','Color','w');
nameList={['G1,N=' num2str(length(ampTotalVinG1t))] ['G2,N=' num2str(length(ampTotalVinG2t))] ...
    ['G3,N=' num2str(length(ampTotalVinG3t))] ['G7,N=' num2str(length(ampTotalVinG7t))] ['G9,N=' num2str(length(ampTotalVinG9t))]};
ampTotalVin={ampTotalVinG1t,ampTotalVinG2t,ampTotalVinG3t,ampTotalVinG7t,ampTotalVinG9t};
errorBarPlotCellArray(ampTotalVin,nameList)
xlim([0.4 5.6])
% ylim([-0.5 2])
title({'Vinculin amplitude'; ['M=' num2str(numel(vinFolder))]})
ylabel({'fluorescence intensity (a.u.)'})
% set(gca,'Unit','pixel'); set(gca,'FontSize',7)
print('-depsc', '-loose', [f1path filesep 'ampTotalGroupsVin.eps'])
print('-dtiff', '-loose', '-r300', [f1path filesep 'ampTotalGroupsVin.tif'])
savefig([f1path filesep 'ampTotalGroupsVin.fig'])
print('-depsc', '-loose', [fS1path filesep 'ampTotalGroupsVin.eps'])
print('-dtiff', '-loose', '-r300', [fS1path filesep 'ampTotalGroupsVin.tif'])
savefig([fS1path filesep 'ampTotalGroupsVin.fig'])
%% ampTotalGroups for five NAs - plotting - talin
figure('Position', pos, 'PaperPositionMode', 'auto','Color','w');
nameList={['G1,N=' num2str(length(ampTotalTalG1t))] ['G2,N=' num2str(length(ampTotalTalG2t))] ...
    ['G3,N=' num2str(length(ampTotalTalG3t))] ['G7,N=' num2str(length(ampTotalTalG7t))] ['G9,N=' num2str(length(ampTotalTalG9t))]};
ampTotalTal={ampTotalTalG1t,ampTotalTalG2t,ampTotalTalG3t,ampTotalTalG7t,ampTotalTalG9t};
errorBarPlotCellArray(ampTotalTal,nameList)
xlim([0.4 5.6])
% ylim([-0.5 2])
title({'Talin amplitude '; ['M=' num2str(numel(talFolder))]})
ylabel({'fluorescence intensity (a.u.)'})
% set(gca,'Unit','pixel'); set(gca,'FontSize',7)
print('-depsc', '-loose', [fS1path filesep 'ampTotalGroupsTal.eps'])
print('-dtiff', '-loose', '-r300', [fS1path filesep 'ampTotalGroupsTal.tif'])
savefig([fS1path filesep 'ampTotalGroupsTal.fig'])
%% ampTotalGroups for five NAs - plotting - paxillin
figure('Position', pos, 'PaperPositionMode', 'auto','Color','w');
nameList={['G1,N=' num2str(length(ampTotalPaxG1t))] ['G2,N=' num2str(length(ampTotalPaxG2t))] ...
    ['G3,N=' num2str(length(ampTotalPaxG3t))] ['G7,N=' num2str(length(ampTotalPaxG7t))] ['G9,N=' num2str(length(ampTotalPaxG9t))]};
ampTotalPax={ampTotalPaxG1t,ampTotalPaxG2t,ampTotalPaxG3t,ampTotalPaxG7t,ampTotalPaxG9t};
errorBarPlotCellArray(ampTotalPax,nameList)
xlim([0.4 5.6])
% ylim([-0.5 2])
title({'Paxillin amplitude '; ['M=' num2str(numel(paxFolder))]})
ylabel({'fluorescence intensity (a.u.)'})
% set(gca,'Unit','pixel'); set(gca,'FontSize',7)
print('-depsc', '-loose', [fS1path filesep 'ampTotalGroupsPax.eps'])
print('-dtiff', '-loose', '-r300', [fS1path filesep 'ampTotalGroupsPax.tif'])
savefig([fS1path filesep 'ampTotalGroupsPax.fig'])
%% LifeTime for five NAs - loading
for ii=1:4
    lifeTimePaxG1{ii}=arrayfun(@(x) nanmean(x.lifeTime),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{1}));
    lifeTimePaxG2{ii}=arrayfun(@(x) nanmean(x.lifeTime),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{2}));
    lifeTimePaxG3{ii}=arrayfun(@(x) nanmean(x.lifeTime),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{3}));
    lifeTimePaxG7{ii}=arrayfun(@(x) nanmean(x.lifeTime),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{7}));
    lifeTimePaxG9{ii}=arrayfun(@(x) nanmean(x.lifeTime),tracksNApaxAll{ii}(idGroupsPaxAll{ii}{9}));
end

for ii=1:6
    lifeTimeVinG1{ii}=arrayfun(@(x) nanmean(x.lifeTime),tracksNAvinAll{ii}(idGroupsVinAll{ii}{1}));
    lifeTimeVinG2{ii}=arrayfun(@(x) nanmean(x.lifeTime),tracksNAvinAll{ii}(idGroupsVinAll{ii}{2}));
    lifeTimeVinG3{ii}=arrayfun(@(x) nanmean(x.lifeTime),tracksNAvinAll{ii}(idGroupsVinAll{ii}{3}));
    lifeTimeVinG7{ii}=arrayfun(@(x) nanmean(x.lifeTime),tracksNAvinAll{ii}(idGroupsVinAll{ii}{7}));
    lifeTimeVinG9{ii}=arrayfun(@(x) nanmean(x.lifeTime),tracksNAvinAll{ii}(idGroupsVinAll{ii}{9}));
end

for ii=1
    lifeTimeTalG1{ii}=arrayfun(@(x) nanmean(x.lifeTime),tracksNAtalAll{ii}(idGroupsTalAll{ii}{1}));
    lifeTimeTalG2{ii}=arrayfun(@(x) nanmean(x.lifeTime),tracksNAtalAll{ii}(idGroupsTalAll{ii}{2}));
    lifeTimeTalG3{ii}=arrayfun(@(x) nanmean(x.lifeTime),tracksNAtalAll{ii}(idGroupsTalAll{ii}{3}));
    lifeTimeTalG7{ii}=arrayfun(@(x) nanmean(x.lifeTime),tracksNAtalAll{ii}(idGroupsTalAll{ii}{7}));
    lifeTimeTalG9{ii}=arrayfun(@(x) nanmean(x.lifeTime),tracksNAtalAll{ii}(idGroupsTalAll{ii}{9}));
end
paxInterest=[1 ];
lifeTimePaxG1t=cell2mat(lifeTimePaxG1(paxInterest)');
lifeTimePaxG2t=cell2mat(lifeTimePaxG2(paxInterest)');
lifeTimePaxG3t=cell2mat(lifeTimePaxG3(paxInterest)');
lifeTimePaxG7t=cell2mat(lifeTimePaxG7(paxInterest)');
lifeTimePaxG9t=cell2mat(lifeTimePaxG9(paxInterest)');

lifeTimeVinG1t=cell2mat(lifeTimeVinG1');
lifeTimeVinG2t=cell2mat(lifeTimeVinG2');
lifeTimeVinG3t=cell2mat(lifeTimeVinG3');
lifeTimeVinG7t=cell2mat(lifeTimeVinG7');
lifeTimeVinG9t=cell2mat(lifeTimeVinG9');

lifeTimeTalG1t=cell2mat(lifeTimeTalG1');
lifeTimeTalG2t=cell2mat(lifeTimeTalG2');
lifeTimeTalG3t=cell2mat(lifeTimeTalG3');
lifeTimeTalG7t=cell2mat(lifeTimeTalG7');
lifeTimeTalG9t=cell2mat(lifeTimeTalG9');
%% lifeTimeGroups for five NAs - plotting - vinculin
figure('Position', pos, 'PaperPositionMode', 'auto','Color','w');
nameList={['G1,N=' num2str(length(lifeTimeVinG1t))] ['G2,N=' num2str(length(lifeTimeVinG2t))] ...
    ['G3,N=' num2str(length(lifeTimeVinG3t))] ['G7,N=' num2str(length(lifeTimeVinG7t))] ['G9,N=' num2str(length(lifeTimeVinG9t))]};
lifeTimeVin={lifeTimeVinG1t,lifeTimeVinG2t,lifeTimeVinG3t,lifeTimeVinG7t,lifeTimeVinG9t};
errorBarPlotCellArray(lifeTimeVin,nameList)
xlim([0.4 5.6])
% ylim([-0.5 2])
title({'Vinculin lifetime'; ['M=' num2str(numel(vinFolder))]})
ylabel({'Lifetime (sec)'})
% set(gca,'Unit','pixel'); set(gca,'FontSize',7)
print('-depsc', '-loose', [f1path filesep 'lifeTimeGroupsVin.eps'])
print('-dtiff', '-loose', '-r300', [f1path filesep 'lifeTimeGroupsVin.tif'])
savefig([f1path filesep 'lifeTimeGroupsVin.fig'])
print('-depsc', '-loose', [fS1path filesep 'lifeTimeGroupsVin.eps'])
print('-dtiff', '-loose', '-r300', [fS1path filesep 'lifeTimeGroupsVin.tif'])
savefig([fS1path filesep 'lifeTimeGroupsVin.fig'])
%% lifeTimeGroups for five NAs - plotting - talin
figure('Position', pos, 'PaperPositionMode', 'auto','Color','w');
nameList={['G1,N=' num2str(length(lifeTimeTalG1t))] ['G2,N=' num2str(length(lifeTimeTalG2t))] ...
    ['G3,N=' num2str(length(lifeTimeTalG3t))] ['G7,N=' num2str(length(lifeTimeTalG7t))] ['G9,N=' num2str(length(lifeTimeTalG9t))]};
lifeTimeTal={lifeTimeTalG1t,lifeTimeTalG2t,lifeTimeTalG3t,lifeTimeTalG7t,lifeTimeTalG9t};
errorBarPlotCellArray(lifeTimeTal,nameList)
xlim([0.4 5.6])
% ylim([-0.5 2])
title({'Talin lifetime'; ['M=' num2str(numel(talFolder))]})
ylabel({'Lifetime (sec)'})
% set(gca,'Unit','pixel'); set(gca,'FontSize',7)
print('-depsc', '-loose', [fS1path filesep 'lifeTimeGroupsTal.eps'])
print('-dtiff', '-loose', '-r300', [fS1path filesep 'lifeTimeGroupsTal.tif'])
savefig([fS1path filesep 'lifeTimeGroupsTal.fig'])
%% lifeTimeGroups for five NAs - plotting - paxillin
figure('Position', pos, 'PaperPositionMode', 'auto','Color','w');
nameList={['G1,N=' num2str(length(lifeTimePaxG1t))] ['G2,N=' num2str(length(lifeTimePaxG2t))] ...
    ['G3,N=' num2str(length(lifeTimePaxG3t))] ['G7,N=' num2str(length(lifeTimePaxG7t))] ['G9,N=' num2str(length(lifeTimePaxG9t))]};
lifeTimePax={lifeTimePaxG1t,lifeTimePaxG2t,lifeTimePaxG3t,lifeTimePaxG7t,lifeTimePaxG9t};
errorBarPlotCellArray(lifeTimePax,nameList)
xlim([0.4 5.6])
% ylim([-0.5 2])
title({'Paxillin lifetime'; ['M=' num2str(numel(paxFolder))]})
ylabel({'Lifetime (sec)'})
% set(gca,'Unit','pixel'); set(gca,'FontSize',7)
print('-depsc', '-loose', [fS1path filesep 'lifeTimeGroupsPax.eps'])
print('-dtiff', '-loose', '-r300', [fS1path filesep 'lifeTimeGroupsPax.tif'])
savefig([fS1path filesep 'lifeTimeGroupsPax.fig'])
%% Fig 3 Heterogeneity in force time-course in G1 adhesions and time-lag quantification method
%% 2a G1 with adhesion channel and force channel
pos(3:4) = [600 600];
unifiedWidthNanometer = 300*MDtal.pixelSize_;
widthPixelTalin = unifiedWidthNanometer/MDtal.pixelSize_;
CurrentFrameTal=100;

figure('Position', pos, 'PaperPositionMode', 'auto','Color','w');
axes('Position', [0 5/6 1/6 1/6]); % 
imshow(imcomplement(mean(talImgStack(:,:,CurrentFrameTal-1:CurrentFrameTal+1),3)),[-900 -100])
% set(gca,'XLim',[30 30+widthPixelTalin],'YLim',[30 30+widthPixelTalin])
set(gca,'XLim',[64 64+widthPixelTalin],'YLim',[200 200+widthPixelTalin])
hold on
line([64+10 64+10+round(5000/MDtal.pixelSize_)],[200+widthPixelTalin-20 200+widthPixelTalin-20],'LineWidth',2,'Color','k')
text(64+10, 200+widthPixelTalin-20-20,'5 um','Color','k','Fontsize',7)
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

% CurrentFrame = 197;
idCurrent = arrayfun(@(x) x.startingFrameExtra<=CurrentFrameTal & x.endingFrameExtra>=CurrentFrameTal,tracksNAtal);
markerSize=2;
% plot(arrayfun(@(x) x.xCoord(CurrentFrame),tracksNAtal(idCurrent & (idGroupsTal{1} | idGroupsTal{7}))),arrayfun(@(x) x.yCoord(CurrentFrame),tracksNAtal(idCurrent & (idGroupsTal{1} | idGroupsTal{7}))),'go','MarkerSize',markerSize)
plot(arrayfun(@(x) x.xCoord(CurrentFrameTal),tracksNAtal(idCurrent & (idGroupsTal{1}))),arrayfun(@(x) x.yCoord(CurrentFrameTal),tracksNAtal(idCurrent & (idGroupsTal{1}))),'go','MarkerSize',markerSize)
xmat = cell2mat(arrayfun(@(x) x.xCoord(1:CurrentFrameTal),tracksNAtal(idCurrent & idGroupsTal{1}),'UniformOutput',false));
ymat = cell2mat(arrayfun(@(x) x.yCoord(1:CurrentFrameTal),tracksNAtal(idCurrent & idGroupsTal{1}),'UniformOutput',false));
plot(xmat',ymat','g','linewidth',0.25)
plot(arrayfun(@(x) x.xCoord(CurrentFrameTal),tracksNAtal(idCurrent & (idGroupsTal{2}))),arrayfun(@(x) x.yCoord(CurrentFrameTal),tracksNAtal(idCurrent & (idGroupsTal{2}))),'bo','MarkerSize',markerSize)

axes('Position', [0 4/6 1/6 1/6]); % 
imshow(mean(talForceStack(:,:,CurrentFrameTal-1:CurrentFrameTal+1),3),[20 400]), colormap(gca,mycmap),hold on
set(gca,'XLim',[64 64+widthPixelTalin],'YLim',[200 200+widthPixelTalin])
% plot(arrayfun(@(x) x.xCoord(CurrentFrameTal),tracksNAtal(idCurrent & (idGroupsTal{1} | idGroupsTal{7}))),arrayfun(@(x) x.yCoord(CurrentFrameTal),tracksNAtal(idCurrent & (idGroupsTal{1} | idGroupsTal{7}))),'go','MarkerSize',markerSize)
plot(arrayfun(@(x) x.xCoord(CurrentFrameTal),tracksNAtal(idCurrent & (idGroupsTal{1}))),arrayfun(@(x) x.yCoord(CurrentFrameTal),tracksNAtal(idCurrent & (idGroupsTal{1}))),'go','MarkerSize',markerSize)
plot(xmat',ymat','g','linewidth',0.25)
plot(arrayfun(@(x) x.xCoord(CurrentFrameTal),tracksNAtal(idCurrent & (idGroupsTal{2}))),arrayfun(@(x) x.yCoord(CurrentFrameTal),tracksNAtal(idCurrent & (idGroupsTal{2}))),'wo','MarkerSize',markerSize)
%% Vinculin
axes('Position', [0 3/6 1/6 1/6]); % 
CurrentFrame = 250;
imshow(imcomplement(mean(vinImgStack(:,:,CurrentFrame-2:CurrentFrame+2),3)),[-4e4 -3800])
widthPixelVinculin = unifiedWidthNanometer/MDvin5.pixelSize_;
set(gca,'XLim',[100 100+widthPixelVinculin],'YLim',[32 32+widthPixelVinculin])
% set(gca,'XLim',[81 81+widthPixelVinculin],'YLim',[72 72+widthPixelVinculin])
hold on
idCurrent = arrayfun(@(x) x.startingFrameExtra<=CurrentFrame & x.endingFrameExtra>=CurrentFrame,tracksNAvin);
% plot(arrayfun(@(x) x.xCoord(CurrentFrame),tracksNAvin(idCurrent & (idGroupsVin5{1} | idGroupsVin5{7}))),arrayfun(@(x) x.yCoord(CurrentFrame),tracksNAvin(idCurrent & (idGroupsVin5{1} | idGroupsVin5{7}))),'go','MarkerSize',markerSize)
plot(arrayfun(@(x) x.xCoord(CurrentFrame),tracksNAvin(idCurrent & (idGroupsVin5{1}))),arrayfun(@(x) x.yCoord(CurrentFrame),tracksNAvin(idCurrent & (idGroupsVin5{1}))),'go','MarkerSize',markerSize)
plot(arrayfun(@(x) x.xCoord(CurrentFrame),tracksNAvin(idCurrent & (idGroupsVin5{2}))),arrayfun(@(x) x.yCoord(CurrentFrame),tracksNAvin(idCurrent & (idGroupsVin5{2}))),'bo','MarkerSize',markerSize)


% axes('Position', [0.33 .6 .15 .17]); % 
axes('Position', [0 2/6  1/6 1/6]); % 
imshow(mean(vinForceStack(:,:,CurrentFrame-2:CurrentFrame+2),3),[20 1500]), colormap(gca,mycmap),hold on
set(gca,'XLim',[100 100+widthPixelVinculin],'YLim',[32 32+widthPixelVinculin])
% set(gca,'XLim',[81 81+widthPixelVinculin],'YLim',[72 72+widthPixelVinculin])
plot(arrayfun(@(x) x.xCoord(CurrentFrame),tracksNAvin(idCurrent & (idGroupsVin5{1}))),arrayfun(@(x) x.yCoord(CurrentFrame),tracksNAvin(idCurrent & (idGroupsVin5{1}))),'go','MarkerSize',markerSize)
plot(arrayfun(@(x) x.xCoord(CurrentFrame),tracksNAvin(idCurrent & (idGroupsVin5{2}))),arrayfun(@(x) x.yCoord(CurrentFrame),tracksNAvin(idCurrent & (idGroupsVin5{2}))),'wo','MarkerSize',markerSize)
% plot(arrayfun(@(x) x.xCoord(CurrentFrame),tracksNAvin(idCurrent & (idGroupsVin5{1} | idGroupsVin5{7}))),arrayfun(@(x) x.yCoord(CurrentFrame),tracksNAvin(idCurrent & (idGroupsVin5{1} | idGroupsVin5{7}))),'go','MarkerSize',markerSize)
%% Life time histogram (testing) - get all cells
% G1 histogram
figure, 
subplot(3,2,1), histogram(lifeTimePaxG1t), title('Pax G1')
subplot(3,2,2), histogram(lifeTimePaxG2t), title('Pax G2')
subplot(3,2,3), histogram(lifeTimeVinG1t), title('Vin G1')
subplot(3,2,4), histogram(lifeTimeVinG2t), title('Vin G2')
subplot(3,2,5), histogram(lifeTimeTalG1t), title('Tal G1')
subplot(3,2,6), histogram(lifeTimeTalG2t), title('Tal G2')
%% Paxillin
% axes('Position', [0.66 .8 .15 .17]); % 
if ~exist('tracksNApax','var')
    tic
    tracksNApax = load([paxFolder{4} filesep 'data' filesep 'tracksNA.mat'],'tracksNA');
    tracksNApax = tracksNApax.tracksNA; toc
end
if ~exist('idGroupsPax','var')
    idGroupsPaxStruct = load([paxFolder{4} filesep 'data' filesep 'idsClassified.mat'],'idGroup1','idGroup2','idGroup3',...
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

CurrentFrame = 275;
axes('Position', [0 1/6 1/6 1/6]); % 
imshow(imcomplement(mean(paxImgStack(:,:,CurrentFrame-2:CurrentFrame+2),3)),[-800 -100])
widthPixelPaxilln = unifiedWidthNanometer/MDpax.pixelSize_;
set(gca,'XLim',[30 30+widthPixelPaxilln],'YLim',[30 30+widthPixelPaxilln])
hold on
idCurrent = arrayfun(@(x) x.startingFrameExtra<=CurrentFrame & x.endingFrameExtra>=CurrentFrame,tracksNApax);
plot(arrayfun(@(x) x.xCoord(CurrentFrame),tracksNApax(idCurrent & (idGroupsPax{1}))),arrayfun(@(x) x.yCoord(CurrentFrame),tracksNApax(idCurrent & (idGroupsPax{1}))),'go','MarkerSize',markerSize)
plot(arrayfun(@(x) x.xCoord(CurrentFrame),tracksNApax(idCurrent & (idGroupsPax{2}))),arrayfun(@(x) x.yCoord(CurrentFrame),tracksNApax(idCurrent & (idGroupsPax{2}))),'bo','MarkerSize',markerSize)
% plot(arrayfun(@(x) x.xCoord(CurrentFrame),tracksNApax(idCurrent & (idGroupsPax{1} |  idGroupsPax{7}))),arrayfun(@(x) x.yCoord(CurrentFrame),tracksNApax(idCurrent & (idGroupsPax{1} |  idGroupsPax{7}))),'go','MarkerSize',markerSize)

% axes('Position', [0.66 .6 .15 .17]); % 
axes('Position', [0 0 1/6 1/6]); % 
imshow(mean(paxForceStack(:,:,295:300),3),[20 1000]), colormap(gca,mycmap),hold on
set(gca,'XLim',[30 30+widthPixelPaxilln],'YLim',[30 30+widthPixelPaxilln])
plot(arrayfun(@(x) x.xCoord(CurrentFrame),tracksNApax(idCurrent & (idGroupsPax{1}))),arrayfun(@(x) x.yCoord(CurrentFrame),tracksNApax(idCurrent & (idGroupsPax{1}))),'go','MarkerSize',markerSize)
plot(arrayfun(@(x) x.xCoord(CurrentFrame),tracksNApax(idCurrent & (idGroupsPax{2}))),arrayfun(@(x) x.yCoord(CurrentFrame),tracksNApax(idCurrent & (idGroupsPax{2}))),'wo','MarkerSize',markerSize)
% plot(arrayfun(@(x) x.xCoord(CurrentFrame),tracksNApax(idCurrent & (idGroupsPax{1} |  idGroupsPax{7}))),arrayfun(@(x) x.yCoord(CurrentFrame),tracksNApax(idCurrent & (idGroupsPax{1} |  idGroupsPax{7}))),'go','MarkerSize',markerSize)
% %% recalculate lifeTime
% tracksNAtal = recalculateLifeTimeTracks(tracksNAtal);
%% intensity and force profiles-talin - now this is in supplementary
% axes('Position', [0.21 .85 .11 .12]); % 
pos(3:4) = [600 600];
figure('Position', pos, 'PaperPositionMode', 'auto','Color','w');
axes('Position', [130/600 530/600 65/600 65/600]); % 
% plotIntensityForce(tracksNAtal(idGroupsTal{1} | idGroupsTal{7}),[],false,false,'UseCurrentAxis',true,'Source',{'lifeTime'},'plotCohorts',true,'tInterval',MDtal.timeInterval_,'prePostFrames',10)
plotIntensityForce(tracksNAtal(idGroupsTal{1}),[],false,false,'UseCurrentAxis',true,'Source',{'ampTotal'},'plotCohorts',true,'tInterval',MDtal.timeInterval_,'prePostFrames',10)
xlim([-5 84])
% set(findobj(gca,'Type','Text'),'FontSize',7)
% axes('Position', [0.21 .65 .11 .12]); % 
axes('Position', [130/600 430/600 65/600 65/600]); % 
% plotIntensityForce(tracksNAtal(idGroupsTal{1} | idGroupsTal{7}),[],false,false,'UseCurrentAxis',true,'Source',{'forceMag'},'plotCohorts',true,'tInterval',MDtal.timeInterval_,'prePostFrames',10)
plotIntensityForce(tracksNAtal(idGroupsTal{1}),[],false,false,'UseCurrentAxis',true,'Source',{'forceMag'},'plotCohorts',true,'tInterval',MDtal.timeInterval_,'prePostFrames',10)
xlim([-5 84])
%% Individual force with average - supplementary
axes('Position', [230/600 530/600 65/600 65/600]); % 
% plotIntensityForce(tracksNAtal(idGroupsTal{1} | idGroupsTal{7}),[],false,false,'UseCurrentAxis',true,'Source',{'ampTotal'},...
%     'plotCohorts',true,'tInterval',MDtal.timeInterval_,'prePostFrames',10,'onlyFirstMode',true)
plotIntensityForce(tracksNAtal(idGroupsTal{1}),[],false,false,'UseCurrentAxis',true,'Source',{'ampTotal'},...
    'plotCohorts',true,'tInterval',MDtal.timeInterval_,'prePostFrames',10,'onlyFirstMode',true)
% plotIntensityForce(tracksG1,[],false,false,'UseCurrentAxis',true,'Source',{'ampTotal'},...
%     'plotCohorts',true,'tInterval',MDtal.timeInterval_,'prePostFrames',10,'onlyFirstMode',true)
xlim([-15 32])
axes('Position', [230/600 430/600 65/600 65/600]); % 
% plotIntensityForce(tracksNAtal(idGroupsTal{1} | idGroupsTal{7}),[],false,false,'UseCurrentAxis',true,'Source',{'forceMag'},...
%     'plotCohorts',true,'tInterval',MDtal.timeInterval_,'prePostFrames',10,'onlyFirstMode',true)
plotIntensityForce(tracksNAtal(idGroupsTal{1}),[],false,false,'UseCurrentAxis',true,'Source',{'forceMag'},...
    'plotCohorts',true,'tInterval',MDtal.timeInterval_,'prePostFrames',10,'onlyFirstMode',true)
xlim([-15 32])
ylim([-15 310])
%% Force-transmitting-normalized - talin
forceTransmittingTal = arrayfun(@(x) ~isempty(x.forceTransmitting) & x.forceTransmitting,tracksNAtal,'UniformOutput',false);
% forceTransmittingTal = arrayfun(@(x) ~isempty(x.forceTransmitting) & x.forceTransmitting,tracksG1,'UniformOutput',false);
idEmptyFTtal=cellfun(@isempty,forceTransmittingTal);
forceTransmittingTal(idEmptyFTtal)={false};
forceTransmittingTal = cell2mat(forceTransmittingTal);
% axes('Position', [0.21 .45 .11 .12]); % only force-transmitting
axes('Position', [530/600 530/600 65/600 65/600]); % only force-transmitting
plotIntensityForce(tracksNAtal((idGroupsTal{1}) & forceTransmittingTal),[],false,false,'UseCurrentAxis',true,...
    'Source',{'ampTotal'},'plotCohorts',true,'tInterval',MDtal.timeInterval_,'prePostFrames',10,'yNormalization',false)
% plotIntensityForce(tracksNAtal((idGroupsTal{1} | idGroupsTal{7}) & forceTransmittingTal),[],false,false,'UseCurrentAxis',true,...
%     'Source',{'ampTotal'},'plotCohorts',true,'tInterval',MDtal.timeInterval_,'prePostFrames',10,'yNormalization',true)
xlim([-5 84])
axes('Position', [530/600 430/600 65/600 65/600]); % only force-transmitting
plotIntensityForce(tracksNAtal((idGroupsTal{1}) & forceTransmittingTal),[],false,false,'UseCurrentAxis',true,...
    'Source',{'forceMag'},'plotCohorts',true,'tInterval',MDtal.timeInterval_,'prePostFrames',10,'yNormalization',false)
% plotIntensityForce(tracksNAtal((idGroupsTal{1} | idGroupsTal{7}) & forceTransmittingTal),[],false,false,'UseCurrentAxis',true,...
%     'Source',{'forceMag'},'plotCohorts',true,'tInterval',MDtal.timeInterval_,'prePostFrames',10,'yNormalization',true)
xlim([-5 84])
ylim([30 75])
%% saving and close
fS3path = '/home2/shan/ownCloud/documents/NascentAdhesionForce/Manuscript/Figures/FigS3heterogeneity';
if ~exist(fS3path,'dir')
    mkdir(fS3path)
end
print('-depsc', '-loose', [fS3path filesep 'talinForceHeterogeneity.eps'])
print('-dtiff', '-loose', '-r300', [fS3path filesep 'talinForceHeterogeneity.tif'])
savefig([fS3path filesep 'talinForceHeterogeneity.fig'])
close
%% Force-transmitting-normalized - talin
forceTransmittingTal = arrayfun(@(x) ~isempty(x.forceTransmitting) & x.forceTransmitting,tracksNAtal,'UniformOutput',false);
% forceTransmittingTal = arrayfun(@(x) ~isempty(x.forceTransmitting) & x.forceTransmitting,tracksG1,'UniformOutput',false);
idEmptyFTtal=cellfun(@isempty,forceTransmittingTal);
forceTransmittingTal(idEmptyFTtal)={false};
forceTransmittingTal = cell2mat(forceTransmittingTal);
% axes('Position', [0.21 .45 .11 .12]); % only force-transmitting
axes('Position', [130/600 530/600 65/600 65/600]); % 
% axes('Position', [530/600 530/600 65/600 65/600]); % only force-transmitting
plotIntensityForce(tracksNAtal((idGroupsTal{1}) & forceTransmittingTal),[],false,false,'UseCurrentAxis',true,...
    'Source',{'ampTotal'},'plotCohorts',true,'tInterval',MDtal.timeInterval_,'prePostFrames',10,'yNormalization',false)
% plotIntensityForce(tracksNAtal((idGroupsTal{1} | idGroupsTal{7}) & forceTransmittingTal),[],false,false,'UseCurrentAxis',true,...
%     'Source',{'ampTotal'},'plotCohorts',true,'tInterval',MDtal.timeInterval_,'prePostFrames',10,'yNormalization',true)
xlim([-10 150])
ylim([251 460])
axes('Position', [130/600 430/600 65/600 65/600]); % only force-transmitting
plotIntensityForce(tracksNAtal((idGroupsTal{1}) & forceTransmittingTal),[],false,false,'UseCurrentAxis',true,...
    'Source',{'forceMag'},'plotCohorts',true,'tInterval',MDtal.timeInterval_,'prePostFrames',10,'yNormalization',false)
% plotIntensityForce(tracksNAtal((idGroupsTal{1} | idGroupsTal{7}) & forceTransmittingTal),[],false,false,'UseCurrentAxis',true,...
%     'Source',{'forceMag'},'plotCohorts',true,'tInterval',MDtal.timeInterval_,'prePostFrames',10,'yNormalization',true)
xlim([-10 150])
ylim([20 100])
%% talin - G2 time series
axes('Position', [230/600 530/600 65/600 65/600]); % 
% axes('Position', [530/600 530/600 65/600 65/600]); % only force-transmitting
plotIntensityForce(tracksNAtal((idGroupsTal{2}) ),[],false,false,'UseCurrentAxis',true,...
    'Source',{'ampTotal'},'plotCohorts',true,'tInterval',MDtal.timeInterval_,'prePostFrames',10,'yNormalization',false)
% plotIntensityForce(tracksNAtal((idGroupsTal{1} | idGroupsTal{7}) & forceTransmittingTal),[],false,false,'UseCurrentAxis',true,...
%     'Source',{'ampTotal'},'plotCohorts',true,'tInterval',MDtal.timeInterval_,'prePostFrames',10,'yNormalization',true)
xlim([-10 150])
ylim([251 460])
axes('Position', [230/600 430/600 65/600 65/600]); % only force-transmitting
plotIntensityForce(tracksNAtal((idGroupsTal{2}) ),[],false,false,'UseCurrentAxis',true,...
    'Source',{'forceMag'},'plotCohorts',true,'tInterval',MDtal.timeInterval_,'prePostFrames',10,'yNormalization',false)
% plotIntensityForce(tracksNAtal((idGroupsTal{1} | idGroupsTal{7}) & forceTransmittingTal),[],false,false,'UseCurrentAxis',true,...
%     'Source',{'forceMag'},'plotCohorts',true,'tInterval',MDtal.timeInterval_,'prePostFrames',10,'yNormalization',true)
xlim([-10 150])
ylim([20 100])

%% intensity and force profiles - vinculin
% axes('Position', [130/600 330/600 65/600 65/600]); % 
% % plotIntensityForce(tracksNAvin(idGroupsVin5{1} | idGroupsVin5{7}),[],false,false,'UseCurrentAxis',true,'Source',{'ampTotal'},'plotCohorts',true,'tInterval',MDvin.timeInterval_,'prePostFrames',10)
% % set(findobj(gca,'Type','Text'),'FontSize',7)
% axes('Position', [130/600 230/600 65/600 65/600]); % 
% plotIntensityForce(tracksNAvin(idGroupsVin5{1} | idGroupsVin5{7}),[],false,false,'UseCurrentAxis',true,'Source',{'forceMag'},'plotCohorts',true,'tInterval',MDvin.timeInterval_,'prePostFrames',10)
%% only F-T - vinculin
forceTransmittingVin = arrayfun(@(x) ~isempty(x.forceTransmitting) & x.forceTransmitting,tracksNAvin,'UniformOutput',false);
idEmptyFTvin=cellfun(@isempty,forceTransmittingVin);
forceTransmittingVin(idEmptyFTvin)={false};
forceTransmittingVin = cell2mat(forceTransmittingVin);
% axes('Position', [230/600 330/600 65/600 65/600]); % only force-transmitting
axes('Position', [130/600 330/600 65/600 65/600]); % 
% plotIntensityForce(tracksNAvin((idGroupsVin5{1} | idGroupsVin5{7}) & forceTransmittingVin),[],false,false,'UseCurrentAxis',true,...
%     'Source',{'ampTotal'},'plotCohorts',true,'tInterval',MDvin.timeInterval_,'prePostFrames',10,'yNormalization',true)
plotIntensityForce(tracksNAvin((idGroupsVin5{1}) & forceTransmittingVin),[],false,false,'UseCurrentAxis',true,...
    'Source',{'ampTotal'},'plotCohorts',true,'tInterval',MDvin.timeInterval_,'prePostFrames',10,'yNormalization',false)
xlim([-5 220])
ylim([230 680])
% axes('Position', [230/600 230/600 65/600 65/600]); % only force-transmitting
axes('Position', [130/600 230/600 65/600 65/600]); % 
% plotIntensityForce(tracksNAvin((idGroupsVin5{1} | idGroupsVin5{7}) & forceTransmittingVin),[],false,false,'UseCurrentAxis',true,...
%     'Source',{'forceMag'},'plotCohorts',true,'tInterval',MDvin.timeInterval_,'prePostFrames',10,'yNormalization',true)
plotIntensityForce(tracksNAvin((idGroupsVin5{1}) & forceTransmittingVin),[],false,false,'UseCurrentAxis',true,...
    'Source',{'forceMag'},'plotCohorts',true,'tInterval',MDvin.timeInterval_,'prePostFrames',10,'yNormalization',false)
xlim([-5 220])
ylim([5 305])
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
%% only F-T - paxillin
% axes('Position', [0.87 .45 .11 .12]); % only force-transmitting
forceTransmittingPax = arrayfun(@(x) ~isempty(x.forceTransmitting) & x.forceTransmitting,tracksNApax,'UniformOutput',false);
idEmptyFTpax=cellfun(@isempty,forceTransmittingPax);
forceTransmittingPax(idEmptyFTpax)={false};
forceTransmittingPax = cell2mat(forceTransmittingPax);
% axes('Position', [230/600 130/600 65/600 65/600]); % only force-transmitting
axes('Position', [130/600 130/600 65/600 65/600]); % 
plotIntensityForce(tracksNApax((idGroupsPax{1}) & forceTransmittingPax),[],false,false,'UseCurrentAxis',true,...
    'Source',{'ampTotal'},'plotCohorts',true,'tInterval',MDpax.timeInterval_,'prePostFrames',10,'yNormalization',false)
% axes('Position', [230/600 30/600 65/600 65/600]); % only force-transmitting
xlim([-5 140])
ylim([300 750])
axes('Position', [130/600 30/600 65/600 65/600]); % 
plotIntensityForce(tracksNApax((idGroupsPax{1}) & forceTransmittingPax),[],false,false,'UseCurrentAxis',true,...
    'Source',{'forceMag'},'plotCohorts',true,'tInterval',MDpax.timeInterval_,'prePostFrames',10,'yNormalization',false)
xlim([-5 130])
ylim([10 150])
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

%% Initial time lag - Loading - now I'll do quantification also for G2, assembly rate and force growth rate
splineParamInit=0.5;
preDetecFactor=1/10; %for paxillin because diffuse signal before point-source-like signal

for ii=1:numel(paxFolder)
%     cur_tlagInitPax = load([paxFolder{ii} '/data/timeInitLagsG1.mat'],'firstIncreseTimeIntAgainstForceAll');
%     cur_tlagInitPax = cur_tlagInitPax.firstIncreseTimeIntAgainstForceAll;
        % try to filter G2 based on starting amplitude, peaks and slopes...
    curTracksNApaxG2 = load([paxFolder{ii} filesep 'data' filesep 'tracksG2.mat'],'tracksG2');
    curTracksNApaxG2 = curTracksNApaxG2.tracksG2; 
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
    
    realTypesFromG1 = zeros(size(curStartingAmpG1));
    for kk=40:numel(curTracksNApaxG1)
        curTracksNApaxG1(kk)=readIntensityFromTracks(curTracksNApaxG1(kk),imgMap,1,'extraLength',30,'movieData',curMD);
        curTracksNApaxG1(kk)=readIntensityFromTracks(curTracksNApaxG1(kk),tMap,2,'movieData',curMD);
        curTracksNApaxG1(kk) = calculateFirstIncreaseTimeTracks(curTracksNApaxG1(kk),splineParamInit,preDetecFactor,tInterval);
        hSum=showSingleAdhesionTrackSummary(curMD,curTracksNApaxG1(kk),imgMap,tMap,kk);
        s=(inputdlg('What class is this NA? (1-9): '));
        realTypesFromG1(kk)=str2num(s{1});
        close(hSum)
    end
    tracksG2fromG1 = curTracksNApaxG1(realTypesFromG1==2);
    tracksG1fromG1 = curTracksNApaxG1(realTypesFromG1==1);
    
    realTypesFromG2 = false(size(realG2pax));
    for kk=find(realG2pax)'
        curTracksNApaxG2(kk)=readIntensityFromTracks(curTracksNApaxG2(kk),imgMap,1,'extraLength',30,'movieData',curMD);
        curTracksNApaxG2(kk)=readIntensityFromTracks(curTracksNApaxG2(kk),tMap,2,'movieData',curMD);
        curTracksNApaxG2(kk) = calculateFirstIncreaseTimeTracks(curTracksNApaxG2(kk),splineParamInit,preDetecFactor,tInterval);
        hSum=showSingleAdhesionTrackSummary(curMD,curTracksNApaxG2(kk),imgMap,tMap,kk);
        s=(inputdlg('What class is this NA? (1-9): '));
        realTypesFromG2(kk)=str2num(s{1});
        close(hSum)
    end
    tracksG2fromG2 = curTracksNApaxG1(realTypesFromG2==2);
    tracksG1fromG2 = curTracksNApaxG1(realTypesFromG2==1);
    
    curTracksNApaxG2 = [tracksG2fromG1; tracksG2fromG2];
    curTracksNApaxG1 = [tracksG1fromG1; tracksG1fromG2];
    
    % Get the time measurement
    tInterval = curMD.timeInterval_;
    [~,curTimeDelayG2] = calculateFirstIncreaseTimeTracks(curTracksNApaxG2,splineParamInit,preDetecFactor,tInterval);
    disp(nanmedian(curTimeDelayG2))
    [~,curTimeDelayG1] = calculateFirstIncreaseTimeTracks(curTracksNApaxG1,splineParamInit,preDetecFactor,tInterval);
    disp(nanmedian(curTimeDelayG1))

    InitTimeLagMaturePax{ii} = curTimeDelayG2;
    disp([paxFolder{ii} ' has been analyzed.'])
    
    tlagInitPax{ii} = curTimeDelayG1; %cur_tlagInitPax;

    cur_tlagPeakPax = load([paxFolder{ii} '/data/timePeaksG1.mat'],'peakTimeIntAgainstForceAll');
    cur_tlagPeakPax = cur_tlagPeakPax.peakTimeIntAgainstForceAll;
    tlagPeakPax{ii} = cur_tlagPeakPax;
    
    % Get the assembly rate and force growth rate
    [~,curEarlyAssmRateG1] = arrayfun(@(x) regression(1:3,x.ampTotal(x.startingFrameExtra:x.startingFrameExtra+2)),curTracksNApaxG1);
    [~,curEarlyAssmRateG2] = arrayfun(@(x) regression(1:5,x.ampTotal(x.startingFrameExtra:x.startingFrameExtra+4)),curTracksNApaxG2);
    earlyAssmRateG1All{ii} = curEarlyAssmRateG1; %cur_tlagInitPax;
    earlyAssmRateG2All{ii} = curEarlyAssmRateG2; %cur_tlagInitPax;
    
    [~,earlyForceRateG1] = arrayfun(@(x) regression(1:3,x.forceMag(x.startingFrameExtra:x.startingFrameExtra+2)),curTracksNApaxG1);
    [~,earlyForceRateG2] = arrayfun(@(x) regression(1:3,x.forceMag(x.startingFrameExtra:x.startingFrameExtra+2)),curTracksNApaxG2);
    earlyForceRateG1All{ii} = earlyForceRateG1; %cur_tlagInitPax;
    earlyForceRateG2All{ii} = earlyForceRateG2; %cur_tlagInitPax;

    tracksG2=curTracksNApaxG2; tracksG1=curTracksNApaxG1;
    save([paxFolder{ii} filesep 'data' filesep 'tracksG2real.mat'],'tracksG2');
    save([paxFolder{ii} filesep 'data' filesep 'tracksG1real.mat'],'tracksG1');
    clear curTracksNApaxG2 curTracksNApaxG1 imgMap tMap
end
totalInitTLagPax=[tlagInitPax{1}; tlagInitPax{2}; tlagInitPax{3}; tlagInitPax{4}];
totalPeakTLagPax=[tlagPeakPax{1}; tlagPeakPax{2}; tlagPeakPax{3}; tlagPeakPax{4}];
totalInitTimeLagMaturePax =[InitTimeLagMaturePax{1}; InitTimeLagMaturePax{2}; 
    InitTimeLagMaturePax{3}; InitTimeLagMaturePax{4}];
totalEarlyAssmRateG1Pax=[earlyAssmRateG1All{1}; earlyAssmRateG1All{2}; earlyAssmRateG1All{3}; earlyAssmRateG1All{4}];
totalEarlyAssmRateG2Pax=[earlyAssmRateG2All{1}; earlyAssmRateG2All{2}; earlyAssmRateG2All{3}; earlyAssmRateG2All{4}];
totalEarlyForceRateG1Pax=[earlyForceRateG1All{1}; earlyForceRateG1All{2}; earlyForceRateG1All{3}; earlyForceRateG1All{4}];
totalEarlyForceRateG2Pax=[earlyForceRateG2All{1}; earlyForceRateG2All{2}; earlyForceRateG2All{3}; earlyForceRateG2All{4}];
clear earlyAssmRateG2All earlyAssmRateG1All earlyForceRateG1All earlyForceRateG2All
%% for vinculin 
preDetecFactor=1/10; %for vinculin 
for ii=1:numel(vinFolder)
%     cur_tlagInitVin = load([vinFolder{ii} '/data/timeInitLagsG1.mat'],'firstIncreseTimeIntAgainstForceAll');
%     cur_tlagInitVin = cur_tlagInitVin.firstIncreseTimeIntAgainstForceAll;
%     tlagInitVin{ii} = cur_tlagInitVin;
    % try to filter G2 based on starting amplitude, peaks and slopes...
    curTracksNAvinG2 = load([vinFolder{ii} filesep 'data' filesep 'tracksG2.mat'],'tracksG2');
    curTracksNAvinG2 = curTracksNAvinG2.tracksG2; 
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
    
    realTypesFromG1 = zeros(size(curStartingAmpG1));
    for kk=1:numel(curTracksNAvinG1)
        curTracksNAvinG1(kk)=readIntensityFromTracks(curTracksNAvinG1(kk),imgMap,1,'extraLength',30,'movieData',curMD);
        curTracksNAvinG1(kk)=readIntensityFromTracks(curTracksNAvinG1(kk),tMap,2,'movieData',curMD);
        curTracksNAvinG1(kk) = calculateFirstIncreaseTimeTracks(curTracksNAvinG1(kk),splineParamInit,preDetecFactor,tInterval);
        hSum=showSingleAdhesionTrackSummary(curMD,curTracksNAvinG1(kk),imgMap,tMap,kk);
        s=(inputdlg('What class is this NA? (1-9): '));
        try
            realTypesFromG1(kk)=str2num(s{1});
        catch
            realTypesFromG1(kk)=3;
        end
        close(hSum)
    end
    tracksG2fromG1 = curTracksNAvinG1(realTypesFromG1==2);
    tracksG1fromG1 = curTracksNAvinG1(realTypesFromG1==1);
    
    realTypesFromG2 = zeros(size(realG2vin));
    for kk=find(realG2vin)'
        curTracksNAvinG2(kk)=readIntensityFromTracks(curTracksNAvinG2(kk),imgMap,1,'extraLength',30,'movieData',curMD);
        curTracksNAvinG2(kk)=readIntensityFromTracks(curTracksNAvinG2(kk),tMap,2,'movieData',curMD);
        curTracksNAvinG2(kk) = calculateFirstIncreaseTimeTracks(curTracksNAvinG2(kk),splineParamInit,preDetecFactor,tInterval);
        hSum=showSingleAdhesionTrackSummary(curMD,curTracksNAvinG2(kk),imgMap,tMap,kk);
        s=(inputdlg('What class is this NA? (1-9): '));
        try
            realTypesFromG2(kk)=str2num(s{1});
        catch
            realTypesFromG2(kk)=3;
        end
        close(hSum)
    end
    tracksG2fromG2 = curTracksNAvinG1(realTypesFromG2==2);
    tracksG1fromG2 = curTracksNAvinG1(realTypesFromG2==1);
    
    curTracksNAvinG2 = [tracksG2fromG1; tracksG2fromG2];
    curTracksNAvinG1 = [tracksG1fromG1; tracksG1fromG2];
    
    % Get the time measurement
    tInterval = curMD.timeInterval_;
    [~,curTimeDelayG2] = calculateFirstIncreaseTimeTracks(curTracksNAvinG2,splineParamInit,preDetecFactor,tInterval);
    disp(nanmedian(curTimeDelayG2))
    [~,curTimeDelayG1] = calculateFirstIncreaseTimeTracks(curTracksNAvinG1,splineParamInit,preDetecFactor,tInterval);
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
    tracksG2=curTracksNAvinG2; tracksG1=curTracksNAvinG1;
    save([vinFolder{ii} filesep 'data' filesep 'tracksG2real.mat'],'tracksG2');
    save([vinFolder{ii} filesep 'data' filesep 'tracksG1real.mat'],'tracksG1');
    clear curTracksNAvinG2 curTracksNAvinG1 imgMap tMap
end
totalInitTLagVin=[tlagInitVin{1}; tlagInitVin{2}; tlagInitVin{3}; tlagInitVin{4}; tlagInitVin{5}; tlagInitVin{6}];
totalPeakTLagVin=[tlagPeakVin{1}; tlagPeakVin{2}; tlagPeakVin{3}; tlagPeakVin{4}; tlagPeakVin{5}; tlagPeakVin{6}];
totalInitTimeLagMatureVin =[InitTimeLagMatureVin{1}; InitTimeLagMatureVin{2}; 
    InitTimeLagMatureVin{3}; InitTimeLagMatureVin{4}; InitTimeLagMatureVin{5}; InitTimeLagMatureVin{6}];
totalEarlyAssmRateG1Vin=[earlyAssmRateG1All{1}; earlyAssmRateG1All{2}; earlyAssmRateG1All{3}; earlyAssmRateG1All{4}; earlyAssmRateG1All{5}; earlyAssmRateG1All{6}];
totalEarlyAssmRateG2Vin=[earlyAssmRateG2All{1}; earlyAssmRateG2All{2}; earlyAssmRateG2All{3}; earlyAssmRateG2All{4}; earlyAssmRateG2All{5}; earlyAssmRateG2All{6}];
totalEarlyForceRateG1Vin=[earlyForceRateG1All{1}; earlyForceRateG1All{2}; earlyForceRateG1All{3}; earlyForceRateG1All{4}; earlyForceRateG1All{5}; earlyForceRateG1All{6}];
totalEarlyForceRateG2Vin=[earlyForceRateG2All{1}; earlyForceRateG2All{2}; earlyForceRateG2All{3}; earlyForceRateG2All{4}; earlyForceRateG2All{5}; earlyForceRateG2All{6}];
clear earlyAssmRateG2All earlyAssmRateG1All earlyForceRateG1All earlyForceRateG2All

%% for talin
totalInitTLagTal=[];
totalPeakTLagTal=[];
splineParamInit=0.95;
preDetecFactor=1.5/10; 
for ii=1:1
%     cur_tlagInitTal = load([talFolder{ii} '/data/timeInitLagsG1.mat'],'firstIncreseTimeIntAgainstForceAll');
%     cur_tlagInitTal = cur_tlagInitTal.firstIncreseTimeIntAgainstForceAll;
%     if ii>1
%         cur_tlagInitTal=cur_tlagInitTal*2;
%     end
    % try to filter G2 based on starting amplitude, peaks and slopes...
    curTracksNAtalG2 = load([talFolder{ii} filesep 'data' filesep 'tracksG2.mat'],'tracksG2');
    curTracksNAtalG2 = curTracksNAtalG2.tracksG2; 
    curForceTransmittingG2 = arrayfun(@(x) ~isempty(x.firstIncreseTimeIntAgainstForce), curTracksNAtalG2);
    curStartingAmpG2 = arrayfun(@(x) x.ampTotal(x.startingFrameExtra), curTracksNAtalG2);

    curTracksNAtalG1 = load([talFolder{ii} filesep 'data' filesep 'tracksG1.mat'],'tracksG1');
    curTracksNAtalG1 = curTracksNAtalG1.tracksG1; 
    curTimeDelay = arrayfun(@(x) x.firstIncreseTimeIntAgainstForce, curTracksNAtalG1,'UniformOutput',false);
    curStartingAmpG1 = arrayfun(@(x) x.ampTotal(x.startingFrameExtra), curTracksNAtalG1);
    thresAmp = mean(curStartingAmpG1)+3*std(curStartingAmpG1);

    % filtering for G2
    stateFA = arrayfun(@(x) any(cellfun(@(y) strcmp(y,'FC') || strcmp(y,'FA'), x.state)), curTracksNAtalG2);
    realG2tal = curForceTransmittingG2 & curStartingAmpG2<thresAmp & stateFA;
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
    
    realTypesFromG1 = zeros(size(curStartingAmpG1));
    for kk=1:numel(curTracksNAtalG1)
        curTracksNAtalG1(kk)=readIntensityFromTracks(curTracksNAtalG1(kk),imgMap,1,'extraLength',30,'movieData',curMD);
        curTracksNAtalG1(kk)=readIntensityFromTracks(curTracksNAtalG1(kk),tMap,2,'movieData',curMD);
        curTracksNAtalG1(kk) = calculateFirstIncreaseTimeTracks(curTracksNAtalG1(kk),splineParamInit,preDetecFactor,tInterval);
        hSum=showSingleAdhesionTrackSummary(curMD,curTracksNAtalG1(kk),imgMap,tMap,kk);
        s=(inputdlg('What class is this NA? (1-9): '));
        try
            realTypesFromG1(kk)=str2num(s{1});
        catch
            realTypesFromG1(kk)=3;
        end
        close(hSum)
    end
    tracksG2fromG1 = curTracksNAtalG1(realTypesFromG1==2);
    tracksG1fromG1 = curTracksNAtalG1(realTypesFromG1==1);
    
    realTypesFromG2 = zeros(size(realG2tal));
    for kk=find(realG2tal)'
        curTracksNAtalG2(kk)=readIntensityFromTracks(curTracksNAtalG2(kk),imgMap,1,'extraLength',30,'movieData',curMD);
        curTracksNAtalG2(kk)=readIntensityFromTracks(curTracksNAtalG2(kk),tMap,2,'movieData',curMD);
        curTracksNAtalG2(kk) = calculateFirstIncreaseTimeTracks(curTracksNAtalG2(kk),splineParamInit,preDetecFactor,tInterval);
        hSum=showSingleAdhesionTrackSummary(curMD,curTracksNAtalG2(kk),imgMap,tMap,kk);
        s=(inputdlg('What class is this NA? (1-9): '));
        try
            realTypesFromG2(kk)=str2num(s{1});
        catch
            realTypesFromG2(kk)=3;
        end
        close(hSum)
    end
    tracksG2fromG2 = curTracksNAtalG1(realTypesFromG2==2);
    tracksG1fromG2 = curTracksNAtalG1(realTypesFromG2==1);
    
    curTracksNAtalG2 = [tracksG2fromG1; tracksG2fromG2];
    curTracksNAtalG1 = [tracksG1fromG1; tracksG1fromG2];
    
    % Get the time measurement
    tInterval = curMD.timeInterval_;
    [~,curTimeDelayG1] = calculateFirstIncreaseTimeTracks(curTracksNAtalG1,splineParamInit,preDetecFactor,tInterval);
    disp(nanmedian(curTimeDelayG1))
    [~,curTimeDelayG2] = calculateFirstIncreaseTimeTracks(curTracksNAtalG2,splineParamInit,preDetecFactor,tInterval);
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
    [~,curEarlyAssmRateG2] = arrayfun(@(x) regression(1:5,x.ampTotal(x.startingFrameExtra:x.startingFrameExtra+4)),curTracksNAtalG2);
    earlyAssmRateG1All{ii} = curEarlyAssmRateG1; %cur_tlagInitPax;
    earlyAssmRateG2All{ii} = curEarlyAssmRateG2; %cur_tlagInitPax;
    
    [~,earlyForceRateG1] = arrayfun(@(x) regression(1:3,x.forceMag(x.startingFrameExtra:x.startingFrameExtra+2)),curTracksNAtalG1);
    [~,earlyForceRateG2] = arrayfun(@(x) regression(1:3,x.forceMag(x.startingFrameExtra:x.startingFrameExtra+2)),curTracksNAtalG2);
    earlyForceRateG1All{ii} = earlyForceRateG1; %cur_tlagInitPax;
    earlyForceRateG2All{ii} = earlyForceRateG2; %cur_tlagInitPax;
    tracksG2=curTracksNAtalG2; tracksG1=curTracksNAtalG1;
    save([talFolder{ii} filesep 'data' filesep 'tracksG2real.mat'],'tracksG2');
    save([talFolder{ii} filesep 'data' filesep 'tracksG1real.mat'],'tracksG1');
%     clear curTracksNAtalG2
end
totalInitTLagTal =[tlagInitTal{1}];
totalInitTimeLagMatureTal =[InitTimeLagMatureTal{1}];
% totalInitTLagTal=cell2mat(tlagInitTal);%[tlagInitTal{1}];
% totalPeakTLagTal=tlagPeakTal;%[tlagPeakTal{1}];

totalEarlyAssmRateG1Tal=[earlyAssmRateG1All{1}];
totalEarlyAssmRateG2Tal=[earlyAssmRateG2All{1}];
totalEarlyForceRateG1Tal=[earlyForceRateG1All{1}];
totalEarlyForceRateG2Tal=[earlyForceRateG2All{1}];
%% filtering with high and low limit (-60 and 60)
lowTLag=-80; hiTLag=60;
totalInitTLagTal2 = totalInitTLagTal(totalInitTLagTal>lowTLag & totalInitTLagTal<hiTLag);
totalInitTLagVin2 = totalInitTLagVin(totalInitTLagVin>lowTLag & totalInitTLagVin<hiTLag);
totalInitTLagPax2 = totalInitTLagPax(totalInitTLagPax>lowTLag & totalInitTLagPax<hiTLag);

%% Looking at each cell
medPaxInit=median(totalInitTLagPax2);
medVinInit=median(totalInitTLagVin2);
medTalInit=median(totalInitTLagTal2);

disp([cellfun(@(x) nanmedian(x),tlagInitPax); cellfun(@(x) length(x),tlagInitPax)])
disp([cellfun(@(x) nanmedian(x),tlagInitVin); cellfun(@(x) length(x),tlagInitVin)])
disp([cellfun(@(x) nanmedian(x),tlagInitTal); cellfun(@(x) length(x),tlagInitTal)])

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

%% statistical testing for inital lag
% normality test
[hNormalPax,pNormalPax] = kstest(totalInitTLagPax2);
% hNormalPax = 1
% pNormalPax = 7.8863e-242
% So it's not normal! 
[pPaxVsVin,hPaxVsVin] = ranksum(totalInitTLagPax2,totalInitTLagVin2);
[pPaxVsTal,hPaxVsTal] = ranksum(totalInitTLagPax2,totalInitTLagTal2);
[pVinVsTal,hVinVsTal] = ranksum(totalInitTLagVin2,totalInitTLagTal2);
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
%% Actual Fig 3
% close all
% % positioning
% pos = get(0, 'DefaultFigurePosition');
% pos(3:4) = [300 300];
% figure('Position', pos, 'PaperPositionMode', 'auto','Color','w');
% % axes('Position', [0.05 .55 .90 .4]);
% subplot(2,1,1)
% box-plot
% making it to matrix
[lengthLongest,iLongest]=max([length(totalInitTLagPax2),length(totalInitTLagVin2),length(totalInitTLagTal2)]);
% it was vinculin that was longest
matrixPaxTalVin = NaN(lengthLongest,3);
matrixPaxTalVin(1:length(totalInitTLagPax2),1) = totalInitTLagPax2;
matrixPaxTalVin(1:length(totalInitTLagVin2),2) = totalInitTLagVin2;
matrixPaxTalVin(1:length(totalInitTLagTal2),3) = totalInitTLagTal2;
boxWidth=0.5;

axes('Position',[380/600,325/600,80/600,65/600])
boxplot(matrixPaxTalVin,'orientation','horizontal','whisker',0.5,'notch','on',...
    'labels',{['Pax (N=' num2str(length(totalInitTLagPax2)) ')'],...
    ['Vin (N=' num2str(length(totalInitTLagVin2)) ')'],...
    ['Tal (N=' num2str(length(totalInitTLagTal2)) ')']},...
    'symbol','','widths',boxWidth,'jitter',1,'colors','k')
txt = findobj(gca,'Type','text');
set(txt(3:end),'VerticalAlignment', 'Middle');
set(findobj(gca,'LineStyle','--'),'LineStyle','-')
set(findobj(gca,'tag','Median'),'LineWidth',2)
xlim([-35 30])
xlabel('Time lag in initial rise (s)')
title('non-maturing adhesions (G1)')

% p-values
% findobj(gca,'tag','Median')
% findobj(gca,'tag','Outlier')
% get(findobj(gca,'Group','Pax'),'OutlierValue')
% [x y]=ginput(1);
text(14,1.5,['p=' num2str(pPaxVsVin,2)]) %pax vs vin
text(14,2.5,['p=' num2str(pVinVsTal,2)]) %tal vs vin
text(14,2,['p=' num2str(pPaxVsTal,2)]) %pax vs tal

% median values
text(medPaxInit,1.5,num2str(medPaxInit,'%10.0f'),'HorizontalAlignment','center') %pax vs vin
text(medVinInit,2.5,num2str(medVinInit,'%10.0f'),'HorizontalAlignment','center') %tal vs vin
text(medTalInit,3.5,num2str(medTalInit,'%10.0f'),'HorizontalAlignment','center') %pax vs tal

set(findobj(gca,'Type','Text'),'FontSize',7)
set(gca,'FontSize',7)

% title('Time lag in initial rises')
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

%% Maturing adhesions - t_init - plotting
[lengthLongest,iLongest]=max([length(totalInitTimeLagMaturePax2),length(totalInitTimeLagMatureVin2),length(totalInitTimeLagMatureTal2)]);
% it was vinculin that was longest
matrixInitTimeLagG2 = NaN(lengthLongest,3);
matrixInitTimeLagG2(1:length(totalInitTimeLagMaturePax2),1) = totalInitTimeLagMaturePax2;
matrixInitTimeLagG2(1:length(totalInitTimeLagMatureVin2),2) = totalInitTimeLagMatureVin2;
matrixInitTimeLagG2(1:length(totalInitTimeLagMatureTal2),3) = totalInitTimeLagMatureTal2;
% axes('Position',[0.08,0.45,0.40,0.18])
axes('Position', [380/600 225/600 80/600 65/600]); % 
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
