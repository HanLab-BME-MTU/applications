function [ output_args ] = coritcalCalculations(projList,edgeDist,distCutoff,jitterParam, makeMovie,subRoiFilename,filterByAngle)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% distCutoff = displacement cut off for high frequency motion 
% jitterParm: microtubules that  move less than this absolute distance will
% be considered jitter not true lateral movement % default = 0.7 um 
% make movie to make movie for each 
% filterByAngle (might not want to include information for the whole pole) 

for iProj = 1:numel(projList)
  %  anDir = projList(iProj).anDir; 
%     anDir = strrep(anDir,'Mijung','Mijung_RPE_MitoticSegmented');
   % imDir = projList(iProj).imDir; 
%     imDir = strrep(imDir,'Mijung','Mijung_RPE_MitoticSegmented'); 
anDir = projList{iProj}; 
    x = getFilenameBody(anDir);
    imDir = [x filesep 'images']; 
    [up2,con,cellNum] = getFilenameBody(x);
    [up3,~] = getFilenameBody(up2); 
    
    % load mask internal
    maskExternal = logical(imread([anDir filesep 'roiMaskSeg.tif']));
    maskInternal = logical(imread([anDir filesep 'roiMaskKinetochore.tif']));
    
    [listOfImages] = searchFiles('.tif',[],imDir,0);
    imageName1 = [char(listOfImages(1,2)) filesep char(listOfImages(1,1))];
    img =  double(imread(imageName1));
    
    if isempty(maskInternal)
        maskInternal = roipoly;
    end
    
    mask = maskExternal-maskInternal;
    % get long centroid
    figure('Visible','off')
    imshow(img,[]); 
   props = regionprops(maskInternal,'centroid','orientation');
     centroid = props.Centroid;
%     angle =- props.Orientation;
% 
% for iPt = 1:2 
% h=msgbox('click the centrosome','help');
% uiwait(h);
% h=impoint;
% position = wait(h);
% ptMask = zeros(size(img)); 
% idx = sub2ind(size(img),round(position(2)),round(position(1))); 
% point(iPt,2) = position(2);
% point(iPt,1) = position(1); 
% ptMask(idx) = 1; 
% poles(:,:,iPt) = ptMask; 
% 
% 
% end 
% 
% close gcf
s = load([anDir filesep 'poleInfo.mat']); 
 %s = load([anDir filesep 'subRois_Poles' filesep 'poles.mat']); 
  %poles = s.poles; 
 
%pole1= poles(:,:,1); 
 %pole2 = poles(:,:,2); 
 %[y1,x1] = ind2sub(size(img),find(pole1==1)); 
 %[y2,x2] = ind2sub(size(img),find(pole2==1));
 %point = [x1,y1;x2,y2]; 

  point = s.poleInfo.coords; 
  poles = s.poleInfo.poleMasks;
    mAxis = (point(1,2)-point(2,2))/(point(1,1)-point(2,1));
    angle = atand(mAxis);
    bAxis = point(1,2) - mAxis*point(1,1);
    
   mAxisPer = -1/mAxis;
   anglePer = atand(mAxisPer); 
   bAxisPer = centroid(2)-mAxisPer*centroid(1); 
   
   % full perpendicular ROIS for Assymetry 
     [imL, imW ] = size(mask);
   x = [1:imW]; 
    y = mAxis +bAxis; 
     
    [~,yAll]=meshgrid(1:imW,1:imL);
    yLineAxis=repmat(mAxis.*[1:imW]+bAxis,[imL,1]);
    yPer = repmat(mAxisPer.*[1:imW]+bAxisPer,[imL,1]);
    
     r11 = yAll<=yPer;
     r12 = yAll>yPer;
     
   
    
     % this will be used to find intersection and test for long axis label
    x = [1:imW] ; 
    y = mAxis.*x+bAxis; 
    x = x(y<imL); 
    y = y(y<imL);
    x = x(y>0); 
    y= y(y>0); 
    pixelTest = sub2ind(size(img),ceil(y),x); 
lineMask = zeros(size(mask)); 
lineMask(pixelTest) = 1; 
     
% sort the poles to correspond to one pole or the other. 
% pole test 
wholeMask = imfill(mask); 
test1 = r11& wholeMask; 
test2 = r12& wholeMask;
idxtest1 =  find(test1==1); 
idxpole1 = find(poles(:,:,1)==1); 
overlap= intersect(idxtest1,idxpole1); 
% switch the regions
if isempty(overlap)
   r11 = yAll>yPer; 
   r12 = yAll<=yPer; 
end 
%% Asymmetry SubRois    
%     m1 =  tand(anglePlusMinus + angle);
%     m2 = tand(angle -anglePlusMinus);
%     b1 = centroid(2) - m1*centroid(1);
%     b2 = centroid(2) - m2*centroid(1);
%     
%    
%      
 
%     
%      x = [1:imW]; 
%      y = mAxis +bAxis; 
%     [~,yAll]=meshgrid(1:imW,1:imL);
%     
%     yLineAxis=repmat(mAxis.*[1:imW]+bAxis,[imL,1]);
% %     y1 = repmat(m1.*[1:imW]+b1,[imL,1]);
% %     y2 = repmat(m2.*[1:imW]+b2,[imL,1]);
%     line = yAll==yLineAxis; 
%     
%     r11 = yAll<=y1;
%     r12 = yAll>y1;
%     
%     r21 = yAll<= y2;
%     r22 = yAll>y2;
%     
%     roiSetTest(:,:,1)= r11 & mask & r22;
%     roiSetTest(:,:,2) = r12 & mask & r21;
%     roiSetTest(:,:,3) = r11 & mask & r21;
%     roiSetTest(:,:,4)  = r12 & mask & r22;
%     
%     %First test for long axis label; 
    
    %     
%     countLong = 1; 
%     countShort = 3; 
%     for i = 1:4 
%         CCTest = bwconncomp(roiSetTest(:,:,i)); 
%       First  pixels = CCTest.PixelIdxList{1}; 
%         if ~isempty(intersect(pixelTest,pixels)); 
%            
%             roiSet(:,:,countLong) = roiSetTest(:,:,i); 
%         countLong = countLong+1; 
%         else 
%             roiSet(:,:,countShort) = roiSetTest(:,:,i); 
%             countShort = countShort+1; 
%         end 
%      end 
%     figure('Visible','on');
%     
%     imshow(img,[])
%     hold on
%     plot(x,y,'b');
%      colors = ['r','b','g','y'];
%     for iRoi = 1:4
%         
%         currentSubRoiDir = [anDir filesep 'subRois_Poles' filesep 'sub_' num2str(iRoi)];
%         if isdir(currentSubRoiDir)
%             rmdir(currentSubRoiDir,'s');
%         end
%         mkdir([currentSubRoiDir filesep 'meta'])
%         mkdir([currentSubRoiDir filesep 'feat'])
%         movieInfo1  = [anDir filesep 'feat' filesep 'movieInfo.mat'];
%         movieInfo2 = [currentSubRoiDir filesep 'feat' filesep 'movieInfo.mat'];
%         copyfile(movieInfo1,movieInfo2);
%         roiMask = roiSet(:,:,iRoi);
%         imwrite(roiMask,[currentSubRoiDir filesep 'roiMask.tif']);
%         
%         [y1,x1]=ind2sub([imL,imW],find(roiSet(:,:,iRoi),1)); % first pixel on boundary
%         roiYX = bwtraceboundary(roiSet(:,:,iRoi),[y1,x1],'N');
%         save([currentSubRoiDir filesep 'roiYX.mat'],'roiYX');
%         
%         
%         plot(roiYX(:,2),roiYX(:,1),colors(iRoi));
%         
%         % make weighted mask using distance transform to find position
%         % where text should go
%         weightedRoi=bwdist(~roiSet(:,:,iRoi));
%         [r,c]=find(weightedRoi==max(weightedRoi(:)));
%         text(c(1),r(1),num2str(iRoi),'color','r','fontsize',14);
%      
%         
%    %   plusTipSubRoiExtractTracksFORPACKAGE(currentSubRoiDir,'turnFiguresOn',0);
%         
%     end
%     
%    
%        spy(poles(:,:,1)|poles(:,:,2),25,'y'); 
%        scatter(centroid(1),centroid(2),10,'c'); 
%     if ~isdir([saveDir filesep 'CollectedSubRois' filesep 'Poles']);
%         mkdir([saveDir filesep 'CollectedSubRois' filesep 'Poles']);
%        
%     end
%     saveas(gcf,[saveDir filesep 'CollectedSubRois' filesep 'Poles' filesep ['subRois_Cell_' num2str(cellNum) '.tif']]);
%      saveas(gcf,[anDir filesep 'subRois_Poles' filesep 'subRoiPoles.eps'],'psc2'); 
%      save([anDir filesep 'subRois_Poles' filesep 'poles.mat'],'poles'); 
%      
%      
%     close(gcf)
% %   name = ['HSET_M10P_Db' num2str(cellNum) '_w1DIC_t1.TIF']; 
% %   name = strrep(name,'b0','b'); 
% %    DIC =  double(imread([up2 filesep 'HM_Db_DIC_other_files' filesep 'images' filesep name])); 
%    figure('visible','off');
% %    imshow(DIC,[]); 
%    hold on 
%    for iRoi = 1:4
%          [y1,x1]=ind2sub([imL,imW],find(roiSet(:,:,iRoi),1)); % first pixel on boundary
%         roiYX = bwtraceboundary(roiSet(:,:,iRoi),[y1,x1],'N');
%          plot(roiYX(:,2),roiYX(:,1),colors(iRoi));
%         
%         % make weighted mask using distance transform to find position
%         % where text should go
%         weightedRoi=bwdist(~roiSet(:,:,iRoi));
%         [r,c]=find(weightedRoi==max(weightedRoi(:)));
%         text(c(1),r(1),num2str(iRoi),'color','r','fontsize',14);
%    end
%    if ~isdir([saveDir filesep 'CollectedSubRois' filesep 'PolesDIC'])
%        mkdir([saveDir filesep 'CollectedSubRois' filesep 'PolesDIC'])
%    end 
%   
%    saveas(gcf,[saveDir filesep 'CollectedSubRois' filesep 'PolesDIC' filesep ['subRois_Cell_DIC' num2str(cellNum) '.tif']]); 
%         
%    close(gcf)
%         
%     clear roiSet roiSetTest poles; 
%% SET-UP CORTICAL SUBREGION FOLDERS

%%%%% MAKE SUBROIS %%%%%
% find initial distance transform from cell edge (1 um) 
subRoiEdgeMask = imfill(mask);
distTrans = bwdist(~subRoiEdgeMask); 

forIntersectAll = subRoiEdgeMask;
forIntersectAll(distTrans>1) = 0; 


distTrans = distTrans.*.108; %  convert to microns
subRoiEdgeMask(distTrans>edgeDist )=0;



roiSet(:,:,1) = subRoiEdgeMask & r11; 
roiSet(:,:,2) = subRoiEdgeMask & r12;
forIntersect(:,:,1) = forIntersectAll & r11; 
forIntersect(:,:,2) = forIntersectAll & r12; 


for iSub = 1:2
% set up new subRoiDirectory 
 subRoiEdgeDir = [anDir filesep subRoiFilename filesep 'sub_' num2str(iSub)];
 if ~isdir(subRoiEdgeDir) 
     mkdir(subRoiEdgeDir)
 end 
  mkdir([subRoiEdgeDir filesep 'meta'])
  mkdir([subRoiEdgeDir filesep 'feat'])
         movieInfo1  = [anDir filesep 'feat' filesep 'movieInfo.mat'];
         movieInfo2 = [subRoiEdgeDir filesep 'feat' filesep 'movieInfo.mat'];
         copyfile(movieInfo1,movieInfo2);
         
 imwrite(roiSet(:,:,iSub),[subRoiEdgeDir filesep 'roiMask.tif']); 
end % iSub 

%%% REGION-OVER-VIEW SANITY PLOT %%%
% search for DIC file using list images 
[DICs] = searchFiles('DIC',[],anDir,0);

imgDIC = double(imread([char(DICs(1,2)) filesep char(DICs(1,1))])); 
%

toPlot(:,:,1) = img; 
toPlot(:,:,2) = imgDIC;
names{1} = 'EB_Comet'; 
names{2} = 'DIC'; 
 for iPlot = 1:2 
     figure('Visible','off'); 
 imshow(toPlot(:,:,iPlot),[])
    hold on
    plot(x,y,'r');
    plotScaleBar(100,5,'Label','10um'); 
     colors = ['b','g'];
     for iSub = 1:2
         roiYX = bwboundaries(roiSet(:,:,iSub));
         roiY = roiYX{1}(:,1);
         roiX = roiYX{1}(:,2);
         plot(roiX,roiY,colors(iSub));
         weightedRoi=bwdist(~roiSet(:,:,iSub));
         [r,c]=find(weightedRoi==max(weightedRoi(:)));
         text(c(1),r(1),num2str(iSub),'color','r','fontsize',14);
         
         %calculate p (distance from pole to cortex)
         pixBound = find(forIntersect(:,:,iSub) ==1);     
         edgePt = intersect(pixBound,pixelTest);
         while isempty(edgePt)
             dil = zeros(size(mask));
             dil(pixBound) = 1;
             dil =  imdilate(dil,strel('disk',1));
             pixBound = find(dil==1);
             edgePt = intersect(pixBound,pixelTest);
         end
         if length(edgePt) >1
             edgePt = edgePt(1); % just take the first
         end
         [yIntersect(iSub), xIntersect(iSub)] = ind2sub(size(img),edgePt);
         scatter(xIntersect(iSub),yIntersect(iSub),100,colors(iSub),'filled');
         d = sqrt((point(iSub,1)-xIntersect(iSub)).^2 + (point(iSub,2)-yIntersect(iSub)).^2);
         p(iSub) = d.*0.108;
         text(xIntersect(iSub)+1,yIntersect(iSub)+1,['P = ' num2str(p(iSub),2) ' um'],'color','y');
        
     end % for iSub
     % plot the poles 
     for iPole= 1:2 
     spy(poles(:,:,iPole),colors(iPole),30) 
     end 
      spindleL = sqrt((point(1,1)-point(2,1)).^2 + (point(1,2)-point(2,2)).^2);
      spindleL = spindleL*0.108;
         text(centroid(1)+1,centroid(2)+1,['L = ' num2str(spindleL,2),' um'],'color','y');
      
    
     % save two places 
     saveas(gcf,[anDir filesep subRoiFilename filesep names{iPlot} '.eps'],'psc2'); 
     collectedDir = [up3 filesep 'CollectedSubRois' filesep 'Cortical' filesep 'RegionOverView' filesep  con]; 
 
     
     if ~isdir(collectedDir)
         mkdir(collectedDir)
     end 
saveas(gcf,[collectedDir filesep 'subRois_Cortical_' names{iPlot} '_' con '_' cellNum '.tif']);
close(gcf) 
% save DIC as well 

 end % iPlot
 % save the pole masks for later use likely better way to store but who
 % cares
save([anDir filesep subRoiFilename  filesep 'poles.mat'],'poles'); 

%% EXTRACT INFORMATION FROM CORTICAL REGIONS AND CLASSIFICATION
%         
%    close(gcf)
%         
 %   clear roiSet roiSetTest poles; 
% make subroiDir get xCoord in all tracks and yCoord in all tracks 
for iSub = 1:2
   subRoiEdgeDir = [anDir filesep subRoiFilename filesep 'sub_' num2str(iSub)];
   
projData = plusTipSubRoiExtractTracksFORPACKAGE(subRoiEdgeDir,'fraction',0.01,'extractType','boundCrossIn','turnFiguresOn',0);
projData.stats.spindleLength = spindleL;
%%% CLASSIFICATION OF EVENTS %%%
vel =projData.speedInMicPerMinAllTracks;

if ~isempty(vel)
vel = vel(~projData.startOrEnd); % remove those tracks at the start or the end of the movie
life = projData.insideSecAllTracks; 
life = life(~projData.startOrEnd); 
dist = vel.*life/60; 
% filter tracks based on dist and velocity
 %velCutOff = 1; % in um/sec maintain consistency with NCB paper
idxLat = find(dist>distCutoff); 
idxPer = find(dist<distCutoff);  
 
projData.stats.staticDwell = nanmean(life(idxPer));
projData.stats.lowDispVelMean = nanmean(vel(idxPer)); 


projData.stats.freqCortical = length(vel)/(projData.nFrames)/projData.secPerFrame*60; % convertToNumberPerMinutePerPole
param2 = ['freqCorticalDispGreater' num2str(distCutoff) 'um_All'];
param3 = ['percentDispGreater' num2str(distCutoff) 'um_All'];
projData.stats.(param2) = length(idxLat)/projData.nFrames/projData.secPerFrame*60; 
projData.stats.(param3) = length(idxLat)/length(vel)*100; 





% latVel = latVel/60;  % convert to um/sec 
% idxFast = find(latVel>velCutOff); 

% projData.stats.freqCorticalDispFast = length(idxFast); 
% projData.stats.freqCorticalDispLong = sum(dist>2.5); 

% load the old astral growth data to normalize
s = load([anDir filesep 'subRois_AstralOnly' filesep 'sub_1' filesep  'meta' filesep 'projData.mat']); 
projDataOld = s.projData; 
astralGrowth = projDataOld.stats.growth_speed_mean_INSIDE_REGION ; 
projData.stats.lowDispVelMeanNorm = projData.stats.lowDispVelMean/astralGrowth; 

    
%% start plotting 
xCoordsIn = projData.xCoordInAllTracks; 
yCoordsIn = projData.yCoordInAllTracks; 
xCoordsOut = projData.xCoordOutAllTracks;
yCoordsOut = projData.yCoordOutAllTracks; 

% since lifetime calcs need to take out those that start or end in the
% first/last frame 
xCoordsIn = xCoordsIn(~projData.startOrEnd,:); 
yCoordsIn = yCoordsIn(~projData.startOrEnd,:); 
xCoordsOut = xCoordsOut(~projData.startOrEnd,:); 
yCoordsOut = yCoordsOut(~projData.startOrEnd,:); 


%%% SANITY PLOT: CLASSIFICATION OF HI AND LOW DISPLACEMENT CORTICAL TRACKS%%% 
figure('Visible','off'); 

imshow(toPlot(:,:,2),[]) % plot the DIC
hold on 
 plot(x,y,'c'); % plot the pole axis 

%     
% for iTrack = 1:length(idxFast)
%     plot(xCoordsIn(idxFast(iTrack),:),yCoordsIn(idxFast(iTrack),:),'m'); % replot very fast in magenta
% end 


for iTrack = 1:length(idxPer) 
    plot(xCoordsIn(idxPer(iTrack),:),yCoordsIn(idxPer(iTrack),:),'b'); 
end 

for iTrack = 1:length(xCoordsOut(:,1))
    plot(xCoordsOut(iTrack,:),yCoordsOut(iTrack,:),'y'); 
end 


% plot connections
subTrackLength = length(xCoordsIn(:,1)); 
%  get idx of the first point of each track in the subregion  
firstPtIdx = arrayfun(@(i) find(~isnan(xCoordsIn(i,:)),1,'first'),1:subTrackLength); 

% extract coords where enters the mask  
xFirst = arrayfun(@(i) xCoordsIn(i,firstPtIdx(i)),1:subTrackLength);
yFirst = arrayfun(@(i) yCoordsIn(i,firstPtIdx(i)),1:subTrackLength); 


lastPtIdx = arrayfun(@(i) find(~isnan(xCoordsOut(i,:)),1,'last'),1:subTrackLength); 
xLast = arrayfun(@(i) xCoordsOut(i,lastPtIdx(i)),1:subTrackLength);
yLast= arrayfun(@(i) yCoordsOut(i,lastPtIdx(i)),1:subTrackLength); 

xConnect = [xFirst' xLast']; 
yConnect = [yFirst' yLast']; 

for iTrack = 1:subTrackLength
    plot(xConnect(iTrack,:),yConnect(iTrack,:),'y');
end 




% plot roi boundaries 
roiYX = bwboundaries(roiSet(:,:,iSub));
plot(roiYX{1}(:,2),roiYX{1}(:,1),'g'); 

 for iPole= 1:2 
     spy(poles(:,:,iPole),'c',10) 
 end
 
%%% REMOVE JITTER FROM HIGH DISPLACEMENT CLASSIFICATION %%% 
xMatLat = xCoordsIn(idxLat,:); 
yMatLat = yCoordsIn(idxLat,:); 


firstPtIdxLat =  arrayfun(@(i) find(~isnan(xMatLat(i,:)),1,'first'),1:length(idxLat)); 
lastPtIdxLat = arrayfun(@(i) find(~isnan(xMatLat(i,:)),1,'last'),1:length(idxLat));
xFirstLat = arrayfun(@(i) xMatLat(i,firstPtIdxLat(i)),1:length(idxLat));
xLastLat = arrayfun(@(i) xMatLat(i,lastPtIdxLat(i)),1:length(idxLat)); 

yFirstLat = arrayfun(@(i) yMatLat(i,firstPtIdxLat(i)),1:length(idxLat)); 
yLastLat = arrayfun(@(i) yMatLat(i,lastPtIdxLat(i)),1:length(idxLat)); 
dJitterCheck = sqrt((xFirstLat-xLastLat).^2 + (yFirstLat-yLastLat).^2);
dJitterCheck = dJitterCheck.*projData.pixSizeNm./1000; % convert to um
idxJitter = find(dJitterCheck<jitterParam); 
idxNonJitter =find(dJitterCheck>jitterParam); 

% plot the microtubule trajectories by type  
for iTrack = 1:length(idxNonJitter)
    
    plot(xMatLat(idxNonJitter(iTrack),:),yMatLat(idxNonJitter(iTrack),:),'r');
    hold on
end

velsHi = vel(idxLat); 

param4  = ['freqCorticalDispGreater' num2str(distCutoff) 'um_JitterOnly']; 
param5 = ['percentDispGreater' num2str(distCutoff) 'um_JitterOnly']; 

param6 = ['freqCorticalDispGreater' num2str(distCutoff) 'um_NoJitter']; 
param7 = ['percentDispGreater' num2str(distCutoff) 'um_NoJitter']; 

if ~isempty(idxJitter)
    
% jitter only 
numJitter = length(velsHi(idxJitter));

projData.stats.(param4) = numJitter/projData.nFrames/projData.secPerFrame*60; 
projData.stats.(param5) = numJitter/length(vel)*100; 
% plot jitter 
for iTrack = 1:length(idxJitter)
     plot(xMatLat(idxJitter(iTrack),:),yMatLat(idxJitter(iTrack),:),'m'); % replot jitter in magenta 
end 

 
 


% jitter filtered 
trueVelsHi = velsHi(idxNonJitter); 

projData.stats.(param6) = length(trueVelsHi)/projData.nFrames/projData.secPerFrame*60;  
projData.stats.(param7) = length(trueVelsHi)/length(vel)*100; 
projData.stats.hiDispVelMean = nanmean(trueVelsHi);
projData.stats.hiDispVelMeanNorm = nanmean(trueVelsHi)/astralGrowth;
else 
    projData.stats.(param6) = projData.stats.(param2); % no jitter so keep same as calculated above (no need to modify) 
    projData.stats.(param7) = projData.stats.(param3); 
    projData.stats.hiDispVelMean = nanmean(velsHi); 
    projData.stats.hiDispVelMeanNorm = projData.stats.hiDispVelMean/astralGrowth; 
end 

%%% SAVE SANITY PLOTS %%%

collectedDir2 = [up3 filesep 'CollectedSubRois' filesep 'Cortical' filesep 'LowVsHigh' filesep con]; 
  if ~isdir(collectedDir2)
    mkdir(collectedDir2)
 end 

% save this plot in folder and collected SubRoiFolder
saveas(gcf,[subRoiEdgeDir filesep 'lowVsHighDisp' con '_' cellNum '_pole_' num2str(iSub) '.eps'],'psc2');
saveas(gcf,[collectedDir2 filesep 'lowVsHighDisp' con '_' cellNum '_pole_' num2str(iSub)  '.tif']); 

 
 
%% ANGULAR MEASUREMENTS OF SUBTRACK INITIATION: FROM POLE TO CORTEX/POLE TO SUBTRACK INITIATION 


% calculation distance from pole to cortex 
dToPole = sqrt((point(iSub,1)-xFirst).^2 + (point(iSub,2)-yFirst).^2);
dToPole = dToPole*projData.pixSizeNm/1000; % convert to um 
dToIntersect =  sqrt((xIntersect(iSub)-xFirst).^2 + (yIntersect(iSub)-yFirst).^2);
dToIntersect = dToIntersect*projData.pixSizeNm/1000; 
dToCortex = p(iSub);
angle = acosd((dToCortex^2 + (dToPole).^2 - (dToIntersect).^2)./(2*dToCortex.*dToPole));   

% for now plot
% for i = 1:length(angle)
%     text(xFirst(i),yFirst(i),num2str(angle(i),3),'color','y');
% end 

angleLow = angle(idxPer); 
angleHi = angle(idxLat); 
angleHiTrue = angleHi(idxNonJitter);

projData.stats.meanAngleOfLowDispEvents = nanmean(angleLow); % in degrees
projData.stats.meanAngleOfHiDispEvents = nanmean(angleHiTrue); 


 
  
  


%%
%Make plots color-coded by value 
 close(gcf)
 
 cMapLength=128; cMap=jet(cMapLength);
 colorMapFig = figure('Visible','off');
 imshow(toPlot(:,:,2),[]); 
 hold on 
 plot(roiYX{1}(:,2),roiYX{1}(:,1),'w'); 
 if ~isempty(idxPer)
 % first plot dwells 
 lifeLow = life(idxPer); 
 data = lifeLow; 
 dataRange = 10; 
 data(data>dataRange)=dataRange;
                        mapper=linspace(0,dataRange,cMapLength)';
                        
                        % get closest colormap index for each feature
                        D=createDistanceMatrix(data,mapper);
                        [sD,idx]=sort(abs(D),2);
                     
                        xPlotDwell=xCoordsIn(idxPer,1:end-1)';
                        yPlotDwell=yCoordsIn(idxPer,1:end-1)';
                       
                        for k=1:cMapLength
                            plot(xPlotDwell(:,idx(:,1)==k),yPlotDwell(:,idx(:,1)==k),'color',cMap(k,:),'lineWidth',2);
                        end 
 
 for iTrack = 1:length(idxPer) 
    plot(xCoordsOut(idxPer(iTrack),:),yCoordsOut(idxPer(iTrack),:),'w'); 
 end 
 
% plot connections
 xCoordsInPer = xCoordsIn(idxPer,:); 
 yCoordsInPer = yCoordsIn(idxPer,:); 
firstPtIdxPer =  arrayfun(@(i) find(~isnan(xCoordsInPer(i,:)),1,'first'),1:length(idxPer)); 
xCoordsFirstIn = arrayfun(@(i) xCoordsInPer(i,firstPtIdxPer(i)),1:length(idxPer));
yCoordsFirstIn= arrayfun(@(i) yCoordsInPer(i,firstPtIdxPer(i)),1:length(idxPer)); 

xCoordsOutPer = xCoordsOut(idxPer,:); 
yCoordsOutPer = yCoordsOut(idxPer,:); 
lastPtIdxPer = arrayfun(@(i) find(~isnan(xCoordsOutPer(i,:)),1,'last'),1:length(idxPer)); 
xCoordsLastOut = arrayfun(@(i) xCoordsOutPer(i,lastPtIdxPer(i)),1:length(idxPer));
yCoordsLastOut= arrayfun(@(i) yCoordsOutPer(i,lastPtIdxPer(i)),1:length(idxPer)); 

connectX = [xCoordsFirstIn' xCoordsLastOut']; 
connectY = [yCoordsFirstIn' yCoordsLastOut']; 

for iTrack =1 :length(connectX(:,1))
    plot(connectX(iTrack,:),connectY(iTrack,:),'w'); 
end 
 end %  isempty idxPer   

collectedDir3 = [up3 filesep 'CollectedSubRois' filesep 'Cortical' filesep 'DwellTimeMaps' filesep con];
if ~isdir(collectedDir3) 
    mkdir(collectedDir3) 
end 

saveas(gcf,[subRoiEdgeDir filesep 'DwellMap' con '_' cellNum '_pole_' num2str(iSub) '.eps'],'psc2');
saveas(gcf,[collectedDir3 filesep 'DwellMap' con '_' cellNum '_pole_' num2str(iSub)  '.tif']); 

      delete(findobj(colorMapFig,'Type','line'));
       plot(roiYX{1}(:,2),roiYX{1}(:,1),'w'); 
if ~isempty(idxLat)
% plot lateral velocities 
velFast = vel(idxLat); 

data = velFast; 
dataRange = 20; 
data(data>dataRange)=dataRange; 
 mapper=linspace(0,dataRange,cMapLength)';
                        
                        % get closest colormap index for each feature
                        D=createDistanceMatrix(data,mapper);
                        [sD,idx]=sort(abs(D),2);
                     
                        xPlotVel=xCoordsIn(idxLat,1:end-1)';
                        yPlotVel=yCoordsIn(idxLat,1:end-1)';
                       
                        for k=1:cMapLength
                            plot(xPlotVel(:,idx(:,1)==k),yPlotVel(:,idx(:,1)==k),'color',cMap(k,:),'lineWidth',2);
                        end 
                        
  
 for iTrack = 1:length(idxLat) 
    plot(xCoordsOut(idxLat(iTrack),:),yCoordsOut(idxLat(iTrack),:),'w'); 
    
 end                       
    
 % plot connections
 xCoordsInLat = xCoordsIn(idxLat,:); 
 yCoordsInLat = yCoordsIn(idxLat,:); 
firstPtIdxLat =  arrayfun(@(i) find(~isnan(xCoordsInLat(i,:)),1,'first'),1:length(idxLat)); 
xCoordsFirstIn2 = arrayfun(@(i) xCoordsInLat(i,firstPtIdxLat(i)),1:length(idxLat));
yCoordsFirstIn2= arrayfun(@(i) yCoordsInLat(i,firstPtIdxLat(i)),1:length(idxLat)); 

xCoordsOutLat = xCoordsOut(idxLat,:); 
yCoordsOutLat = yCoordsOut(idxLat,:); 
lastPtIdxLat = arrayfun(@(i) find(~isnan(xCoordsOutLat(i,:)),1,'last'),1:length(idxLat)); 
xCoordsLastOut2 = arrayfun(@(i) xCoordsOutLat(i,lastPtIdxLat(i)),1:length(idxLat));
yCoordsLastOut2= arrayfun(@(i) yCoordsOutLat(i,lastPtIdxLat(i)),1:length(idxLat)); 

 
 connectX2 = [xCoordsFirstIn2' xCoordsLastOut2']; 
connectY2 = [yCoordsFirstIn2' yCoordsLastOut2']; 

for iTrack =1 :length(connectX2(:,1))
    plot(connectX2(iTrack,:),connectY2(iTrack,:),'w'); 
end 
end 
 collectedDir4 = [up3 filesep 'CollectedSubRois' filesep 'Cortical' filesep 'SlideVelMaps' filesep con];
 
 if ~isdir(collectedDir4) 
     mkdir(collectedDir4) ;
 end 
 
 saveas(gcf,[subRoiEdgeDir filesep 'LatVelocity' con '_' cellNum '_pole_' num2str(iSub) '.eps'],'psc2');
saveas(gcf,[collectedDir4 filesep 'LatVelocity' con '_' cellNum '_pole_' num2str(iSub)  '.tif']);      
close(gcf)
else % no tracks make criteria 
    vel = nan;
    dist = nan; 
    life = nan; 
    
    
    projData.stats.staticDwell = nan;
    projData.stats.lowDispVelMean = nan; 
    projData.stats.lowDispVelMeanNorm = nan; 

    projData.stats.freqCortical = 0 ; % no tracks in region in sample time 
    projData.stats.(param2) = 0;
    projData.stats.(param3) = nan; % can't calculate probabilility of end-on to side-on transition if no tracks in region
    projData.stats.(param4) = 0; % freqCorticalDispGreaterXum_JitterOnly
    projData.stats.(param5) = nan; % percentDispGreaterXum_JitterOnly 
    projData.stats.(param6) = 0; % freqCorticalDispGreaterXum_NoJitter
    projData.stats.(param7) = nan; % percentDispGreaterXum_NoJitter
    projData.stats.hiDispVelMean = nan; 
    projData.stats.hiDispVelMeanNorm = nan; 
    projData.stats.meanAngleOfLowDispEvents =nan; % in degrees
    projData.stats.meanAngleOfHiDispEvents = nan; 
    angle = NaN; 
end 


% sanity check 
% imshow(toPlot(:,:,2),[])
% hold on 


projData.corticalDataMat = [vel life dist angle']; 
projData.poleCortexDist = p(iSub); 
% update projData
save([subRoiEdgeDir filesep 'meta' filesep 'projData.mat'],'projData'); 
% e
projData.saveDir = [subRoiEdgeDir filesep 'trackMovie']; 
if ~isdir(projData.saveDir) 
    mkdir(projData.saveDir)
end 

roiYX = bwboundaries(roiSet(:,:,iSub)); 
roiYX = roiYX{1}; 
if makeMovie ==1 
plusTipTrackMovie_new(projData,[],[],roiYX,[],1,3,0,0,[],[]); 
end 
end
clear roiSet forIntersectAll r11 r12 forIntersect toPlot
end




