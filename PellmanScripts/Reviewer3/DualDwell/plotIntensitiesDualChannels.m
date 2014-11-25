function [ output_args ] = plotIntensitiesDualChannels( infoDir, saveDir,makeMovie)
% INPUT: 
% infoDir : the character array with the path name to the directory above
% the two channels 
% saveDir : the character array with the path name to the output directory 
% 
% OUTPUT: makes an output directory 'DwellPlots' with the intensities of 
% each of the two  channels in time right before and after enter small
% dwell window. 

figDir = [saveDir filesep 'DualChannelDwellPlots']; 
if ~isdir(figDir)
    mkdir(figDir)
end 

movDir = [saveDir filesep 'DualChannelMoviePlots']; 
if ~isdir(movDir) 
    mkdir(movDir);
end 



imDir1 = [infoDir filesep 'w2491' filesep 'images' ]; 
imDir2 = [infoDir filesep 'w3561']; 

subRoiDir = [infoDir filesep 'w2491' filesep 'roi_1' filesep 'subRoisCortical' filesep 'sub_1']; 

load([subRoiDir filesep 'dwellMasks']);
load([subRoiDir filesep 'meta' filesep 'projData.mat']);
dataMat = projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix; % load to find which tracks end in term,pause, or shrinkage 
dataMat = dataMat(dataMat(:,9) ~= 0,:); % filter the dataMat to get rid of pause, shrinkage, and undefined gaps. 
% idxEndTerm = dataMat(:,9) == 1; % collect the terminal growth events 
% idxEndPause = dataMat(:,9) == 2; 
% idxEndShrink = dataMat(:,9) == 3; % collect the growth events ending in shrinkage 
% % collect the growth events ending in a pause ; 
% idxEndUndefined = dataMat(:,9) == 4; 

% make all the directories for all the different types of dwells. 
type{1} = 'Terminal'; 
type{2} = 'Pause'; 
type{3} = 'Shrinkage'; 
type{4} = 'Undefined';
% make all the different movie directories and figure directories for 
% dwell events ending in term, shrinkage,pause, and undefined events. 
movieDirs = cellfun(@(x) [movDir filesep x],type,'uniformoutput',0); 
figDirs = cellfun(@(x) [figDir filesep x],type,'uniformoutput',0); 
for i = 1:numel(movieDirs) 
    if ~isdir(movieDirs{i})
        mkdir(movieDirs{i}); 
    end 
end 

for i = 1:numel(figDirs) 
    if ~isdir(figDirs{i})
        mkdir(figDirs{i}) 
    end 
end 

    

xMat = projData.xCoordDwell;
yMat = projData.yCoordDwell;
xMatOut = projData.xCoordDwellOutAllTracks;
yMatOut = projData.yCoordDwellOutAllTracks; 

listOfImages = searchFiles('.tif',[],imDir1,0); 
listOfImages2 = searchFiles('.tif',[],imDir2,0); 
[sortedList,sortNum] = sortUnpaddedList(listOfImages);
[sortedList2,sortNum2] = sortUnpaddedList(listOfImages2); 
for iFrame = 1:length(listOfImages)
 % just load all now so don't have to keep reloading them 
 img1(:,:,iFrame) =  double(imread([char(sortedList(iFrame,1)) filesep char(sortedList(iFrame,2)),...;
        num2str(sortNum(iFrame)) char(sortedList(iFrame,4))]));
 img2(:,:,iFrame) = double(imread( [char(sortedList2(iFrame,1)) filesep char(sortedList2(iFrame,2)),...;
        num2str(sortNum2(iFrame)) char(sortedList2(iFrame,4))]));
end 
trunc = 1;



% check to make sure that the length of the dataMat after remove all gaps
% is the same as the length of the dwell time. 
if numel(framesDwell) ~= size(dataMat,1)
    warning(['Indexing incorrect for' infoDir]); 
end 
    

if trunc == 1
    
% truncate by long dwell times (more than 5 frames) 
time = cellfun(@(x) length(x),framesDwell); % framesDwell should be saved in the dwell masks .mat file 
idxToPlot = find(time>=4); 
%dataMat = dataMat(time>=4,:); % dataMat should now correspond to only those tracks that are in the dwell window for greater than 
% 4 frames 
%endType = dataMat(:,9); % get the identifiers for the termination type. 
else 
    
    idxToPlot = 1:size(dwellMasks,3); 
end 

for i= 1:length(idxToPlot)
    endType = dataMat(idxToPlot(i),9); % get if the track was a 'terminal = 1', a shrinkage = 3, a pause = 2, or an 
    % undefined event (4). 
    iDwell = idxToPlot(i); 
    frameC = framesDwell{iDwell}; 
    
   
    framesToPlot = [frameC(1)-5:frameC(1)-1, frameC, frameC(end)+1: frameC(end)+5];
   framesToPlot =  framesToPlot(framesToPlot>0 & framesToPlot<101);
    %if makeMovie ==1 
    % get the index  
     
        % save in the directory corresponding to the correct file 
        dwellDirC = [movieDirs{endType} filesep 'Dwell' num2str(iDwell,'%02d')]; 
        
        if ~isdir(dwellDirC);
                mkdir(dwellDirC); 
        end 
    %end 
    
                
    for iFrame = 1:length(framesToPlot)
    
        % get all the values from the image 
    img1C = img1(:,:,framesToPlot(iFrame));
    b = min(img1C(:));
    c = max(img1C(:)); 
    img1CN =  (img1C-b)./(c-b); 
    clear b c 
   
    img2C = img2(:,:,framesToPlot(iFrame)); 
    b = min(img2C(:)); 
    c = max(img2C(:)); 
    img2CN = (img2C-b)./(c-b); 
     
    values1(iFrame) =  mean(img1C(logical(dwellMasks(:,:,iDwell)))); 
    values2(iFrame) = mean(img2C(logical(dwellMasks(:,:,iDwell)))); 
    %% MARIA HERE IS WHERE YOU WOULD ADD IF THERE IS A CORRELATION ETC 
   % find the maximum value within the dwell window
    
    
    
    
    error1(iFrame) = std(img1C(logical(dwellMasks(:,:,iDwell))));
    error2(iFrame) = std(img2C(logical(dwellMasks(:,:,iDwell)))); 
  
        if makeMovie == 1
            dwellMask = dwellMasks(:,:,iDwell); 
        [maskCrop,x,y,imgToPlot(:,:,2)] = cropImageBasedOnMask(dwellMask,5, img1CN); % plotEB in green
        [~,~,~,imgToPlot(:,:,1)] = cropImageBasedOnMask(dwellMask,5,img2CN); % 
        imgToPlot(:,:,3) = zeros(size(imgToPlot(:,:,1))); 
        [ny,nx,~] = size(img1); 
        setFigure(nx,ny,'on'); 
        image(imgToPlot)
        hold on 
        roiYX =  bwboundaries(maskCrop); 
        
        cellfun(@(x) plot(x(:,2),x(:,1),'w'),roiYX);
        % cropping dwell differently for each mask to have to update here
        xMatRect = xMat - min(x) +1; 
         yMatRect = yMat -min(y) +1; 
        %
         xCoord = xMatRect(iDwell,framesToPlot(iFrame)); 
         yCoord = yMatRect(iDwell,framesToPlot(iFrame)); 
         scatter(xCoord,yCoord,100,'filled','w'); 
         
         xMatOutRect = xMatOut-min(x)+1; 
         yMatOutRect = yMatOut-min(y)+1; 
         xCoordOut = xMatOutRect(iDwell,framesToPlot(iFrame)); 
         yCoordOut = yMatOutRect(iDwell,framesToPlot(iFrame)); 
         scatter(xCoordOut,yCoordOut,100,'filled','b');
         text(2,2,num2str(framesToPlot(iFrame),'%03d'),'color','w','FontSize',16); 
    saveas(gcf, [dwellDirC filesep 'DwellMovie' num2str(framesToPlot(iFrame),'%03d') '.png']); 
    
    close gcf
    
        end 
    end 
    
    figure('Visible','on')
 %% MODIFIED 20141107
 dwellInt.channel1.values{endType}{i} = values1; % currently index by the number 
 dwellInt.channel2.values{endType}{i} = values2; 
 
 % find the max of the series 
idxMaxSeries = find(values1 == max(values1));

 
 
 % find the max in the dwell window
 valuesDwell = values1(6:5+length(frameC)); 
 idxMaxDwell = find(values1 == max(valuesDwell)); 
  
 
 % ADDED 20141124 
 %Find range 
 range1 = max(values1) - min(values1); 
 range2 = max(values2)- min(values2); 
 dwellInt.channel1.range1{endType}{i} = range1; 
 dwellInt.channel2.range2{endType}{i} = range2; 
 
 
 % test if the max in the dwell window is the max of the entire series 
 % if it isn't discard 
    if idxMaxSeries == idxMaxDwell 
        dwellInt.channel1.maxInWind{endType}{i} =  1;
    else 
        dwellInt.channel1.maxInWind{endType}{i} = 0; 
    end  
    
         
 %%
 
 
 
subplot(2,1,1); 

 errorbar(framesToPlot,values1,error1,'color','g'); 
 hold on 
  scatter(framesToPlot, values1,100,'g','filled'); 
 

 line([frameC(1) frameC(1)],[min(values1)-100 max(values1)+100],'color','k','Linewidth',2); 
 line([frameC(end) frameC(end)],[min(values1)-100 max(values1)+100],'color','k','Linewidth',2); 
 ylabel('Raw Intensity Values (AU) in Dwell Window'); 
 xlabel('Frames (0.85 sec intervals)')
 % plot max value within well in the first channel 
 line([framesToPlot(idxMaxDwell),framesToPlot(idxMaxDwell)],[min(values1)-100 max(values1)+100],'color','r','Linewidth',2); 
 % plot values decay 
  valuesDecay1 = values1(idxMaxDwell:5+length(frameC)); 
  framesDecay = framesToPlot(idxMaxDwell:5+length(frameC)); 
  scatter(framesDecay,valuesDecay1,100','k'); % outline the values used for decay in black 
  slopeValues1 = diff(valuesDecay1);
  meanSlopeValues1 = mean(slopeValues1); 
 dwellInt.channel1.slopeValues{endType}{i} = slopeValues1; 
 dwellInt.channel1.meanSlope{endType}{i} = meanSlopeValues1; 
 
  
  text(framesToPlot(idxMaxDwell),values1(idxMaxDwell),['m = ' num2str(meanSlopeValues1,3)]); 
  
    
 %   
subplot(2,1,2); 
errorbar(framesToPlot,values2,error2,'color','r'); 
hold on 

 scatter(framesToPlot,values2,100,'r','filled'); 
 line([frameC(1) frameC(1)],[min(values2)-100 max(values2)+100],'color','k','Linewidth',2); 
 line([frameC(end) frameC(end)],[min(values2)-100 max(values2)+100],'color','k','Linewidth',2);  
 %line([frameC(idx2) frameC(idx2)]
 ylabel('Raw Intensity Values(AU) in Dwell Window'); 
 xlabel('Frames (0.85 sec intervals'); 

 valuesDecay2 = values2(idxMaxDwell:5+length(frameC)); 
  %framesDecay = framesToPlot(idxMaxDwell:5+length(frameC)); 
  scatter(framesDecay,valuesDecay2,100','k'); % outline the values used for decay in black 
  slopeValues2 = diff(valuesDecay2);
  meanSlopeValues2 = mean(slopeValues2); 
 dwellInt.channel2.slopeValues{endType}{i} = slopeValues2; 
 dwellInt.channel2.meanSlope{endType}{i} = meanSlopeValues2; 
 
  line([framesToPlot(idxMaxDwell),framesToPlot(idxMaxDwell)],[min(values2)-100 max(values2)+100],'color','r','Linewidth',2); 
  text(framesToPlot(idxMaxDwell),values2(idxMaxDwell),['m = ' num2str(meanSlopeValues2,3)]); 
 
 test = (meanSlopeValues1<0 & meanSlopeValues2<0); 
 dwellInt.consistencyTest{endType}{i} = test; 
 if test == 1 
 title('Potential Consistent Decay'); 
 end
 
 
 
 saveas(gcf,[figDirs{endType} filesep 'Dwell_' num2str(iDwell,'%03d') '.fig'])
 saveas(gcf,[figDirs{endType} filesep 'Dwell_' num2str(iDwell,'%03d') '.tif'])

 clear values1 values2 framesToPlot error1 error2 
 close gcf
 
  
  
  
  
     
  
end

save([saveDir filesep 'dwellInt.mat'],'dwellInt');

end 

