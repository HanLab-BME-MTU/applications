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
%subRoiDir = [infoDir filesep 'w2491' filesep 'roi_1' filesep 'subRois_FullEdge_20150131_WithFrameInfo\sub_1']; 

 subRoiDir = [infoDir filesep 'w2491' filesep 'roi_1' filesep 'subRoisCortical' filesep 'sub_1'];
maskDir = [subRoiDir filesep 'masks'];
load([subRoiDir filesep 'dwellMasks']);
load([subRoiDir filesep 'meta' filesep 'projData.mat']);
dataMat = projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix; % load to find which tracks end in term,pause, or shrinkage
dataMat = dataMat(dataMat(:,9) ~= 0,:); % filter the dataMat to get rid of pause, shrinkage, and undefined gaps.
xMat = projData.xCoordDwell;
yMat = projData.yCoordDwell;
xMatOut = projData.xCoordDwellOutAllTracks;
yMatOut = projData.yCoordDwellOutAllTracks;

%%
filterCorticalClass = 0;
%FILTER BY CORTICAL : ADDED 20150201
if filterCorticalClass == 1;
    
    dwells = projData.dwellAllTracks;
    disp = vertcat(projData.dispByFrame{:});
    
    load([subRoiDir filesep 'meta' filesep 'CorticalInfo' filesep 'corticalData.mat' ]);
    idxLatLog = corticalData.params3_60_0pt7_0pt3.idxLatLogical;
    % get the indices for each class (indices are for dataMat filtered)
    %idxEndTerm = dataMat(:,9) == 1; % collect the terminal growth events
    %idxEndPause = dataMat(:,9) == 2;
    %idxEndShrink = dataMat(:,9) == 3; % collect
    
    %percentPause = sum(idxEndPause)./length(idxEndPause)*100;
    
    dispDiscar = 0.3;
    % filter by displacement
    %dwells = dwells(disp>dispDiscard & ~isnan(disp));
    % ok I know this is a bit weird here but filter by these when I do the
    % classifications
    %idxEndTerm = idxEndTerm(disp>dispDiscard & ~isnan(disp));
    %idxEndPause = idxEndPause(disp>dispDiscard & ~isnan(disp));
    %idxEndShrink = idxEndShrink(disp>dispDiscard & ~isnan(disp));
    
    % filter dataMat by displacement filter
    dataMat = dataMat(disp>dispDiscar & ~isnan(disp),:);
    % filter dwellmasks by displacement filter
    dwellMasks = dwellMasks(:,:,disp>dispDiscar & ~isnan(disp));
    xMat = xMat(disp>dispDiscar & ~isnan(disp),:);
    yMat = yMat(disp>dispDiscar & ~isnan(disp),:);
    
    xMatOut = xMatOut(disp>dispDiscar & ~isnan(disp),:);
    yMatOut = yMatOut(disp>dispDiscar& ~isnan(disp),:);
    dwells = dwells(disp>dispDiscar& ~isnan(disp));
    framesDwell = framesDwell(disp>dispDiscar&~isnan(disp));
    % get only end on dwells
    %dwellsEndOnP = dwells(~idxLatLog & idxEndPause);
    %dwellsEndOnST = dwells(~idxLatLog & (idxEndTerm | idxEndShrink));
    
    %filter out lateral transition tracks for both the dataMat and the dwellMasks
    dataMat = dataMat(~idxLatLog,:);
    dwellMasks = dwellMasks(:,:,~idxLatLog);
    xMat = xMat(~idxLatLog,:);
    yMat = yMat(~idxLatLog,:);
    xMatOut = xMatOut(~idxLatLog,:);
    yMatOut = yMatOut(~idxLatLog,:);
    dwells = dwells(~idxLatLog);
    framesDwell = framesDwell(~idxLatLog);
    
end


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





listOfImages = searchFiles('.tif',[],imDir1,0);
listOfImages2 = searchFiles('.tif',[],imDir2,0);
[sortedList,sortNum] = sortUnpaddedList(listOfImages);
[sortedList2,sortNum2] = sortUnpaddedList(listOfImages2);

if exist(maskDir)~=0
    listOfMasks = searchFiles('.tif',[],maskDir,0,'all',1);
    multMask = 1;
else
    
    roiMask = logical(imread([subRoiDir filesep 'roiMask.tif']));
    multMask = 0;
end

for iFrame = 1:size(listOfImages,1)
    % just load all now so don't have to keep reloading them
    img1(:,:,iFrame) =  double(imread([char(sortedList(iFrame,1)) filesep char(sortedList(iFrame,2)),...;
        num2str(sortNum(iFrame)) char(sortedList(iFrame,4))]));
    img2(:,:,iFrame) = double(imread( [char(sortedList2(iFrame,1)) filesep char(sortedList2(iFrame,2)),...;
        num2str(sortNum2(iFrame)) char(sortedList2(iFrame,4))]));
    if multMask == 1;
        roiMask = logical(imread(char(listOfMasks(iFrame,1))));
    end
    maskImg = roiMask.*img2(:,:,iFrame);
    values = maskImg(maskImg~=0);
    [x,noise]  = fitGaussianModeToPDF(values);
%     setAxis
%     hist(values,100,'color','k'); 
%    % bottom= min(values);
%     %top = max(values); 
%     line([x,x],[0 300],'color','k','Linewidth',2); 
%     line([x+noise, x+noise],[0 1000],'color','k','Linewidth',2); 
%     line([x+2*noise,x+2*noise],[0 1000],'color','r','Linewidth',2); 
    meansChannel2(iFrame) =x;
    stdChannel2(iFrame) = noise;
end

% estimate std for each frame


check =0;
[ny,nx] = size(img1(:,:,1));
setFigure(nx,ny,'on');

if check == 1
    hold on
    
    imshow(-img1(:,:,1),[]);
    hold on
    
    arrayfun(@(i) plot(xMat(i,:),yMat(i,:),'color','r'),1:size(xMat,1));
    arrayfun(@(i) plot(xMatOut(i,:),yMatOut(i,:),'color','k'),1:size(xMatOut,1));
    timeFlags = ~isnan(xMat);
    starts =  arrayfun(@(i) find(timeFlags(i,:),1,'first'),1:size(xMat,1));
    
    arrayfun(@(i) text(xMat(i,starts(i)),yMat(i,starts(i)),num2str(i)),1:size(xMat,1)); % the ID
    % arrayfun(@(i) text(xMat(i,starts(i)),yMat(i,starts(i)),num2str(dwells(i,1),3)),1:size(xMat,1));
    saveas(gcf,[saveDir filesep 'testID.fig']);
end
trunc = 0;



% check to make sure that the length of the dataMat after remove all gaps
% is the same as the length of the dwell time.
if numel(framesDwell) ~= size(dataMat,1)
    warning(['Indexing incorrect for' infoDir]);
end


% if trunc == 1
%     
%     % truncate by long dwell times (more than 5 frames)
%     time = cellfun(@(x) length(x),framesDwell); % framesDwell should be saved in the dwell masks .mat file
%     idxToPlot = find(time>=2);
%     %dataMat = dataMat(time>=4,:); % dataMat should now correspond to only those tracks that are in the dwell window for greater than
%     % 4 frames
%     %endType = dataMat(:,9); % get the identifiers for the termination type.
% else
%     
%     idxToPlot = 1:size(dwellMasks,3);
% end
idxToPlot = 25;
for i= 1:length(idxToPlot)
    endType = dataMat(idxToPlot(i),9); % get if the track was a 'terminal = 1', a shrinkage = 3, a pause = 2, or an
    % undefined event (4).
    iDwell = idxToPlot(i);
    frameC = framesDwell{iDwell};
    
    
    framesToPlot = [frameC(1)-5:frameC(1)-1, frameC, frameC(end)+1: frameC(end)+5];
    framesToPlot =  framesToPlot(framesToPlot>0 & framesToPlot<101);
    
    % get the estimated variance of the image intensities for the current
    % frames
    stdC =  mean(stdChannel2(framesToPlot));
    %if makeMovie ==1
    % get the index
    
    % save in the directory corresponding to the correct file
    dwellDirC = [movieDirs{endType} filesep 'Dwell' num2str(iDwell,'%02d')];
    
    if ~isdir(dwellDirC);
        mkdir([dwellDirC filesep 'Dual']);
        mkdir([dwellDirC filesep 'EBChannel']);
        mkdir([dwellDirC filesep 'mCherryChannel']);
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
        if iDwell ==25;
            if makeMovie == 1
                
                dwellMask = dwellMasks(:,:,iDwell);
                [maskCrop,x,y,imgToPlot(:,:,2)] = cropImageBasedOnMask(dwellMask,20, img1CN); % plotEB in green
                [~,~,~,imgToPlot(:,:,1)] = cropImageBasedOnMask(dwellMask,20,img2CN); %
                imgToPlot(:,:,3) = zeros(size(imgToPlot(:,:,1)));
                [ny,nx,~] = size(img1);
                valuesAll = imgToPlot(:);
                valuesAll = valuesAll(valuesAll~=0);
                minVal = min(valuesAll);
                maxVal = max(valuesAll);
                setFigure(nx,ny,'on');
                imshow(imgToPlot,[minVal,maxVal])
                set(gcf,'Color',[1 1 1]); 
                hold 
                roiYX =  bwboundaries(maskCrop);
                
                cellfun(@(x) plot(x(:,2),x(:,1),'w','LineWidth',1),roiYX);
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
                text(5,5,num2str(framesToPlot(iFrame)),'color','w','FontSize',30);
                saveas(gcf, [dwellDirC filesep 'Dual' filesep 'DwellMovie' num2str(framesToPlot(iFrame),'%03d') '.png']);
                saveas(gcf, [dwellDirC filesep 'Dual' filesep 'DwellMovie' num2str(framesToPlot(iFrame),'%03d') '.fig']); 
                saveas(gcf,[dwellDirC filesep 'Dual' filesep 'DwellMovie' num2str(framesToPlot(iFrame),'%03d') '.eps'],'psc2'); 
                
                close gcf
                
                
                % plot other channel in folder
                % get limits
                chName{1} = 'mCherryChannel';
                chName{2} = 'EBChannel';
                for iCh = 1:2
                    setFigure(nx,ny,'on');
                    
                    
                    imshow(imgToPlot(:,:,iCh),[minVal,maxVal]);
                    %imshow(imgToPlot(:,:,1),[minVal,maxVal])
                    hold on
                    roiYX =  bwboundaries(maskCrop);
                    
                    cellfun(@(x) plot(x(:,2),x(:,1),'y','LineWidth',1),roiYX);
                    scatter(xCoordOut,yCoordOut,100,'filled','r');
                    scatter(xCoord,yCoord,100,'filled','y');
                    %text(5,5,num2str(framesToPlot(iFrame)),'color','w','FontSize',16);
                    set(gcf,'Color',[1,1,1]); 
                    saveas(gcf, [dwellDirC filesep chName{iCh} filesep 'DwellMovie' num2str(framesToPlot(iFrame),'%03d') '.png']);
                    saveas(gcf, [dwellDirC filesep chName{iCh} filesep 'DwellMovie' num2str(framesToPlot(iFrame),'%03d') '.fig']);
                     saveas(gcf, [dwellDirC filesep chName{iCh} filesep 'DwellMovie' num2str(framesToPlot(iFrame),'%03d') '.eps'],'psc2');
                    close gcf
                end % for iCh
            end
        end %if
    end
    
    figure('Visible','off')
    %% MODIFIED 20141107
    dwellInt(i).channel1.values = values1; % currently index by the number
    dwellInt(i).channel2.values = values2;
    dwellInt(i).channel2.stdNoise = stdC;
    % find the max of the series
    idxMaxSeries = find(values1 == max(values1));
    
    
    
    % find the max in the dwell window
    valuesDwell = values1(6:5+length(frameC));
    idxMaxDwell = find(values1 == max(valuesDwell));
    
    
    % ADDED 20141124
    %Find range
    range1 = max(values1) - min(values1);
    range2 = max(values2)- min(values2);
    enoughSig = range2>2*stdC;
    dynamicRange1 = max(values1)/min(values1);
    dynamicRange2 = max(values2)/min(values2);
    
    dwellInt(i).channel1.range1 = range1;
    dwellInt(i).channel2.range2 = range2;
    dwellInt(i).channel1.dynamicRange = dynamicRange1;
    dwellInt(i).channel2.dynamicRange = dynamicRange2;
    
    dwellInt(i).endType = endType;
    dwellInt(i).enoughSig = enoughSig;
    
    % test if the max in the dwell window is the max of the entire series
    % if it isn't discard
    if idxMaxSeries == idxMaxDwell
        dwellInt(i).channel1.maxInWind =  1;
    else
        dwellInt(i).channel1.maxInWind = 0;
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
    if length(idxMaxDwell)>1
        idxMaxDwell = idxMaxDwell(1);
    end
    % plot max value within well in the first channel
    line([framesToPlot(idxMaxDwell),framesToPlot(idxMaxDwell)],[min(values1)-100 max(values1)+100],'color','r','Linewidth',2);
    % plot values decay
    valuesDecay1 = values1(idxMaxDwell:5+length(frameC));
    framesDecay = framesToPlot(idxMaxDwell:5+length(frameC));
    scatter(framesDecay,valuesDecay1,100','k'); % outline the values used for decay in black
    slopeValues1 = diff(valuesDecay1);
    meanSlopeValues1 = mean(slopeValues1);
    dwellInt(i).channel1.slopeValues = slopeValues1;
    dwellInt(i).channel1.meanSlope = meanSlopeValues1;
    
    
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
    dwellInt(i).channel2.slopeValues = slopeValues2;
    dwellInt(i).channel2.meanSlope = meanSlopeValues2;
    
    line([framesToPlot(idxMaxDwell),framesToPlot(idxMaxDwell)],[min(values2)-100 max(values2)+100],'color','r','Linewidth',2);
    text(framesToPlot(idxMaxDwell),values2(idxMaxDwell),['m = ' num2str(meanSlopeValues2,3)]);
    
    test = (meanSlopeValues1<0 & meanSlopeValues2<0);
    dwellInt(i).consistencyTest = test;
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

