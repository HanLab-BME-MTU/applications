function [ output_args ] = wrapperThresholding( projList,threshMethod )
%Small Wrapper to perform Rosin Thresholding To Obtain An Approximate 
% Cell Edge. Perform Thresholding on the intensities of the entire image. 
%INPUT: projList a N Movie structure array with fields 
%      .imDir :  (string) path to the image directory
%      .anDir : (string) path to the roi directory
%  threshMethod: segmentation method, currently supports minAfterFirstMax
%  and Rosin Thresholding but others can be added. 
% OUTPUT: Creates a file 'masks' in the projList(iProj).anDir containing 
%        .tif files for the entire movie. Note no smoothing of the edge is
%        performed however the largest conneected component is chosen. 
%        Inside the masks folder a subFolder overlays is likewise created
%        where the mask outline is plotted over the original image for each
%        frame. 
if nargin<1 || isempty(projList) 
    [file ,path] = uigetfile(pwd,'Please Select a projList.mat file'); 
     load([path filesep file]); 
end 
saveSumMask = 1; % flag to create the summation mask, currently just default
% always on
if nargin<2|| isempty(threshMethod)
      % ask user which type to use 
        names = {'Rosin','MinMax','Otsu'}; 
stopFlag = 0; 
while stopFlag == 0
    
    [s,v] = listdlg('PromptString',{'Please Select'; 'a Threshold Method';'To Compare Methods' ; 'Select All'},...
        'ListString',names);
    
    if v == 1
      if length(s) == 1
        threshMethod = names{s} ;
        close gcf
       stopFlag = 1;
      else 
          % get first image 
          % do calc for both overlay 
          % 
            listOfImages = searchFiles('.tif',[],[projList(1).imDir],0,'all',1);
          
           
               
        img = double(imread(listOfImages{1}));
        imgFiltered = filterGauss2D(img,1);
        maxSignal = max(imgFiltered(:)); 
        minSignal = min(imgFiltered(:)); 
        %normalize image between 0 and 1
imgFilteredNorm = (imgFiltered - minSignal)/(maxSignal - minSignal);
      
      
        [levelRosin]  = thresholdRosin(imgFilteredNorm);
     
   [levelMaxMin] = thresholdFluorescenceImage(imgFilteredNorm);
   
   [levelOtsu] = graythresh(imgFilteredNorm);
        
        maskR = imgFilteredNorm>=levelRosin; % get mask
        maskM = imgFilteredNorm>=levelMaxMin; 
        maskO = imgFilteredNorm>=levelOtsu;
        % clean it up
        maskR = imfill(maskR,'holes');
        maskM = imfill(maskM,'holes'); 
        maskO= imfill(maskO,'holes'); 
        % get largest CC
        maskR = logical(getLargestCC(maskR));
        % mask = imopen(mask,strel('disk',3));
        maskM = logical(getLargestCC(maskM));
        % maskDir = [projList{iProj} filesep 'masks'];
        maskO = logical(getLargestCC(maskO)); 
        
        roiYXR = bwboundaries(maskR);
        roiYXM = bwboundaries(maskM); 
        roiYXO = bwboundaries(maskO); 
       figure;
        subplot(1,2,1)
        imshow(-img,[]);
        hold on
         cellfun(@(x) plot(x(:,2),x(:,1),'b'),roiYXR);
         cellfun(@(x) plot(x(:,2),x(:,1),'m'),roiYXM);
         cellfun(@(x) plot(x(:,2),x(:,1),'g'),roiYXO); 
         text(5,10,'Rosin','color','b'); 
         text(5,25,'MinMax','color','m'); 
         text(5,40,'Otsu','color','g'); 
          subplot(1,2,2)
          hist(imgFilteredNorm(:),100)
          xlabel('Normalized Img Intensity')
          ylabel('Number of Pixels'); 
         
          hold on 
          [x,~] = hist(imgFilteredNorm(:),100); 
       
          line([levelRosin,levelRosin],[0,max(x)],'color','b'); 
             line([levelMaxMin,levelMaxMin],[0 max(x)],'color','m'); 
             line([levelOtsu,levelOtsu],[0 max(x)],'color','g'); 
              ms =  msgbox('Examine Thresholding: Zoom if Needed, Click OK When Finished');
            uiwait(ms)
      end     
          
    else
        x = questdlg('Try Again?');
        if (strcmpi(x,'no') || strcmpi(x,'cancel'))
            stopFlag =1;
       error('User Exited: No Analysis Performed'); 
        
        end 
    end
end
end 

for iProj  =1 : numel(projList)
    display(['Performing ' threshMethod 'Thresholding for ' projList(iProj).imDir])
    % get images for each movie
    %upOne = upDirectory([projList{iProj,1}],1);
    upOne = upDirectory(projList(iProj).imDir,1);
    % listOfImages = searchFiles('.tif',[],[upOne filesep 'images'],0);
    listOfImages = searchFiles('.tif',[],[projList(iProj).imDir]);
    
%     [sortedImages, sortednum]  = sortUnpaddedList(listOfImages);
    sortedImages = listOfImages; 
    % make masks
    
    for iMask = 1:size(sortedImages,1);
        % load
%         fileNameIm = [char(sortedImages(iMask,1)) filesep char(sortedImages(iMask,2)),...;
%             num2str(sortednum(iMask)) char(sortedImages(iMask,4))];
          fileNameIm = [listOfImages{iMask,2} filesep listOfImages{iMask,1}]; 

        img = double(imread(fileNameIm));
        imgFiltered = filterGauss2D(img,1);
        maxSignal = max(imgFiltered(:)); 
        minSignal = min(imgFiltered(:)); 
        %normalize image between 0 and 1
imgFilteredNorm = (imgFiltered - minSignal)/(maxSignal - minSignal);
      if strcmpi(threshMethod,'Rosin')
      
        [level]  = thresholdRosin(imgFilteredNorm);
      elseif strcmpi(threshMethod,'MinMax')
   [level] = thresholdFluorescenceImage(imgFilteredNorm);
      else 
          [level] = graythresh(imgFilteredNorm);
      end 
        
        mask = imgFilteredNorm>=level; % get mask
        
        
        % clean it up
        mask = imfill(mask,'holes');
        % get largest CC
        mask = logical(getLargestCC(mask));
        % mask = imopen(mask,strel('disk',3));
        mask = logical(getLargestCC(mask));
        % maskDir = [projList{iProj} filesep 'masks'];
        maskDir = [upOne filesep 'roi_1' filesep 'masks'];
        
        if~ isdir(maskDir)
            %rmdir(maskDir,'s')
            mkdir(maskDir)
        end
        
        if ~isdir([maskDir filesep 'overlays'])
%            rmdir([maskDir filesep 'overlays'],'s'); 
            
            mkdir([maskDir filesep 'overlays']);
        end
        roiYX = bwboundaries(mask);
        
        figure('visible','off')
        imshow(-img,[]);
        hold on
        cellfun(@(x) plot(x(:,2),x(:,1),'r'),roiYX);
        saveas(gcf,[maskDir filesep 'overlays' filesep threshMethod 'Overlay' num2str(iMask, '%03d') '.tif' ]);
        
        imwrite(mask,[ maskDir filesep 'mask' threshMethod num2str(iMask,'%03d') '.tif']);
        
        if saveSumMask ==1
            if iMask == 1
                [ny,nx] = size(img);
                masks = zeros(ny,nx,iMask);
            else
                masks(:,:,iMask)= mask;
            end
        end % if saveSumMask
        clear mask
        close gcf
    end % for iMask
    
if saveSumMask == 1
    sumMask = sum(masks,3);
    [ny,nx] = size(img); 
    setFigure(nx,ny,'off'); 
    imagesc(sumMask)
    colorbar
    saveas(gcf,[projList(iProj).anDir filesep 'sumMask.tif']);
    close gcf
end
end % for iProj

% add small option to make a sum mask


end