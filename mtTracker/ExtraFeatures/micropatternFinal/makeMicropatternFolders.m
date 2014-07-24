function [ groupListSubRoiIndProj, subRoiDirs ] = makeMicropatternFolders(currentProjDir,projNum,roiSet,maskParams,groupName,toExtract)
%
% For each roiSet makes the respective subDirs and saves a combinedMaskFig
% to that folder
% 
%
%
% 
% start
[imL imW] = size(roiSet(:,:,1));  
numRois = numel(roiSet(1,1,:)); 
subRoiDirs = cell(numRois,1);
groupListSubRoiIndProj = cell(numRois,2); 
cMap=hsv(9);

% file management: start at level of whole cell roi
upOne = getFilenameBody(currentProjDir); 
upTwo = getFilenameBody(upOne); 
%upThree = getFilenameBody(upTwo); 
s1 = num2str(maskParams.numWindows); 
s2 = num2str(maskParams.windowSize); 
if maskParams.subRegions == 1
    s3 = 'subRegions';
else 
    s3 = '';
end ; 
s4 = toExtract; 
saveDir = [currentProjDir filesep 'SUBROIS_' s1 '_' s2 '_umWindows_' char(s3) '_' char(s4)]; 
collectDir = [upTwo filesep 'subRoiMaskSummary' filesep s1 '_' s2 '_umWindows_' char(s3)];   
if ~isdir([upTwo filesep 'subRoiMaskSummary' ])
    mkdir([upTwo filesep 'subRoiMaskSummary']); 
end 
if ~isdir(collectDir) 
    mkdir(collectDir)
end 

% make the subRoi directory 
    if ~isdir(saveDir) 
        mkdir(saveDir); 
    end 
   


% get representive image for plotting
imDir = [upOne filesep 'images']; 
[listOfImages] = searchFiles('.tif',[],imDir,0);
firstImage = [char(listOfImages(1,2)) filesep char(listOfImages(1,1))];
img = double(imread(firstImage));    
img1=(img-min(img(:)))./(max(img(:))-min(img(:)));
img=repmat(img1,[1,1,3]);


for iRegion = 1:numRois    
      if iRegion == 1
            % get whole cell mask
            roiMask=imread([currentProjDir filesep 'roiMaskSeg.tif']);
            roiMask = logical(roiMask); 
            roiYX=load([currentProjDir filesep 'roiYXSeg.mat'],'roiYX'); roiYX=roiYX.roiYX;

            % plot cell outline on image
            compImage=figure('Visible','off');
            imshow(repmat(roiMask,[1 1 3]).*img)
            hold on
            plot(roiYX(:,2),roiYX(:,1),'w');
      end 
           
 
    
   
    
    roiDir = [saveDir filesep 'sub_' num2str(iRegion)]; 
    if ~isdir(roiDir) 
       mkdir(roiDir)
    end 
    
    maskSubRoi = roiSet(:,:,iRegion); 
    
    if sum(maskSubRoi(:)) == 0 % bad mask need to break loop 
     groupListSubRoiIndProj = []; 
     % flag for bad mask
    
     
        return
    end 
        
        
     [y1,x1]=ind2sub([imL,imW],find(maskSubRoi,1)); % first pixel on boundary
     roiYXSubRoi = bwtraceboundary(maskSubRoi,[y1,x1],'N'); % get all pixels on boundary
   
     % save subRegion Masks
            imwrite(maskSubRoi,[roiDir filesep 'roiMask.tif']); 
            save([roiDir filesep 'roiYXSubRoi'],'roiYXSubRoi'); 
    
    % store directories for output
    subRoiDirs{iRegion} = roiDir; 
    
    groupListSubRoiIndProj{iRegion,2} = roiDir; 
    groupListSubRoiIndProj{iRegion,1} = [char(groupName) '_sub' num2str(iRegion)]; 
    
    if ~isdir([roiDir filesep 'feat']) 
        mkdir([roiDir filesep 'feat']) 
    end ; 
    movieInfo1  = [currentProjDir filesep 'feat' filesep 'movieInfo.mat'];  
    movieInfo2 = [roiDir filesep 'feat' filesep 'movieInfo.mat']; 
    copyfile(movieInfo1,movieInfo2); 
        
    
    
    if ~isdir([roiDir filesep 'meta'])
        mkdir([roiDir filesep 'meta'])
    end 
    
     plot(roiYXSubRoi(:,2),roiYXSubRoi(:,1),'Color',cMap(max(1,mod(iRegion,10)),:));

        % make weighted mask using distance transform to find position
        % where text should go
        weightedRoi=bwdist(~maskSubRoi);
        [r,c]=find(weightedRoi==max(weightedRoi(:)));
        text(c(1),r(1),num2str(iRegion),'color','r','fontsize',14);

end 
      % save composite image
     indxStr1 = num2str(numRois);
     saveas(compImage,[saveDir filesep 'subROIs_' indxStr1 '.fig'])
     saveas(compImage,[saveDir filesep 'subROIs_' indxStr1 '.tif'])
     saveas(compImage,[collectDir filesep 'subROIs_' num2str(projNum) '.tif']); 
     close(compImage)
     save([collectDir filesep 'roiSet','roiSet']); 
end % 






