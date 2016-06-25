function [ projList ] = mitoticPerformManualFunctions(projList)

%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes her
if (nargin<1 || isempty(projList))
 [filename,path] =  uigetfile(pwd,'Please Select a projList.mat File to Load');
  load([path filesep filename]); 
end 
%% first Make kinetochore mask 
for iProj = 1:numel(projList)
[listOfImages] = searchFiles('.tif',[],projList(iProj).imDir,0);
fileNameIm = [char(listOfImages(1,2)) filesep char(listOfImages(1,1))];
img = double(imread(fileNameIm))./((2^12)-1);
if exist([projList(iProj).anDir filesep 'roiMaskKinetochore.tif'],'file'); 
maskInternal =  logical(imread(([projList(iProj).anDir filesep 'roiMaskKinetochore.tif']))); 
else 

h = msgbox('please make a spindle roiMask'); 
uiwait(h); 
figure; 
 [maskInternal,polyXcoord,polyYcoord]=roipoly(img);
 up1 = getFilenameBody(char(listOfImages(1,2)));
 imwrite(maskInternal,[up1 filesep 'roi_1' filesep 'roiMaskKinetochore.tif'],'tif'); 
  close(gcf)       
end 
%% test for bad masks 
%maskExternal = logical(imread([projList(iProj).anDir filesep 'roiMaskSeg.tif']));
%upOne = upDirectory(projList(iProj).anDir,1); 
if exist([projList(iProj).anDir filesep 'badMaskFlag.mat']) == 0; 
maskExternal = double(imread([projList(iProj).anDir filesep 'masks' filesep 'maskRosin001.tif']));
  mask = maskExternal - maskInternal; 
%mask = maskExternal;
figure; 
imshow(img,[]); 
hold on
     roiYX = bwboundaries(mask);
     cellfun(@(x) plot(x(:,2),x(:,1),'y'),roiYX);
     reply1 = questdlg('Is this a good mask?') ; 
     if strcmpi(reply1,'no') 
         
         badMaskFlag = 0; % 
     elseif strcmpi(reply1,'yes')
         badMaskFlag = 1; 
     end 
      projList(iProj).badMaskFlag = badMaskFlag; 
         save([projList(iProj).anDir filesep 'badMaskFlag.mat'],'badMaskFlag'); % kindof stupid data struct but whatever
else 
    display(['Mask has been checked for ' projList(iProj).anDir 'Skipping']) ; 
end 
%% start pole documentation 
if exist([projList(iProj).anDir filesep 'poleInfo.mat'])==0; 
poleMaskTotal = zeros(size(img)); 
clickPole = 1; 
ptCount =1;
while clickPole == 1
    
    imshow(img,[]);
    hold on 
    roiYX = bwboundaries(mask);
    cellfun(@(x) plot(x(:,2),x(:,1),'y'),roiYX);
    spy(poleMaskTotal,'r',50); % plot the poles thus far ; 
    reply2 = questdlg('Document New Centrosome?');
   
   
    if strcmpi(reply2,'yes')
%         imshow(img,[])
%         hold on 
%         cellfun(@(x) plot(x(:,2),x(:,1),'y'),roiYX);
%         spy(poleMaskTotal,'r',50); % plot the poles thus far ;
        h=impoint;
        position = wait(h);
        idx = sub2ind(size(img),round(position(2)),round(position(1)));
        coords(ptCount,2) = position(2);
        coords(ptCount,1) = position(1);
        temp = zeros(size(img)); 
        temp(idx) = 1; 
        poleMasksInd(:,:,ptCount) = temp ;
        poleMaskTotal(idx) = 1; % add points to mask 
        ptCount = ptCount+1; 
    elseif (strcmpi(reply2,'no') && ptCount ==1) 
        coords(ptCount,2) = nan; 
        coords(ptCount,1) = nan; 
        poleInfo.coords = coords; 
        poleInfo.numPoles = 0 ; 
        save([projList(iProj).anDir filesep 'poleInfo.mat'],'poleInfo');
        clickPole = 0; 
    else
        % finish up and save
        if ~isempty(coords)
          poleInfo.coords = coords; % save coordinations of  poles 
          poleInfo.poleMasksInd = poleMasksInd; 
          poleInfo.poleMaskTotal = poleMaskTotal; 
          numPoles = length(coords(:,1)); 
            poleInfo.numPoles = numPoles; % save info regarding the number of poles
            if numPoles > 2 
                poleInfo.mult = 1; 
            else poleInfo.mult =0; 
            end 
        
        save([projList(iProj).anDir filesep 'poleInfo.mat'],'poleInfo');    
                
        end % isemtpy coords
        clear coords poleMasksInd poleMaskTotal poleInfo 
      
        clickPole = 0; % stop iteration
    end % if strcmpi
        
end % while 
else 
    display(['Centrosomes have been documented for ' projList(iProj).anDir ': Skipping'])
end % if exist
% 
% [~,name,num] = upDirectory(projList(iProj).anDir,2); 
% % check DIC 
% [filename,path] = uigetfile('.tif',['Select a DIC file for ' name num]); 
% 
% img = double(imread([path filesep filename])); 
%    imshow(img,[]) ; 
%    
%    
%    list{1} = 'Anaphase'; 
%    list{2} = 'Metaphase'; 
%    ok = 0; 
%    while ok ==0       
%    [reply1,ok] = listdlg('ListString',list) ; 
%    end 
%    cellStage = list{reply1}; 
% save([projList(iProj).anDir filesep 'cellStage.mat'],'cellStage')% for now just save as a simple string. 
% 
% 
% 
  close(gcf)
end        
   
   
       
         
end

     
         
        

 


    



