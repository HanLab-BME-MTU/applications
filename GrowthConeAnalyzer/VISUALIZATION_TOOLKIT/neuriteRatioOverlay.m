function [ output_args ] = neuriteRatioOverlay(ratioDir,maskDir)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% get the limits for the movie 
listOfRatFiles = searchFiles('.mat',[],ratioDir,0); 
listOfMasks = searchFiles('.tif',[],maskDir,0); 
upOne = upDirectory(ratioDir,1); 

saveDir = [upOne filesep 'ratio_Movie']; 
if ~isdir(saveDir) 
    mkdir(saveDir) 
end 

for iFrame = 1:length(listOfRatFiles) 
    load([listOfRatFiles{iFrame,2} filesep listOfRatFiles{iFrame,1}]); 
 mask=    logical(imread([listOfMasks{iFrame,2} filesep listOfMasks{iFrame,1}])); 
 values{iFrame} = currRatio(mask); 
 
end 

% calculate clims for whole movie 
climMax = prctile(vertcat(values{:}),95); 
climMin = prctile(vertcat(values{:}),5); 
clims = [climMin, climMax]; 
for iFrame = 1:length(listOfRatFiles) 
     load([listOfRatFiles{iFrame,2} filesep listOfRatFiles{iFrame,1}]); 
     mask=    double(imread([listOfMasks{iFrame,2} filesep listOfMasks{iFrame,1}])); 
     ratImg = currRatio.*mask; 
     [ny,nx] = size(mask); 
     setFigure(ny,nx,'on'); 
     
     h = imagesc(ratImg,clims); 
     alphaMask = true(size(ratImg)); 
   alphaMask(ratImg ==0 | isnan(ratImg)) = false;
   set(h, 'alphaData', alphaMask,'alphaDataMapping','none'); 
   
   
  
    
   colorbar
   saveas(gcf,[saveDir filesep 'Ratio' num2str(iFrame,'%03d') '.png']); 
   if iFrame == 3 
       saveas(gcf,[saveDir filesep 'Ratio' num2str(iFrame,'%03d') '.fig']); 
   end 
     close gcf
     
end 





end

