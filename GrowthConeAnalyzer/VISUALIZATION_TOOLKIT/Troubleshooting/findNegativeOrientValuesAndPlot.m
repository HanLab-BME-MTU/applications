function [ output_args ] = findNegativeOrientValuesAndPlot( movieData )
%UNTITLED4 Summary of this function goes he
%   Detailed explanation goes here
load([movieData.outputDirectory_ filesep 'filopodia_fits' filesep ...
    'Filopodia_Fits_Channel_1' filesep 'analInfoTestSave.mat']); 
% get the frames with the bad filo orientation 
badOrientFlag = arrayfun(@(x)  sum(vertcat(analInfo(x).filoInfo.orientation)<0)>0,1:length(analInfo)-1) ;

idxFramesToCheck = find(badOrientFlag);

for iFrame = 1:length(idxFramesToCheck) 
    % load image 
    cFrame = idxFramesToCheck(iFrame) ; 
   img = double(imread([movieData.getChannelPaths{1} filesep movieData.getImageFileNames{1}{cFrame}])); 
   [ny,nx] = size(img); 
   setFigure(nx,ny,'on'); 
   imshow(-img,[]); 
   hold on 
   
    
    % get the problem filo 
    
    filoInfo = analInfo(cFrame).filoInfo; 
    orients = vertcat(filoInfo(:).orientation); 
    filoInfoIdx = filoInfo(orients<0); 
    mask =  analInfo(cFrame).masks.neuriteEdge; 
    roiYX = bwboundaries(mask); 
    cellfun(@(x) plot(x(:,2),x(:,1),'y'),roiYX); 
    hold on 
    GCAVisualsMakeOverlaysFilopodia(filoInfoIdx,[ny,nx],1,1,'r',1);
    GCAVisualsPlotFilopodiaIDs(filoInfo,find(orients<0),'k'); 
    
    
  saveas(gcf,['problemFrame' num2str(cFrame,'%03d') '.fig']);   
    close gcf
end 


end

