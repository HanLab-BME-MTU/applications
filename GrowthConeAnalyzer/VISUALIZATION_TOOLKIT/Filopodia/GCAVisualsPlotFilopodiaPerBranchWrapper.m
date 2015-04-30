function [ output_args ] = GCAVisualsPlotFilopodiaPerBranchWrapper(projList)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for iProj = 1:numel(projList)
    load([projList{iProj} filesep 'movieData.mat']); 
    load([MD.outputDirectory_ filesep 'filopodia_reconstruct' filesep 'Filopodia_Reconstruct_Channel_1'...
        filesep 'analInfoTestSave.mat']); 
    
    saveDir = [MD.outputDirectory_ filesep 'OUTPUT_MOVIES' filesep 'BranchGroups']; 
   % imgNames = MD.getImageFilenames{1}; 
    if ~isdir(saveDir) 
        mkdir(saveDir)
    end 
    
    for iFrame = 1: numel(analInfo)-1 
        filoInfo = analInfo(iFrame).filoInfo; 
        img = double(imread([MD.getChannelPaths{1} filesep MD.getImageFileNames{1}{iFrame}])); 
       
        [ny,nx] = size(img); 
        setFigure(nx,ny,'off'); 
        imshow(-img,[]); 
        hold on 
        GCAVisualsPlotFilopodiaPerBranchGroup(filoInfo,MD.imSize_); 
        saveas(gcf,[saveDir filesep num2str(iFrame,'%03d') '.png']); 
        close gcf
        
    end 
end 

end

