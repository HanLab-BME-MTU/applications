function [ output_args ] = helperMoveLocalThreshPatchSizeParam(projList )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


for iProj =1 :size(projList,1)
    oldFile = [projList{iProj,1} filesep 'ANALYSIS' filesep 'neurite_body_masks' filesep ...
        'Neurite_Body_Masks_Channel_1' filesep 'params.mat'];
    if exist(oldFile)~=0
        load(oldFile);
        patchSize = p.patchSize;
    else 
       patchSize=  inputdlg(['Input PatchSize for' projList{iProj,2}] ,'PatchSize',1,{'75'}); 
       patchSize = str2double(patchSize); 
    end
       
        
        newDir =  [projList{iProj,1} filesep 'GrowthConeAnalyzer'  filesep...
            'SegmentationPackage' filesep 'StepsToReconstruct' filesep 'III_veilStem_reconstruction' ...
            filesep 'Channel_1' filesep 'LocalThreshPatchSizeTest'];
        if ~isdir(newDir)
            mkdir(newDir)
        end
        
        
        save([newDir filesep 'manualPatchSizeSelect.mat'],'patchSize');
        
    
    
end
end






