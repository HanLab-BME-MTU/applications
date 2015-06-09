function [ output_args ] = helperMoveLocalThreshPatchSizeParam(projList )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


for iProj =1 :numel(projList)
    oldFile = [projList{iProj} filesep 'ANALYSIS' filesep 'neurite_body_masks' filesep ...
        'Neurite_Body_Masks_Channel_1' filesep 'params.mat'];
    if exist(oldFile)~=0
        load(oldFile);
        patchSize = p.patchSize;
        
        newDir =  [projList{iProj} filesep 'GrowthConeAnalyzer'  filesep...
            'SegmentationPackage' filesep 'StepsToReconstruct' filesep 'III_veilStem_reconstruction' ...
            filesep 'Channel_1' filesep 'LocalThreshPatchSizeTest'];
        if ~isdir(newDir)
            mkdir(newDir)
        end
        
        
        save([newDir filesep 'manualPatchSizeSelect.mat'],'patchSize');
        
    end
    
end
end






