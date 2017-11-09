function [ output_args ] = GCAfindNeuriteLengthWrapper(projList)
%
for iProj = 1: numel(projList)
    load([projList{iProj} filesep 'ANALYSIS' filesep 'movieData.mat']);
    % check for veilstem refinement folder 
    veilStemRefine = [MD.outputDirectory_ filesep 'neurite_veilStem_refinements'...
        filesep 'analInfoWithScaleSave.mat'];
    if exist(veilStemRefine,'file') ~=0
        display('Loading veil/stem masks from neurite_veilStem_refinements folder');
        load(veilStemRefine);
        
        skip = 0;% don't skip the calc
    else % check for other veil Stem path
        veilStemPath= [MD.outputDirectory_ filesep 'neurite_body_masks' filesep ...
            'Neurite_Body_Masks_Channel_1' filesep 'analInfoTestSave.mat'];
        
        
        if exist(veilStemPath,'file') ~=0
            display('Loading Veil/Stem Masks from neurite_body_masks folder');
            load(veilStemPath)
            
            skip = 0;% don't skip
        else
            display(['No Veil/Stem Analysis Found for ' projList{iProj} ...
                'Skipping']);
            skip = 1;
        end
    end
    
    if exist([MD.outputDirectory_ filesep 'PARAMETER_EXTRACTION_20150316\'...
            'GlobalFunctional\neurite_outgrowth_measurements' filesep 'neuriteLengthOutput.mat']) ~=0
        display(['NeuriteLength.mat Already Run for ' projList{iProj} ': Skipping']);
        skip = 1; 
    else
        
        if skip ==0 
        
        imDir = MD.getChannelPaths{1};
        saveDir = MD.outputDirectory_;
        GCAfindNeuriteLength20150307(analInfo,saveDir,imDir,'k');
        end 
    end
    
    
    
end



end

