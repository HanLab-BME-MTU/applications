function this_MD = vimscreen_forceMDAddProcessSave(this_MD)
% Function of single image filament segmentation with input this_MD from other
%               successfully segmented movie for the parameters in this_MD

% Input:      this_MD:              the movieData for data to be segmented
% Created 2015.1 by Liya Ding, Matlab R2013a


%% get the image dir and make a this_MD just for this image
ROOT_DIR = this_MD.outputDirectory_;

nChannel = numel(this_MD.channels_);
nProcess = numel(this_MD.processes_);
nFrame = this_MD.nFrames_;
nPackage = length(this_MD.packages_);

default_parameter_cell= vimscreen_default_parameters(this_MD);
set_parameter_cell = default_parameter_cell;

% look for exisitng filament package
indexFilamentPackage = 0;
for i = 1 : nPackage
    if(strcmp(this_MD.packages_{i}.getName,'FilamentAnalysis')==1)
        indexFilamentPackage = i;
        break;
    end
end

nProcess = numel(this_MD.processes_);
nPackage = length(this_MD.packages_);

% if there is no filament analysis package,
% add filament analysis package
if(indexFilamentPackage == 0)
    this_MD.addPackage(FilamentAnalysisPackage(this_MD));
    mkdir([this_MD.outputDirectory_,filesep,'FilamentAnalysisPackage']);
    indexFilamentPackage = 0;
    for i = 1 : nPackage
        if(strcmp(this_MD.packages_{i}.getName,'FilamentAnalysis')==1)
            indexFilamentPackage = i;
            break;
        end
    end
end

% start process by process

    %%     % check if there is each of the process
    indexCellSegProcess = 0;
    for i = 1 : numel(this_MD.processes_)
        if(strcmp(this_MD.processes_{i}.getName,'Thresholding')==1)
            indexCellSegProcess = i;
            break;
        end
    end
    
    indexCellRefineProcess = 0;
    for i = 1 : numel(this_MD.processes_)
        if(strcmp(this_MD.processes_{i}.getName,'Mask Refinement')==1)
            indexCellRefineProcess = i;
            break;
        end
    end
    
    indexFlattenProcess = 0;
    for i = 1 : numel(this_MD.processes_)
        if(strcmp(this_MD.processes_{i}.getName,'Image Flatten')==1)
            indexFlattenProcess = i;
            break;
        end
    end
    
    indexSteerabeleProcess = 0;
    for i = 1 : numel(this_MD.processes_)
        if(strcmp(this_MD.processes_{i}.getName,'Steerable filtering')==1)
            indexSteerabeleProcess = i;
            break;
        end
    end
    
    indexFilamentSegmentationProcess = 0;
    for i = 1 : numel(this_MD.processes_)
        if(strcmp(this_MD.processes_{i}.getName,'Filament Segmentation')==1)
            indexFilamentSegmentationProcess = i;
            break;
        end
    end
    
    %%  % 1 threshold
    
    %%   % when there is no input this_MD, and no existing processes
    % use the default for all of the processes
    
     default_Params = set_parameter_cell{1};
   
    if(indexCellSegProcess==0)
        this_MD.addProcess(ThresholdProcess(this_MD,'FilamentAnalysisPackage',default_Params));
    end
    
    
    %%  % 2 mask refine
    % with the default of mask refinement
 default_Params = set_parameter_cell{2};
     
    if(indexCellRefineProcess==0)
        this_MD.addProcess(MaskRefinementProcess(this_MD,'FilamentAnalysisPackage',default_Params));
    end
    
    
    %%   % 3 image flatten
    
 default_Params = set_parameter_cell{3};
     
    if(indexFlattenProcess==0)
        this_MD.addProcess(ImageFlattenProcess(this_MD,'funParams',default_Params));
    end
    
   
    %% % 4 steerable filter
    default_Params = set_parameter_cell{4};
 
    
    if(indexSteerabeleProcess==0)
         this_MD.addProcess(SteerableFilteringProcess(this_MD,'funParams',default_Params));
    end
    
     
    %%
    % 5 filament segmentation
    
    default_Params = set_parameter_cell{5};
 
    if(indexFilamentSegmentationProcess==0)
        this_MD.addProcess(FilamentSegmentationProcess(this_MD,'funParams',default_Params));
    end
    
    %% after adding processes, set the paths
        %%     % check if there is each of the process
    indexCellSegProcess = 0;
    for i = 1 : numel(this_MD.processes_)
        if(strcmp(this_MD.processes_{i}.getName,'Thresholding')==1)
            indexCellSegProcess = i;
            break;
        end
    end
    
    indexCellRefineProcess = 0;
    for i = 1 : numel(this_MD.processes_)
        if(strcmp(this_MD.processes_{i}.getName,'Mask Refinement')==1)
            indexCellRefineProcess = i;
            break;
        end
    end
    
    indexFlattenProcess = 0;
    for i = 1 : numel(this_MD.processes_)
        if(strcmp(this_MD.processes_{i}.getName,'Image Flatten')==1)
            indexFlattenProcess = i;
            break;
        end
    end
    
    indexSteerabeleProcess = 0;
    for i = 1 : numel(this_MD.processes_)
        if(strcmp(this_MD.processes_{i}.getName,'Steerable filtering')==1)
            indexSteerabeleProcess = i;
            break;
        end
    end
    
    indexFilamentSegmentationProcess = 0;
    for i = 1 : numel(this_MD.processes_)
        if(strcmp(this_MD.processes_{i}.getName,'Filament Segmentation')==1)
            indexFilamentSegmentationProcess = i;
            break;
        end
    end
    
     if(isempty(this_MD.processes_{indexCellSegProcess}.outFilePaths_{1}))
        this_MD.processes_{indexCellSegProcess}.setOutFilePaths({[this_MD.outputDirectory_,filesep,'FilamentAnalysisPackage',filesep,'masks',filesep,'Channel1'],...
        [this_MD.outputDirectory_,filesep,'FilamentAnalysisPackage',filesep,'masks',filesep,'Channel2'],...
        [this_MD.outputDirectory_,filesep,'FilamentAnalysisPackage',filesep,'masks',filesep,'Channel3']});
     end
     
     if(isempty(this_MD.processes_{indexCellRefineProcess}.outFilePaths_{1}))
        this_MD.processes_{indexCellRefineProcess}.setOutFilePaths({[this_MD.outputDirectory_,filesep,'FilamentAnalysisPackage',filesep,'refined_masks',filesep,'Channel1'],...
        [this_MD.outputDirectory_,filesep,'FilamentAnalysisPackage',filesep,'refined_masks',filesep,'Channel2'],...
        [this_MD.outputDirectory_,filesep,'FilamentAnalysisPackage',filesep,'refined_masks',filesep,'Channel3']});
    end
     
    if(isempty(this_MD.processes_{indexFlattenProcess}.outFilePaths_{1}))
        this_MD.processes_{indexFlattenProcess}.setOutFilePaths({[this_MD.outputDirectory_,filesep,'FilamentAnalysisPackage',filesep,'ImageFlatten',filesep,'Channel1'],...
        [this_MD.outputDirectory_,filesep,'FilamentAnalysisPackage',filesep,'ImageFlatten',filesep,'Channel2'],...
        [this_MD.outputDirectory_,filesep,'FilamentAnalysisPackage',filesep,'ImageFlatten',filesep,'Channel3']});
    end
    
    if(isempty(this_MD.processes_{indexSteerabeleProcess}.outFilePaths_{1}))
        this_MD.processes_{indexSteerabeleProcess}.setOutFilePaths({[this_MD.outputDirectory_,filesep,'FilamentAnalysisPackage',filesep,'SteerableFiltering',filesep,'Channel1'],...
        [this_MD.outputDirectory_,filesep,'FilamentAnalysisPackage',filesep,'SteerableFiltering',filesep,'Channel2'],...
        ''});
    end    
     
    if(isempty(this_MD.processes_{indexFilamentSegmentationProcess}.outFilePaths_{1}))
        this_MD.processes_{indexFilamentSegmentationProcess}.setOutFilePaths({[this_MD.outputDirectory_,filesep,'FilamentAnalysisPackage',filesep,'FilamentSegmentation',filesep,'Channel1'],...
        [this_MD.outputDirectory_,filesep,'FilamentAnalysisPackage',filesep,'FilamentSegmentation',filesep,'Channel2'],...
        ''});
    end
    
    
    this_MD.save();
    
