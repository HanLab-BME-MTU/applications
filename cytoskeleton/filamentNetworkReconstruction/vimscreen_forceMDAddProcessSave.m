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
end

% start process by process

    %%     % check if there is each of the process
    indexCellSegProcess = 0;
    for i = 1 : nProcess
        if(strcmp(this_MD.processes_{i}.getName,'Thresholding')==1)
            indexCellSegProcess = i;
            break;
        end
    end
    
    indexCellRefineProcess = 0;
    for i = 1 : nProcess
        if(strcmp(this_MD.processes_{i}.getName,'Mask Refinement')==1)
            indexCellRefineProcess = i;
            break;
        end
    end
    
    indexFlattenProcess = 0;
    for i = 1 : nProcess
        if(strcmp(this_MD.processes_{i}.getName,'Image Flatten')==1)
            indexFlattenProcess = i;
            break;
        end
    end
    
    indexSteerabeleProcess = 0;
    for i = 1 : nProcess
        if(strcmp(this_MD.processes_{i}.getName,'Steerable filtering')==1)
            indexSteerabeleProcess = i;
            break;
        end
    end
    
    indexFilamentSegmentationProcess = 0;
    for i = 1 : nProcess
        if(strcmp(this_MD.processes_{i}.getName,'Filament Segmentation')==1)
            indexFilamentSegmentationProcess = i;
            break;
        end
    end
    
    %%  % 1 threshold
    
    %%   % when there is no input this_MD, and no existing processes
    % use the default for all of the processes
    
    default_Params=[];
    
    default_Params.ChannelIndex = 3;
    default_Params.GaussFilterSigma=2;
    default_Params.ThresholdValue=[];
    default_Params.OutputDirectory = [this_MD.outputDirectory_,filesep,'FilamentAnalysisPackage',filesep,'thres'];
    default_Params.MethodIndx=4;
    default_Params.ProcessIndex=[];
    default_Params.MaxJump=10;
    default_Params.BatchMode=0;
    
    
    if(indexCellSegProcess==0)
        this_MD.addProcess(ThresholdProcess(this_MD,'FilamentAnalysisPackage',default_Params));
    end
    
    
    %%  % 2 mask refine
    % with the default of mask refinement
    
    default_Params=[];
    
    default_Params.ChannelIndex = 3;
    default_Params.MaskCleanUp = true;
    default_Params.MinimumSize = 10;
    default_Params.ClosureRadius = 3;
    default_Params.ObjectNumber = 50;
    default_Params.FillHoles = true;
    default_Params.EdgeRefinement = false; %This off by default because it sort of sucks, and is slow.
    default_Params.SegProcessIndex = []; %No default.
    default_Params.OutputDirectory = [this_MD.outputDirectory_,filesep,'FilamentAnalysisPackage',filesep,'refined_masks'];
    default_Params.OpeningRadius = 0;
    default_Params.SuppressBorder = true;
    default_Params.MaxEdgeAdjust = []; %Use refineMaskEdges.m function defaults for these settings
    default_Params.MaxEdgeGap = [];
    default_Params.PreEdgeGrow = [];
    default_Params.BatchMode = false;
    
    if(indexCellRefineProcess==0)
        this_MD.addProcess(MaskRefinementProcess(this_MD,'FilamentAnalysisPackage',default_Params));
    end
    
    
    %%   % 3 image flatten
    
    default_Params=[];
    
    default_Params.ChannelIndex= [1 2];
    default_Params.GaussFilterSigma= 0.2000;
    default_Params.method_ind= 3;
    default_Params.imageflattening_mode= 2;
    default_Params.TimeFilterSigma= 0;
    default_Params.Sub_Sample_Num= 1;
    default_Params.outputDir= [this_MD.outputDirectory_,filesep,'FilamentAnalysisPackage',filesep,'ImageFlatten'];
    default_Params.stat.low_005_percentile= 0;
    default_Params.stat.high_995_percentile=65535;
    default_Params.stat.center_value_int=300;
    
    if(indexFlattenProcess==0)
        this_MD.addProcess(ImageFlattenProcess(this_MD,'funParams',default_Params));
    end
    
   
    %% % 4 steerable filter
    default_Params=[];
    
    default_Params.ChannelIndex= [1 2];
    default_Params.BaseSteerableFilterSigma= 1;
    default_Params.Levelsofsteerablefilters= 2;
    default_Params.Sub_Sample_Num= 1;
    default_Params.ImageFlattenFlag= 2;
    
    
    if(indexSteerabeleProcess==0)
         this_MD.addProcess(SteerableFilteringProcess(this_MD,'funParams',default_Params));
    end
    
     
    %%
    % 5 filament segmentation
    
    default_Params.ChannelIndex = [1 2];
    default_Params.StPace_Size = [3 3 3];
    default_Params.StPatch_Size = [21 21 21];
    default_Params.st_lowerbound_localthresholding = [90 90 90];
    default_Params.IntPace_Size = [3 3 3];
    default_Params.IntPatch_Size = [21 21 21];
    default_Params.int_lowerbound_localthresholding = [90 90 90];
    default_Params.Combine_Way{1} = 'geo_based_no_GM';
    default_Params.Combine_Way{2} = 'geo_based_no_GM';
    default_Params.Combine_Way{3} = 'geo_based_no_GM';
    default_Params.Cell_Mask_ind = [5 5 1];
    default_Params.Whole_movie_ind = [1 1 2];
    default_Params.Whole_movie_stat_cell{1} = [];
    default_Params.Whole_movie_stat_cell{2} = [];
    default_Params.Whole_movie_stat_cell{3} = [];
    default_Params.Rerun_WholeMovie = 0;
    default_Params.VIF_Outgrowth_Flag = 0;
    default_Params.Sub_Sample_Num = 1;
    default_Params.F_classifier{1} = [];
    default_Params.F_classifier{2} = [];
    default_Params.F_classifier{3} = [];
    default_Params.Classifier_Type_ind = [1 1 1];
    default_Params.training_sample_number = [30 30 30];
    default_Params.LengthThreshold = [4 4 4];
    default_Params.CurvatureThreshold = [0.1 0.1 0.1];
    default_Params.CoefAlpha = [1.8 1.8 1.8];
    default_Params.IternationNumber = [0 0 2];
    default_Params.CannyHigherThreshold = [80 80 80];
    default_Params.CannyLowerThreshold = [80 80 80];
    default_Params.nofiguredisruption = 1;
    default_Params.savestepfigures = 0;
    default_Params.showdetailmessages = 0;
    default_Params.channel_specific = 0;
    
    if(indexFilamentSegmentationProcess==0)
        this_MD.addProcess(FilamentSegmentationProcess(this_MD,'funParams',default_Params));
    end
    
    
    nProcess = numel(this_MD.processes_);
    nPackage = length(this_MD.packages_);


    %% after adding processes, set the paths
        %%     % check if there is each of the process
    indexCellSegProcess = 0;
    for i = 1 : nProcess
        if(strcmp(this_MD.processes_{i}.getName,'Thresholding')==1)
            indexCellSegProcess = i;
            break;
        end
    end
    
    indexCellRefineProcess = 0;
    for i = 1 : nProcess
        if(strcmp(this_MD.processes_{i}.getName,'Mask Refinement')==1)
            indexCellRefineProcess = i;
            break;
        end
    end
    
    indexFlattenProcess = 0;
    for i = 1 : nProcess
        if(strcmp(this_MD.processes_{i}.getName,'Image Flatten')==1)
            indexFlattenProcess = i;
            break;
        end
    end
    
    indexSteerabeleProcess = 0;
    for i = 1 : nProcess
        if(strcmp(this_MD.processes_{i}.getName,'Steerable filtering')==1)
            indexSteerabeleProcess = i;
            break;
        end
    end
    
    indexFilamentSegmentationProcess = 0;
    for i = 1 : nProcess
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
    
