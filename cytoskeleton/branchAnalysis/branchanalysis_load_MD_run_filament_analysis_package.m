function branchanalysis_load_MD_run_filament_analysis_package(this_MD, Parameter_MD, varargin)
% Function of single image filament segmentation with input this_MD from other
%               successfully segmented movie for the parameters in this_MD

% Input:      this_MD:              the movieData for data to be segmented
%             Parameter_MD:    a loaded this_MD with good segmentation parameters
%                              if none, put [], so a default setting will be used
%                              with (1) Otsu with smoothing 1 (2) mask
%                              refine with 1 object, (3) image flatten with
%                              square, (4) steerable filter with [1 2], (5)
%                              segmentation with geo based alpha=2
%             whole_movie_filename: optional input, if given, the whole
%                                movie statistics will be loaded from this filename, which
%                                will be used in the filament segmentsation; by default,
%                                nothing is input, so whole movie stat will be calculated,
%                                well, in this case, without any meaningful effect

% Created 2014.12 by Liya Ding, Matlab R2012b

ip = inputParser;
ip.addRequired('this_MD', @(x) isa(x,'MovieData'));
ip.addRequired('Parameter_MD',@(x) isa(x,'MovieData') | isempty(x));
ip.addOptional('whole_movie_filename', ' ', @ischar);

ip.parse(this_MD,Parameter_MD,varargin{:});
whole_movie_filename= ip.Results.whole_movie_filename;

%% get the image dir and make a this_MD just for this image
ROOT_DIR = this_MD.outputDirectory_;

nChannel = numel(this_MD.channels_);
nProcess = numel(this_MD.processes_);
nFrame = this_MD.nFrames_;
nPackage = length(this_MD.packages_);

% look for exisitng filament package
indexSegmentPackage = 0;
for i = 1 : nPackage
    if(strcmp(this_MD.packages_{i}.getName,'Segmentation')==1)
        indexSegmentPackage = i;
        break;
    end
end


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


% if there is no segment package,
% add segment package
if(indexSegmentPackage == 0)
    this_MD.addPackage(SegmentationPackage(this_MD));
    mkdir([this_MD.outputDirectory_,filesep,'SegmentationPackage']);
end

% if there is no filament analysis package,
% add filament analysis package
if(indexFilamentPackage == 0)
    this_MD.addPackage(FilamentAnalysisPackage(this_MD));
    mkdir([this_MD.outputDirectory_,filesep,'FilamentAnalysisPackage']);
end

% start process by process

if(~isempty(Parameter_MD))
    %% % if there is input Parameter_MD
    
    for iPro =  1 : numel(Parameter_MD.processes_)
        % 1 threshold
        
        if(strcmp(Parameter_MD.processes_{iPro}.getName, 'Thresholding'))
            given_Params = Parameter_MD.processes_{iPro}.funParams_;
            given_Params.OutputDirectory = [this_MD.outputDirectory_,filesep,'FilamentAnalysisPackage',filesep,'thres'];
            
            this_MD.addProcess(ThresholdProcess(this_MD,'default_Params',given_Params));
            this_MD = thresholdMovie(this_MD,this_MD.processes_{iPro}.funParams_);
            
        else
            % 2 mask refine
            if(strcmp(Parameter_MD.processes_{iPro}.getName, 'Mask Refinement'))
                given_Params = Parameter_MD.processes_{iPro}.funParams_;
                given_Params.OutputDirectory = [this_MD.outputDirectory_,filesep,'FilamentAnalysisPackage',filesep,'maskrefine'];
                
                this_MD.addProcess(MaskRefinementProcess(this_MD,'default_Params',given_Params));
                this_MD = refineMovieMasks(this_MD,this_MD.processes_{iPro}.funParams_);
%                 
                
            else
                % 3 image flatten
                if(strcmp(Parameter_MD.processes_{iPro}.getName, 'Image Flatten'))
                    given_Params = Parameter_MD.processes_{iPro}.funParams_;
                    
                    this_MD.addProcess(ImageFlattenProcess(this_MD));
                    this_MD = image_flatten(this_MD,given_Params);
                    
                else
                    % 4 steerable filter
                    if(strcmp(Parameter_MD.processes_{iPro}.getName, 'Steerable filtering'))
                        given_Params = Parameter_MD.processes_{iPro}.funParams_;
                        this_MD.addProcess(SteerableFilteringProcess(this_MD,'default_Params',given_Params));
                        this_MD = steerable_filter_forprocess(this_MD,given_Params);
                        
                    else
                        % 5 filament segmentation
                        if(strcmp(Parameter_MD.processes_{iPro}.getName, 'Filament Segmentation'))
                            given_Params = Parameter_MD.processes_{iPro}.funParams_;
                            this_MD.addProcess(FilamentSegmentationProcess(this_MD,'default_Params',given_Params));
                            
                            % if user give the whole movie stat, use this
                            if(~ismember('whole_movie_filename',ip.UsingDefaults))
                                this_MD = filament_segmentation(this_MD,given_Params,whole_movie_filename);
                            else
                                this_MD = filament_segmentation(this_MD,given_Params);
                            end
                        end
                    end
                end
            end
        end
    end
else
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
    
    default_Params.ChannelIndex = [1 2];
    default_Params.GaussFilterSigma=2;
    default_Params.ThresholdValue=[];
    default_Params.OutputDirectory = [this_MD.outputDirectory_,filesep,'SegmentationPackage',filesep,'thres'];
    default_Params.MethodIndx=4;
    default_Params.ProcessIndex=[];
    default_Params.MaxJump=10;
    default_Params.BatchMode=0;
    
    
    if(indexCellSegProcess==0)
        this_MD.addProcess(ThresholdProcess(this_MD,'SegmentationPackage',default_Params));
    end
    
%      this_MD = thresholdMovie(this_MD,default_Params);
    
    %%  % 2 mask refine
    % with the default of mask refinement
    
    default_Params=[];
    
    default_Params.ChannelIndex = [1 2];
    default_Params.MaskCleanUp = true;
    default_Params.MinimumSize = 10;
    default_Params.ClosureRadius = 3;
    default_Params.ObjectNumber = 20;
    default_Params.FillHoles = true;
    default_Params.EdgeRefinement = false; %This off by default because it sort of sucks, and is slow.
    default_Params.SegProcessIndex = []; %No default.
    default_Params.OutputDirectory = [this_MD.outputDirectory_,filesep,'SegmentationPackage',filesep,'refined_masks'];
    default_Params.OpeningRadius = 0;
    default_Params.SuppressBorder = true;
    default_Params.MaxEdgeAdjust = []; %Use refineMaskEdges.m function defaults for these settings
    default_Params.MaxEdgeGap = [];
    default_Params.PreEdgeGrow = [];
    default_Params.BatchMode = false;
    
    if(indexCellRefineProcess==0)
        this_MD.addProcess(MaskRefinementProcess(this_MD,'SegmentationPackage',default_Params));
    end
    
%      this_MD = refineMovieMasks(this_MD,default_Params);
     
    %%   % 3 image flatten
    
    default_Params=[];
    
    default_Params.ChannelIndex= [2];
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
    
    this_MD = image_flatten(this_MD,default_Params);
    
    %% % 4 steerable filter
    default_Params=[];
    
    default_Params.ChannelIndex= 2;
    default_Params.BaseSteerableFilterSigma= 1;
    default_Params.Levelsofsteerablefilters= 2;
    default_Params.Sub_Sample_Num= 1;
    default_Params.ImageFlattenFlag= 2;
    
    
    if(indexSteerabeleProcess==0)
         this_MD.addProcess(SteerableFilteringProcess(this_MD,'funParams',default_Params));
    end
    
    this_MD = steerable_filter_forprocess_continue(this_MD,default_Params);
    
    %%
    % 5 filament segmentation
    
    default_Params.ChannelIndex = 2;
    default_Params.StPace_Size = [3 3];
    default_Params.StPatch_Size = [21 21];
    default_Params.st_lowerbound_localthresholding = [90 90] ;
    default_Params.IntPace_Size = [3 3];
    default_Params.IntPatch_Size = [21 21 ];
    default_Params.int_lowerbound_localthresholding = [90 90 ];
    default_Params.Combine_Way{1} = 'geo_based_GM';
    default_Params.Combine_Way{2} = 'geo_based_GM';
    default_Params.Cell_Mask_ind = [5 5];
    default_Params.Whole_movie_ind = [1 1];
    default_Params.Whole_movie_stat_cell{1} = [];
    default_Params.Whole_movie_stat_cell{2} = [];
    default_Params.Rerun_WholeMovie = 0;
    default_Params.VIF_Outgrowth_Flag = 0;
    default_Params.Sub_Sample_Num = 1;
    default_Params.F_classifier{1} = [];
    default_Params.F_classifier{2} = [];
    default_Params.Classifier_Type_ind = [1 1];
    default_Params.training_sample_number = [30 30];
    default_Params.LengthThreshold = [4 4];
    default_Params.CurvatureThreshold = [0.1 0.1];
    default_Params.CoefAlpha = [1.8 1.7];
    default_Params.IternationNumber = [2 2];
    default_Params.CannyHigherThreshold = [80 80];
    default_Params.CannyLowerThreshold = [80 80];
    default_Params.nofiguredisruption = 1;
    default_Params.savestepfigures = 0;
    default_Params.showdetailmessages = 0;
    default_Params.channel_specific = 0;
    
    if(indexFilamentSegmentationProcess==0)
        this_MD.addProcess(FilamentSegmentationProcess(this_MD,'funParams',default_Params));
    end
    
    % if user give the whole movie stat, use this
    if(~ismember('whole_movie_filename',ip.UsingDefaults))
        this_MD = filament_segmentation_continue(this_MD,default_Params,whole_movie_filename);
    else
        this_MD = filament_segmentation_continue(this_MD,default_Params);
    end
    
end
