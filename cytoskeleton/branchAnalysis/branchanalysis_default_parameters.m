function default_parameter_set = branchanalysis_default_parameters(this_MD)
% Set default parameters for branch analysis filament segmentation package, always for all the 5 processes, in order

% intialize
default_parameter_set=cell(1,5);


%%  % 1 threshold


default_Params=[];

default_Params.ChannelIndex = [1 2];
default_Params.GaussFilterSigma=2;
default_Params.ThresholdValue=[];
default_Params.OutputDirectory = [this_MD.outputDirectory_,filesep,'SegmentationPackage',filesep,'thres'];
default_Params.MethodIndx=4;
default_Params.ProcessIndex=[];
default_Params.MaxJump=10;
default_Params.BatchMode=0;

default_parameter_set{1} = default_Params;

%%  % 2 mask refine

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

default_parameter_set{2} = default_Params;

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


default_parameter_set{3} = default_Params;


%% % 4 steerable filter
default_Params=[];

default_Params.ChannelIndex= 2;
default_Params.BaseSteerableFilterSigma= 1;
default_Params.Levelsofsteerablefilters= 2;
default_Params.Sub_Sample_Num= 1;
default_Params.ImageFlattenFlag= 2;

default_parameter_set{4} = default_Params;

%% % 5 filament segmentation
default_Params=[];

default_Params.ChannelIndex = 2;
default_Params.StPace_Size = [3 3];
default_Params.StPatch_Size = [21 21];
default_Params.st_lowerbound_localthresholding = [90 90] ;
default_Params.IntPace_Size = [3 3];
default_Params.IntPatch_Size = [21 21 ];
default_Params.int_lowerbound_localthresholding = [90 90 ];
default_Params.Combine_Way{1} = 'geo_based_GM';
default_Params.Combine_Way{2} = 'geo_based_GM';
% mask as 6, meaning only running filament segmentation for marked cell area
default_Params.Cell_Mask_ind = [6 6];
default_Params.Whole_movie_ind = [1 1];
default_Params.Whole_movie_stat_cell{1} = [];
default_Params.Whole_movie_stat_cell{2} = [];
default_Params.Rerun_WholeMovie = 1;
default_Params.VIF_Outgrowth_Flag = 0;
default_Params.Sub_Sample_Num = 1;
default_Params.F_classifier{1} = [];
default_Params.F_classifier{2} = [];
default_Params.Classifier_Type_ind = [1 1];
default_Params.training_sample_number = [30 30];
default_Params.LengthThreshold = [4 4];
default_Params.CurvatureThreshold = [0.1 0.1];
default_Params.CoefAlpha = [1.3 1.3];
default_Params.IternationNumber = [2 2];
default_Params.CannyHigherThreshold = [80 80];
default_Params.CannyLowerThreshold = [80 80];
default_Params.nofiguredisruption = 1;
default_Params.savestepfigures = 0;
default_Params.showdetailmessages = 0;
default_Params.channel_specific = 0;

default_parameter_set{5} = default_Params;

