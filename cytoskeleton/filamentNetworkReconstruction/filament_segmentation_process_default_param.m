function default_Params = filament_segmentation_process_default_param(MD)
% Set default parameters
default_Params=[];
default_Params.StPace_Size = 3;
default_Params.StPatch_Size  = 21;
default_Params.st_lowerbound_localthresholding  = 90; % default 90%
default_Params.IntPace_Size = 3;
default_Params.IntPatch_Size  = 21;
default_Params.int_lowerbound_localthresholding  = 90; % default 90%
default_Params.Combine_Way = 'geo_based_GM';
default_Params.Whole_movie_ind = 2;
default_Params.Rerun_WholeMovie = 0;
default_Params.VIF_Outgrowth_Flag = 0;
default_Params.Sub_Sample_Num = 1;
default_Params.Classifier_Type_ind=1;
default_Params.training_sample_number=30;
default_Params.LengthThreshold=4;
default_Params.CurvatureThreshold=0.1;
default_Params.IternationNumber=2;
default_Params.CannyHigherThreshold=80;
default_Params.CannyLowerThreshold=80;
default_Params.nofiguredisruption = 1;
default_Params.savestepfigures = 0;
default_Params.showdetailmessages = 0;
default_Params.channel_specific=0;

default_Params.CoefAlpha = 2;
default_Params.Cell_Mask_ind = 1;
default_Params.ChannelIndex = 1:numel(MD.channels_);

default_Params.F_classifier = cell(1,max(default_Params.ChannelIndex));
default_Params.Whole_movie_stat_cell = cell(1,max(default_Params.ChannelIndex));

