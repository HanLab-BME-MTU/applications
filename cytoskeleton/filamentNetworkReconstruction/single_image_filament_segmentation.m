function [current_seg,current_seg_orientation, current_model,RGB_seg_orient_heat_map,MAX_st_res,nms] ...
    = single_image_filament_segmentation(filename, Parameter_MD, varargin)
% Function of single image filament segmentation with input MD from other
%               successfully segmented movie for the parameters in MD

% Input:      filename:        the filename of the image to be segmented
%             Parameter_MD:    a loaded MD with good segmentation parameters
%                              if none, put [], so a default setting will be used
%                              with (1) Otsu with smoothing 1 (2) mask
%                              refine with 1 object, (3) image flatten with
%                              square, (4) steerable filter with [1 2], (5)
%                              segmentation with geo based alpha=2
%             pick_channel:    optional input, which channel to use in case there are more
%                               than one channel in the MD, default 1
%             keep_steps:      optional input, only ofr debugging, if 1,
%                               the steps are kept, default 0
%             output_dir:      optional input, if given, the results will
%                                be saved in this dir; if not given, will be in the image
%                                folder
%             whole_movie_filename: optional input, if given, the whole
%                                movie statistics will be loaded from this filename, which
%                                will be used in the filament segmentsation; by default,
%                                nothing is input, so whole movie stat will be calculated,
%                                well, in this case, without any meaningful effect

% Output:     current_model: the filament-by-filament model
%             current_seg: the black-white segmentation results
%             current_seg_orientation: the orienation of segmented pixels
%             RGB_seg_orient_heat_map: heat map for display
%             MAX_st_res: is the steerable filtering result
%             nms: the non-maximum-surpress ST
%             These variables will be save in the result .mat file in the
%             name of the image name plus _filament_seg_results.mat.
%             The bw segmentation and the heatmap is saved as images.

% Created 05 2014 by Liya Ding, Matlab R2012b

ip = inputParser;
ip.addRequired('filename',@ischar);
ip.addRequired('Parameter_MD',@(x) isa(x,'MovieData') | isempty(x));
ip.addOptional('pick_channel', 1,@isnumeric);
ip.addOptional('keep_steps', 0,@isnumeric);
ip.addOptional('output_dir', ' ', @ischar);
ip.addOptional('whole_movie_filename', ' ', @ischar);

ip.parse(filename,Parameter_MD,varargin{:});
pick_channel = ip.Results.pick_channel;
keep_steps = ip.Results.keep_steps;
output_dir = ip.Results.output_dir;
whole_movie_filename= ip.Results.whole_movie_filename;

%% get the image dir and make a MD just for this image
% first parse the input image path and output path
filename = GetFullPath(filename);
index_1 = find(filename=='/' | filename=='\');
ROOT_DIR = filename(1:max(index_1));
ROOT_DIR(index_1)=filesep;

file_image_full_name = filename(max(index_1)+1:end);
index_2 = find(file_image_full_name=='.');
file_image_only_name = file_image_full_name(1:max(index_2)-1);

% if there is user input for output_dir, use that as main working dir
if(~ismember('output_dir',ip.UsingDefaults))
    output_dir = GetFullPath(output_dir);
    ROOT_DIR = output_dir;
    if (~exist(ROOT_DIR,'dir'))
        mkdir(ROOT_DIR);
    end
end

% this is the major working dir, where all the steps will be
% by default this folder will be deleted, unless keep_steps==1
working_ROOT_DIR = [ROOT_DIR, filesep,'working',filesep];
if (~exist(working_ROOT_DIR,'dir'))
    mkdir(working_ROOT_DIR);
end


% start to get channels
channels_obj_one = [];
channels_obj_all=[];

% if there is input MD, copy the number of channels, other wise just 1
% channel
if(~isempty(Parameter_MD))
    NChannel = numel(Parameter_MD.channels_);
else
    NChannel=1;
    pick_channel=1;
end

% duplicate channels if the MD has more than 1
for iC = 1 : NChannel
    
    % copy the image into a folder for movieData structure
    imagefolder_DIR = [working_ROOT_DIR,'single_run_image',num2str(iC)];
    if (~exist(imagefolder_DIR,'dir'))
        mkdir(imagefolder_DIR);
    end
    
    % make copy of the image to channel folders
    copyfile(filename,[imagefolder_DIR,filesep,'temp',num2str(iC),'.tif']);
    
    % use this folder for a channel
    channels_obj_one = Channel([working_ROOT_DIR,'single_run_image',num2str(iC)]);
    channels_obj_one.getImageFileNames();
    channels_obj_one.sanityCheck();
    
    % stack the channels
    channels_obj_all = [channels_obj_all channels_obj_one];
end

% build the MD
this_MD = MovieData(channels_obj_all,working_ROOT_DIR);
this_MD.setPath(working_ROOT_DIR);
this_MD.setFilename('movieData.mat')
this_MD.sanityCheck();

this_MD.addPackage(FilamentAnalysisPackage);

% start process by process
% if there is input Parameter_MD
if(~isempty(Parameter_MD))
    for iPro =  1 : numel(Parameter_MD.processes_)
        % 1 threshold
        
        if(strcmp(Parameter_MD.processes_{iPro}.getName, 'Thresholding'))
            given_Params = Parameter_MD.processes_{iPro}.funParams_;
            given_Params.OutputDirectory = [this_MD.outputDirectory_,filesep,'thres'];
            
            given_Params.ChannelIndex = pick_channel;
            this_MD.addProcess(ThresholdProcess(this_MD,'funParams',given_Params));
            this_MD = thresholdMovie(this_MD,this_MD.processes_{iPro}.funParams_);
            
            
        else
            % 2 mask refine
            if(strcmp(Parameter_MD.processes_{iPro}.getName, 'Mask Refinement'))
                given_Params = Parameter_MD.processes_{iPro}.funParams_;
                given_Params.OutputDirectory = [this_MD.outputDirectory_,filesep,'maskrefine'];
                given_Params.ChannelIndex = pick_channel;
                this_MD.addProcess(MaskRefinementProcess(this_MD,'funParams',given_Params));
                this_MD = refineMovieMasks(this_MD,this_MD.processes_{iPro}.funParams_);
                
                
            else
                % 3 image flatten
                if(strcmp(Parameter_MD.processes_{iPro}.getName, 'Image Flatten'))
                    given_Params = Parameter_MD.processes_{iPro}.funParams_;
                    given_Params.ChannelIndex = pick_channel;
                    this_MD.addProcess(ImageFlattenProcess(this_MD));
                    this_MD = image_flatten(this_MD,given_Params);
                    
                    
                    
                else
                    % 4 steerable filter
                    if(strcmp(Parameter_MD.processes_{iPro}.getName, 'Steerable filtering'))
                        given_Params = Parameter_MD.processes_{iPro}.funParams_;
                        given_Params.ChannelIndex = pick_channel;
                        this_MD.addProcess(SteerableFilteringProcess(this_MD,'funParams',given_Params));
                        this_MD = steerable_filter_forprocess(this_MD,given_Params);
                        
                        STOutputDir = [this_MD.processes_{iPro}.outFilePaths_{pick_channel}];
                        
                        load([STOutputDir,'/steerable_temp',num2str(pick_channel),'.tif.mat'],...
                            'MAX_st_res','nms');
                        
                        
                    else
                        % 5 filament segmentation
                        if(strcmp(Parameter_MD.processes_{iPro}.getName, 'Filament Segmentation'))
                            given_Params = Parameter_MD.processes_{iPro}.funParams_;
                            given_Params.ChannelIndex = pick_channel;
                            this_MD.addProcess(FilamentSegmentationProcess(this_MD,'funParams',given_Params));
                            
                            if(given_Params.Cell_Mask_ind(pick_channel)==4)
                                given_Params.Cell_Mask_ind(pick_channel)=1;
                            end
                            
                            % if user give the whole movie stat, use this
                            if(~ismember('whole_movie_filename',ip.UsingDefaults))
                               this_MD = filament_segmentation(this_MD,given_Params,whole_movie_filename);
                            else
                                this_MD = filament_segmentation(this_MD,given_Params);
                            end
                            
                            DataOutputDir = [this_MD.processes_{iPro}.outFilePaths_{pick_channel},'/DataOutput'];
                            
                            load([DataOutputDir,'/steerable_vote_temp',num2str(pick_channel),'.tif.mat'],...
                                'RGB_seg_orient_heat_map', ...
                                'current_seg', ...
                                'current_model', 'current_seg_orientation');                            
                            
                        end
                    end
                end
            end
        end
    end
else
    
    % when there is no input MD use the default for all of the processes
    % 1 threshold
    dafault_Params=[];
    dafault_Params.ChannelIndex=1;
    dafault_Params.GaussFilterSigma=2;
    dafault_Params.ThresholdValue=[];
    dafault_Params.OutputDirectory = [this_MD.outputDirectory_,filesep,'thres'];
    dafault_Params.MethodIndx=3;
    dafault_Params.ProcessIndex=[];
    dafault_Params.MaxJump=10;
    dafault_Params.BatchMode=0;    
    
    this_MD.addProcess(ThresholdProcess(this_MD,'funParams',dafault_Params));
    this_MD = thresholdMovie(this_MD,dafault_Params);
    % 2 mask refine
    % with the default of mask refinement
    this_MD.addProcess(MaskRefinementProcess(this_MD));
    this_MD = refineMovieMasks(this_MD);
    
    % 3 image flatten
    this_MD.addProcess(ImageFlattenProcess(this_MD));
    this_MD = image_flatten(this_MD);
    % 4 steerable filter
    this_MD.addProcess(SteerableFilteringProcess(this_MD));
    this_MD = steerable_filter_forprocess(this_MD);
    
    STOutputDir = [this_MD.processes_{4}.outFilePaths_{pick_channel}];
    load([STOutputDir,'/steerable_temp',num2str(pick_channel),'.tif.mat'],...
        'MAX_st_res','nms');
    
    % 5 filament segmentation
    
    this_MD.addProcess(FilamentSegmentationProcess(this_MD));
    default_Params = this_MD.processes_{5}.funParams_;
    default_Params.Cell_Mask_ind(pick_channel)=1;
    default_Params.CoefAlpha = 1.8;
            
    this_MD = filament_segmentation(this_MD, default_Params);
    
    DataOutputDir = [this_MD.processes_{5}.outFilePaths_{pick_channel},'/DataOutput'];
    
    load([DataOutputDir,'/steerable_vote_temp',num2str(pick_channel),'.tif.mat'],...
        'RGB_seg_orient_heat_map', ...
        'current_seg', ...
        'current_model', 'current_seg_orientation');
    
end

% save the output
save([ROOT_DIR, filesep, file_image_only_name,'_filament_seg_results.mat'],...
    'RGB_seg_orient_heat_map', ...
    'MAX_st_res','nms', 'current_seg', ...
    'current_model', 'current_seg_orientation');

% save two images, bw segmentation and heatmap with orientation
imwrite(current_seg,[ROOT_DIR, filesep, file_image_only_name,'_seg.tif']);
imwrite(RGB_seg_orient_heat_map,[ROOT_DIR, filesep, file_image_only_name,'_heat.tif']);

if(keep_steps==0)
    % delete everything in the working folder
    rmdir(working_ROOT_DIR,'s');
end

