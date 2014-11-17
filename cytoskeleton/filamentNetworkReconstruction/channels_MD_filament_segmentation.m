function this_MD = channels_MD_filament_segmentation(image_folder_names, Parameter_MD, varargin)
% Function of image filament segmentation with input MD from other
%               successfully segmented movie for the parameters in MD
%               this is run in the exact same way as in GUI strcture

% Input:      image_folder:    the folder of images to be segmented(one
%                              channel by one channel
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

% Output:     this_MD: the movieData object, everything is saved as
%                           in GUI run

% Created 05 2014 by Liya Ding, Matlab R2012b

ip = inputParser;
ip.addRequired('image_folder_names',@(x) ischar(x) | iscell(x));
ip.addRequired('Parameter_MD',@(x) isa(x,'MovieData') | isempty(x));
ip.addOptional('pick_channel', 0,@isnumeric);
ip.addOptional('keep_steps', 0,@islogic);
ip.addOptional('output_dir', ' ', @ischar);
ip.addOptional('whole_movie_filename', ' ', @ischar);

ip.parse(image_folder_names,Parameter_MD,varargin{:});
pick_channel = ip.Results.pick_channel;
keep_steps = ip.Results.keep_steps;
output_dir = ip.Results.output_dir;
whole_movie_filename= ip.Results.whole_movie_filename;



%% get the image dir and make a MD just for this image
% first parse the input image path and output path
if(~iscell(image_folder_names))
    % if this is not a cell(so it is a string), put it into a cell
    % structure
    temp = image_folder_names;
    image_folder_names = cell(1,1);
    image_folder_names{1}=temp;
end

% if user didn't input what channel to work on, use all the channels
% well, this sounds pretty obvious
if(pick_channel==0)
    pick_channel = 1:numel(image_folder_names);
end


for iC = 1 : numel(image_folder_names)
    image_folder_names{iC} = GetFullPath(image_folder_names{iC});
end

filename = image_folder_names{1};
index_1 = find(filename=='/' | filename=='\');
filename(index_1)=filesep;

if(filename(end)==filesep)
    filename = filename(1:end-1);
end

index_1 = find(filename==filesep);
ROOT_DIR = filename(1:max(index_1));

% if there is user input for output_dir, use that as main working dir
if(min(output_dir==' ')==0)
    output_dir = GetFullPath(output_dir);
    ROOT_DIR = output_dir;
    if (~exist(ROOT_DIR,'dir'))
        mkdir(ROOT_DIR);
    end
end

% this is the major working dir, where all the steps will be
% by default this folder will be deleted, unless keep_steps==1
working_ROOT_DIR = [ROOT_DIR, filesep,'FilamentSegmentationResults',filesep];
if (~exist(working_ROOT_DIR,'dir'))
    mkdir(working_ROOT_DIR);
end

% define the channel number
NChannel = numel(image_folder_names);

% if there is input parameter MD, then match the number of channels by
% repeating the last channel
if(~isempty(Parameter_MD))
    PM_NChannel = numel(Parameter_MD.channels_);
    if(PM_NChannel>NChannel)
        for iC = NChannel+1:PM_NChannel
            image_folder_names{iC} = image_folder_names{NChannel};
        end
    end
end

% and then get the number of channel again
NChannel = numel(image_folder_names);

% start to get channels
channels_obj_one = [];
channels_obj_all=[];

for iC = 1 : NChannel
    
    % use each folder for a channel
    channels_obj_one = Channel(image_folder_names{iC});
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
                                                
                    else
                        % 5 filament segmentation
                        if(strcmp(Parameter_MD.processes_{iPro}.getName, 'Filament Segmentation'))
                            given_Params = Parameter_MD.processes_{iPro}.funParams_;
                            given_Params.ChannelIndex = pick_channel;
                            this_MD.addProcess(FilamentSegmentationProcess(this_MD,'funParams',given_Params));
                            
                            % if there is only one channel, use this
                            % channel segmentation
                            if (numel(pick_channel)==1)
                                if(given_Params.Cell_Mask_ind(pick_channel)==4)
                                    given_Params.Cell_Mask_ind(pick_channel)=1;
                                end
                            end
                            
                            % if user give the whole movie stat, use this
                            if( min(whole_movie_filename==' ')==0 )
                                this_MD = filament_segmentation(this_MD,given_Params,whole_movie_filename);
                            else
                                this_MD = filament_segmentation(this_MD,given_Params);
                            end
                            
                            
                                for iC = pick_channel(:)
                                    brief_results = [working_ROOT_DIR,filesep,'Channel_',num2str(iC),'_HeatResults/'];
                                    if(exist(brief_results,'dir')==0)
                                        mkdir(brief_results);
                                    end
                                    HeatOutputDir = [this_MD.processes_{iPro}.outFilePaths_{iC},'/HeatOutput/Enh'];
                                    
                                    copyfile([HeatOutputDir,'/segment_heat_*.tif'], brief_results);
                                end
                           
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
    dafault_Params.ChannelIndex=1:NChannel;
    dafault_Params.GaussFilterSigma=1;
    dafault_Params.ThresholdValue=[];
    dafault_Params.OutputDirectory = [this_MD.outputDirectory_,filesep,'thres'];
    dafault_Params.MethodIndx=3;
    dafault_Params.ProcessIndex=[];
    dafault_Params.MaxJump=0;
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
    
    % 5 filament segmentation
    
    this_MD.addProcess(FilamentSegmentationProcess(this_MD));
    default_Params = this_MD.processes_{5}.funParams_;
    default_Params.Cell_Mask_ind(pick_channel)=1;
    default_Params.CoefAlpha = 1.8;
            
    this_MD = filament_segmentation(this_MD, default_Params);
    
    for iC = pick_channel(:)
        brief_results = [working_ROOT_DIR,filesep,'Channel_',num2str(iC),'_HeatResults/'];
        if(exist(brief_results,'dir')==0)
            mkdir(brief_results);
        end
        HeatOutputDir = [this_MD.processes_{5}.outFilePaths_{iC},'/HeatOutput/Enh'];
        
        copyfile([HeatOutputDir,'/segment_heat_*.tif'], brief_results);
    end
end

