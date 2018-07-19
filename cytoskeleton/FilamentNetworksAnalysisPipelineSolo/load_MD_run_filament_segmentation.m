function this_MD = load_MD_run_filament_segmentation(this_MD, Parameter_MD, varargin)
% Function of filament segmentation with input this_MD from other
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

ip.addOptional('run_with_new_param',  0, @isnumeric);
ip.addOptional('input_parameter_set', [], @(x) iscell(x) | isempty(x));
ip.addOptional('whole_movie_filename', [], @(x) ischar(x) | isempty(x));

ip.parse(this_MD,Parameter_MD,varargin{:});
whole_movie_filename= ip.Results.whole_movie_filename;
run_with_new_param = ip.Results.run_with_new_param;
input_parameter_set = ip.Results.input_parameter_set;

% set from the function for default, which is the function to update if
% later changed default
default_parameter_cell= default_parameter_filament_segmentation(this_MD);

set_parameter_cell = default_parameter_cell;

if(~isempty(input_parameter_set))
    for i = 1 : numel(input_parameter_set)
        if(~isempty(input_parameter_set{i}))
%              set_parameter_cell{i} = parameter_replace_n_merge(set_parameter_cell{i}, input_parameter_set{i});
            set_parameter_cell{i} = input_parameter_set{i};
        end
    end    
end
    
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
nPackage = numel(this_MD.packages_);

% if there is no filament analysis package,
% add filament analysis package
if(indexFilamentPackage == 0)
    this_MD.addPackage(FilamentAnalysisPackage(this_MD));
    
    
    for i = 1 : numel(this_MD.packages_)
        if(strcmp(this_MD.packages_{i}.getName,'FilamentAnalysis')==1)
            indexFilamentPackage = i;
            break;
        end
    end
        
    mkdir([this_MD.outputDirectory_,filesep,'FilamentAnalysisPackage']);
end


 %     % check if there is each of the process
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



nProcess = numel(this_MD.processes_);
nPackage = numel(this_MD.packages_);



% start process by process

if(~isempty(Parameter_MD))
    %% % if there is input Parameter_MD
    
    
    for iPro =  1 : numel(Parameter_MD.processes_)
        % 1 threshold
        
        if(strcmp(Parameter_MD.processes_{iPro}.getName, 'Thresholding'))
            given_Params = Parameter_MD.processes_{iPro}.funParams_;
            given_Params.OutputDirectory = [this_MD.outputDirectory_,filesep,'FilamentAnalysisPackage',filesep,'thres'];
            
            
            if(indexCellSegProcess==0)
                this_MD.addProcess(ThresholdProcess(this_MD,'default_Params',given_Params));
                this_MD = thresholdMovie(this_MD,this_MD.processes_{iPro}.funParams_);
            else
                this_MD = thresholdMovie(this_MD,given_Params);
            end          
            
            
        else
            % 2 mask refine
            if(strcmp(Parameter_MD.processes_{iPro}.getName, 'Mask Refinement'))
                given_Params = Parameter_MD.processes_{iPro}.funParams_;
                given_Params.OutputDirectory = [this_MD.outputDirectory_,filesep,'FilamentAnalysisPackage',filesep,'maskrefine'];
                
                if(indexCellRefineProcess==0)                    
                    this_MD.addProcess(MaskRefinementProcess(this_MD,'default_Params',given_Params));                    
                    this_MD = refineMovieMasks(this_MD,this_MD.processes_{iPro}.funParams_);
                else
                    this_MD = refineMovieMasks(this_MD,given_Params);
                end
                
            else
                % 3 image flatten
                if(strcmp(Parameter_MD.processes_{iPro}.getName, 'Image Flatten'))
                    given_Params = Parameter_MD.processes_{iPro}.funParams_;
                    
                    if(indexFlattenProcess==0)
                        this_MD.addProcess(ImageFlattenProcess(this_MD));
                    end
                    
                    this_MD = image_flatten(this_MD,given_Params);
                    
                else
                    % 4 steerable filter
                    if(strcmp(Parameter_MD.processes_{iPro}.getName, 'Steerable filtering'))
                        given_Params = Parameter_MD.processes_{iPro}.funParams_;
                        if(indexSteerabeleProcess==0)                            
                            this_MD.addProcess(SteerableFilteringProcess(this_MD,'default_Params',given_Params));
                        end
                        this_MD = steerable_filter_forprocess(this_MD,given_Params);
                        
                    else
                        % 5 filament segmentation
                        if(strcmp(Parameter_MD.processes_{iPro}.getName, 'Filament Segmentation'))
                            given_Params = Parameter_MD.processes_{iPro}.funParams_;
                            
                            if(indexFilamentSegmentationProcess==0)
                            this_MD.addProcess(FilamentSegmentationProcess(this_MD,'default_Params',given_Params));
                            end
                            
                            % if user give the whole movie stat, use this
                            if(~ismember('whole_movie_filename',ip.UsingDefaults))
                                this_MD = filament_segmentation(this_MD,given_Params,whole_movie_filename);
                            else
                                this_MD = filament_segmentation(this_MD,given_Params);
                            end
                            %save MD to disk
                            this_MD.save();
                        end
                        %save MD to disk
                        this_MD.save();
                    end
                    %save MD to disk
                    this_MD.save();
                end
            end
        end
    end
    
else % else for if there is a input parameter MD
   
    
    %%  % 1 threshold
    
    %%   % when there is no input this_MD, and no existing processes
    % use the default for all of the processes
     
    default_Params = set_parameter_cell{1};
           
    if(indexCellSegProcess==0)
        this_MD.addProcess(ThresholdProcess(this_MD,'FilamentAnalysisPackage',default_Params));
        
        for i = 1 : numel(this_MD.processes_)
            if(strcmp(this_MD.processes_{i}.getName,'Thresholding')==1)
                indexCellSegProcess = i;
                break;
            end
        end
        
    end
    
    % when default parameter say no channel to run, skip this step    
    if(~isempty(default_Params.ChannelIndex))
     this_MD = thresholdMovie(this_MD,default_Params);
    end
    
    %%  % 2 mask refine
    % with the default of mask refinement
    default_Params = set_parameter_cell{2};
    
    if(indexCellRefineProcess==0)
        this_MD.addProcess(MaskRefinementProcess(this_MD,'FilamentAnalysisPackage',default_Params));
        
        for i = 1 : numel(this_MD.processes_)
            if(strcmp(this_MD.processes_{i}.getName,'Mask Refinement')==1)
                indexCellRefineProcess = i;
                break;
            end
        end
    end
    
    % when default parameter say no channel to run, skip this step    
    if(~isempty(default_Params.ChannelIndex))    
     this_MD = refineMovieMasks(this_MD,default_Params);
    end
    
    %%   % 3 image flatten
     default_Params = set_parameter_cell{3};
   
    if(indexFlattenProcess==0)
        this_MD.addProcess(ImageFlattenProcess(this_MD,'funParams',default_Params));
        
        for i = 1 : numel(this_MD.processes_)
            if(strcmp(this_MD.processes_{i}.getName,'Image Flatten')==1)
                indexFlattenProcess = i;
                break;
            end
        end
    end
    
    if(run_with_new_param ==1)
        this_MD = image_flatten(this_MD,default_Params);
    else
        this_MD = image_flatten_norepeating(this_MD,default_Params);
    end
    %save MD to disk
    this_MD.save();
    
    %% % 4 steerable filter
    
    default_Params = set_parameter_cell{4};
       
    if(indexSteerabeleProcess==0)
         this_MD.addProcess(SteerableFilteringProcess(this_MD,'funParams',default_Params));
         
         for i = 1 : numel(this_MD.processes_)
             if(strcmp(this_MD.processes_{i}.getName,'Steerable filtering')==1)
                 indexSteerabeleProcess = i;
                 break;
             end
         end
         
    end
    
    if(run_with_new_param ==1)
        this_MD = steerable_filter_forprocess(this_MD,default_Params);
    else
        this_MD = steerable_filter_forprocess_continue(this_MD,default_Params);
    end
    
    %save MD to disk
    this_MD.save();
        
    %%
    % 5 filament segmentation
    default_Params = set_parameter_cell{5};
      
    if(indexFilamentSegmentationProcess==0)
        this_MD.addProcess(FilamentSegmentationProcess(this_MD,'funParams',default_Params));
        
        for i = 1 : numel(this_MD.processes_)
            if(strcmp(this_MD.processes_{i}.getName,'Filament Segmentation')==1)
                indexFilamentSegmentationProcess = i;
                break;
            end
        end
    end
    
    
    %% run the filament segmentation with/out replacement, with/out whole movie filename
      
    % if user give the whole movie stat, use this
    if(~ismember('whole_movie_filename',ip.UsingDefaults))
        % run_with_new_param  is 1, means the user want replace the
        % filament segmentation results
        if(run_with_new_param ==1)
            this_MD = filament_segmentation(this_MD,default_Params,whole_movie_filename);
        else
            this_MD = filament_segmentation_continue(this_MD,default_Params,whole_movie_filename);
        end
    else        
        if(run_with_new_param ==1)
            this_MD = filament_segmentation(this_MD,default_Params);
        else
            this_MD = filament_segmentation_continue(this_MD,default_Params);
        end
    end    
    
    %save MD to disk
    this_MD.save();
end


%Finally, save MD to disk
this_MD.save();