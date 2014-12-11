function load_MD_run_filament_analysis_package(this_MD, Parameter_MD, varargin)
% Function of single image filament segmentation with input MD from other
%               successfully segmented movie for the parameters in MD

% Input:      this_MD:              the movieData for data to be segmented
%             Parameter_MD:    a loaded MD with good segmentation parameters
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

%% get the image dir and make a MD just for this image
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
            
            this_MD.addProcess(ThresholdProcess(this_MD,'funParams',given_Params));
            this_MD = thresholdMovie(this_MD,this_MD.processes_{iPro}.funParams_);
            
        else
            % 2 mask refine
            if(strcmp(Parameter_MD.processes_{iPro}.getName, 'Mask Refinement'))
                given_Params = Parameter_MD.processes_{iPro}.funParams_;
                given_Params.OutputDirectory = [this_MD.outputDirectory_,filesep,'FilamentAnalysisPackage',filesep,'maskrefine'];
                
                this_MD.addProcess(MaskRefinementProcess(this_MD,'funParams',given_Params));
                this_MD = refineMovieMasks(this_MD,this_MD.processes_{iPro}.funParams_);
                
                
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
                        this_MD.addProcess(SteerableFilteringProcess(this_MD,'funParams',given_Params));
                        this_MD = steerable_filter_forprocess(this_MD,given_Params);
                        
                    else
                        % 5 filament segmentation
                        if(strcmp(Parameter_MD.processes_{iPro}.getName, 'Filament Segmentation'))
                            given_Params = Parameter_MD.processes_{iPro}.funParams_;
                            this_MD.addProcess(FilamentSegmentationProcess(this_MD,'funParams',given_Params));
                            
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
    
    if(nProcess==0)
        %%        
        % when there is no input MD, and no existing processes
        % use the default for all of the processes
        
        % 1 threshold
        
%         this_MD.addProcess(ThresholdProcess(this_MD));
         this_MD.addProcess(ThresholdProcess(this_MD,'FilamentAnalysisPackage'));
        this_MD = thresholdMovie(this_MD);
        
        % 2 mask refine
        % with the default of mask refinement
         this_MD.addProcess(MaskRefinementProcess(this_MD,'FilamentAnalysisPackage'));
%         this_MD.addProcess(MaskRefinementProcess(this_MD));
        this_MD = refineMovieMasks(this_MD);
        
        % 3 image flatten
%         this_MD.addProcess(ImageFlattenProcess(this_MD,['FilamentAnalysisPackage', filesep,'ImageFlatten']));
        this_MD.addProcess(ImageFlattenProcess(this_MD));
        this_MD = image_flatten(this_MD);
        
        % 4 steerable filter
%         this_MD.addProcess(SteerableFilteringProcess(this_MD,['FilamentAnalysisPackage', filesep,'SteerableFiltering']));
        this_MD.addProcess(SteerableFilteringProcess(this_MD));
        this_MD = steerable_filter_forprocess(this_MD);
        
        % 5 filament segmentation
%         this_MD.addProcess(FilamentSegmentationProcess(this_MD,['FilamentAnalysisPackage', filesep,'FilamentSegmentation']));
        this_MD.addProcess(FilamentSegmentationProcess(this_MD));
        
        default_Params = filament_segmentation_process_default_param(this_MD);
        default_Params.Cell_Mask_ind = 5;
        default_Params.CoefAlpha = 1.8;
        default_Params.ChannelIndex = 1:2;
        
        this_MD = filament_segmentation(this_MD, default_Params);
        
    else
        %% or, just run it with its existing processes;
                
        for iPro =  1 : numel(this_MD.processes_)
            % 1 threshold
            
            if(strcmp(this_MD.processes_{iPro}.getName, 'Thresholding'))
                this_MD = thresholdMovie(this_MD);
            else
                % 2 mask refine
                if(strcmp(this_MD.processes_{iPro}.getName, 'Mask Refinement'))
                    this_MD = refineMovieMasks(this_MD);                    
                    
                else
                    % 3 image flatten
                    if(strcmp(this_MD.processes_{iPro}.getName, 'Image Flatten'))
                        this_MD = image_flatten(this_MD);
                        
                    else
                        % 4 steerable filter
                        if(strcmp(this_MD.processes_{iPro}.getName, 'Steerable filtering'))
                            this_MD = steerable_filter_forprocess(this_MD);
                            
                        else
                            % 5 filament segmentation
                            if(strcmp(this_MD.processes_{iPro}.getName, 'Filament Segmentation'))
                                
                                default_Params = filament_segmentation_process_default_param(this_MD);
                                default_Params.Cell_Mask_ind = 5;
                                default_Params.CoefAlpha = 1.8;
                                default_Params.ChannelIndex = 1:2;
   
                                % if user give the whole movie stat, use this
                                if(~ismember('whole_movie_filename',ip.UsingDefaults))
                                    this_MD = filament_segmentation(this_MD,default_Params,whole_movie_filename);
                                else
                                    this_MD = filament_segmentation(this_MD,default_Params);
                                end
                            end
                        end
                    end
                end
            end
        end
        %%
        
    end
    
end



