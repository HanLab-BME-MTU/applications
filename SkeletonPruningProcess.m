classdef SkeletonPruningProcess < ImageAnalysisProcess
%A process class for pruning mask skeletons in every frame of a 3D movie
%using prune3DMovieSkeletonBranches.m
%
% Because this function has multiple inputs per directory, the
% inputFilePaths field will be a 3xM cell array, where M is the number of
% channels. The rows contain the input as follows:
%   Row 1 - The input mask directories
%   Row 2 - The input skeleton graph directories
%   Row 3 - The input mask geometry directories
%
% Hunter Elliott
% 4/2011
    properties (SetAccess = protected, GetAccess = public)
        
        
    end
    
    methods (Access = public)
    
     % ------ Constructor -----%
        function obj = SkeletonPruningProcess(owner,funParams)
                                              
            if nargin == 0
                super_args = {};
            else                
                
                super_args{1} = owner;
                super_args{2} = 'SkeletonPruning';
                super_args{3} = @prune3DMovieSkeletonBranches;                               
                
                if nargin < 2 || isempty(funParams)                                       
                    
                    %----Defaults----%      
                                        
                    funParams.SkelProcessIndex = [];%No default 
                    funParams.PruneParam = [];%For storing pruneSkeletonGraph.m-specific parameters
                    funParams.OutputDirectory = ...
                        [owner.outputDirectory_  filesep 'pruned_skeleton_graphs'];                    
                    funParams.BatchMode = false;                                                      
                                    
                    
                end
                
                super_args{4} = funParams;    
                                
            end
            
            obj = obj@ImageAnalysisProcess(super_args{:});
        end
        
        % ------ Set / Get Methods ----- %
        
        %Uses generic ImageAnalysisProcess methods
        
        % ----- Output Check/Load Methods ----- %
                 
                              
        function skelGraph = loadChannelOutput(obj,iChan,iFrame)
            
            if nargin < 3 || isempty(iChan) || isempty(iFrame)
                error('You must specify a channel number and frame number!')                
            elseif ~obj.checkChanNum(iChan)
                error('Invalid channel number!')
            elseif iFrame < 1 || iFrame > obj.owner_.nFrames_
                error('Invalid frame number!')
            end
            
            if ~obj.checkChannelOutput(iChan);
                error('The specified channel does not have valid skeleton graph files!')
            end
                        
            skelNames = dir([obj.outFilePaths_{iChan} filesep '*.mat']);
            
            skelGraph = load([obj.outFilePaths_{iChan} filesep skelNames(iFrame).name]);
            
            fNames = fieldnames(skelGraph);
            
            if numel(fNames) ~= 3 || ~isfield(skelGraph,'vertices') || ...
                    ~isfield(skelGraph,'edges') || ...
                    ~isfield(skelGraph,'edgePaths')
                error('Invalid skeleton graph file!')                                            
            end
                
            
        end


    end
    
end