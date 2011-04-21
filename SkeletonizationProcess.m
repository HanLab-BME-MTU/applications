classdef SkeletonizationProcess < ImageProcessingProcess
% A process class for performing skeletonization on a movie
%     
% This clas inherits ImageProcessingProcess, but is actually a hybrid image
% processing/analysis function, because the output skeletons are images,
% but the skeleton graphs are not.
%
% Hunter Elliott
% 12/2010    
%
    properties (SetAccess = protected, GetAccess = public)
        
        
    end
    
    methods (Access = public)
        
        
        % ------ Constructor -----%
        function obj = SkeletonizationProcess(owner,funParams)
                                              
            if nargin == 0
                super_args = {};
            else                
                
                super_args{1} = owner;
                super_args{2} = 'Skeletonization';
                super_args{3} = @skeletonize3DMovieMasks;                               
                
                if nargin < 2 || isempty(funParams)                                       
                    
                    %----Defaults----%      
                                        
                    funParams.ChannelIndex = 1; 
                    funParams.GetGraph = true;
                    funParams.ClearBoundary = true;
                    funParams.OutputDirectory = ...
                        [owner.outputDirectory_  filesep 'skeletonization'];                    
                    funParams.BatchMode = false;                                                      
                                    
                end
                
                super_args{4} = funParams;    
                                
            end
            
            obj = obj@ImageProcessingProcess(super_args{:});
        end
        
        % ------ Set / Get Methods ----- %
        
        function setOutGraphPath(obj,iChan,graphDir)
            
            if nargin < 3 || isempty(iChan) || isempty(graphDir)
                error('You must input a channel number and directory!')
            elseif ~exist(graphDir,'dir')
                error('The directory specified does not exist!')
            end
            
            %Second row is skeleton graph directories.
            obj.outFilePaths_{2,iChan} = graphDir;
        
        
        end
       
        % ----- Output Check/Load Methods ----- %
         
        %overload the default to support 3D images.
        function outIm = loadOutImage(obj,iChan,iFrame)
            
            if nargin < 3 || isempty(iChan) || isempty(iFrame)
                error('You must specify a frame and channel number!')
            end            
            
            if length(iChan) > 1 || length(iFrame) > 1
                error('You can only specify 1 image to load!')
            end
            
            if ~obj.checkFrameNum(iFrame)
                error('Invalid frame number!')
            end
            
            %get the image names
            imNames = getOutImageFileNames(obj,iChan);
            
            outIm = stackRead([obj.outFilePaths_{1,iChan} ...
                filesep imNames{1}{iFrame}]);
            
        end
        
      
        
        function OK = checkChannelSkeletonGraphs(obj,iChan)
            
            if nargin < 2 || isempty(iChan)
                iChan = 1:numel(obj.owner_.channels_);
            elseif iChan < 1 || iChan > numel(obj.owner_.channels_) || numel(iChan) > 1
                error('Invalid channel number!')
            end          
            nChanCheck = numel(iChan);
            OK = false(1,nChanCheck);
            for j = 1:nChanCheck
                if size(obj.outFilePaths_,1) > 1 && ...
                        ~isempty(obj.outFilePaths_{2,iChan}) && ...
                        numel(dir([obj.outFilePaths_{2,iChan}  filesep '*.mat'])) ...
                            == obj.owner_.nFrames_;

                    OK(j) = true;
                else
                    OK(j) = false;
                end                                    
            end
        end
        
        
        function skelGraph = loadSkeletonGraph(obj,iChan,iFrame)
            
            if nargin < 3 || isempty(iChan) || isempty(iFrame)
                error('You must specify a channel number and frame number!')                
            elseif ~obj.checkChanNum(iChan)
                error('Invalid channel number!')
            elseif iFrame < 1 || iFrame > obj.owner_.nFrames_
                error('Invalid frame number!')
            end
            
            if ~obj.checkChannelSkeletonGraphs(iChan);
                error('The specified channel does not have valid skeleton graph files!')
            end
                        
            skelNames = dir([obj.outFilePaths_{2,iChan} filesep '*.mat']);
            
            skelGraph = load([obj.outFilePaths_{2,iChan} filesep skelNames(iFrame).name]);
            
            fNames = fieldnames(skelGraph);
            
            if numel(fNames) ~= 3 || ~isfield(skelGraph,'vertices') || ...
                    ~isfield(skelGraph,'edges') || ...
                    ~isfield(skelGraph,'edgePaths')
                error('Invalid skeleton graph file!')                                            
            end
                
            
        end
            
    end
    
    
end
        
    
    