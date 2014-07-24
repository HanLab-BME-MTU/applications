classdef SkeletonTrackingProcess < ImageAnalysisProcess
%A process class for tracking pruned mask skeletons throughout a 3D movie
%using track3DMovieSkeleton.m
%
% Hunter Elliott
% 7/2012

    properties (SetAccess = protected, GetAccess = public)
        
        
    end
    
    methods (Access = public)
    
     % ------ Constructor -----%
        function obj = SkeletonTrackingProcess(owner,funParams)
                                              
            if nargin == 0
                super_args = {};
            else                
                
                super_args{1} = owner;
                super_args{2} = 'SkeletonTracking';
                super_args{3} = @track3DMovieSkeleton;                               
                
                if nargin < 2 || isempty(funParams)                                                           
                    funParams = SkeletonTrackingProcess.getDefaultParams(owner);
                end
                
                super_args{4} = funParams;                                    
            end
            
            obj = obj@ImageAnalysisProcess(super_args{:});
        end
        
        % ------ Set / Get Methods ----- %
        
        %Uses generic ImageAnalysisProcess methods
        
        % ----- Output Check/Load Methods ----- %
                 
        function status = checkChannelOutput(obj,iChan)
            
            if nargin < 2 || isempty(iChan) || ~obj.checkChanNum(iChan)
                error('You must specify a valid channel number!')
            end
                        
            %Not channel-specific, only one file/channel
            status =  exist(obj.outFilePaths_{iChan},'file') ~= 0;
            
        end
                              
        function skelTracks = loadChannelOutput(obj,iChan)
            
            if nargin < 2 || isempty(iChan)
                error('You must specify a channel number and frame number!')                
            elseif ~obj.checkChanNum(iChan)
                error('Invalid channel number!')            
            end
            
            if ~obj.checkChannelOutput(iChan);
                error('The specified channel does not have a valid skeleton tracking file!')
            end
                                                
            skelTracks = load(obj.outFilePaths_{iChan});
            skelTracks = skelTracks.skelTracks;            
            
        end
        
        function sanityCheck(obj)
            
        end
    end
    
    methods(Static)
        
        function getName
        end
        
        function h = GUI()   
            h = [];
        end
        function procClasses = getConcreteClasses()
            procClasses = [];
        end        

        function funParams = getDefaultParams(owner)
            %----Defaults----%

            funParams.ChannelIndex = [];%Still using stupid hard-coded processingChannel for non-channel specific processing
            funParams.SkelProcessIndex = [];%No default 
            funParams.TrackParam.MaxDisp = 20e3*60;%Maximum displacement, in nm/s (chaosen so it would allow all displacements from bobs hand-tracking results)
            funParams.OutputDirectory = ...
                [owner.outputDirectory_  filesep 'pruned_skeleton_tracking'];                    
            funParams.BatchMode = false;                                                      

        end

    end
    
end