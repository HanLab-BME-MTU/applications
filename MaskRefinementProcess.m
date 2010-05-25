classdef RefineMaskProcess < SegmentationProcess
    %INHERITANCE OF SEGMENTATIONPROCESS IS TEMPORARY - NEED TO MOVE TO OWN
    %CLASS. HLE  

    methods(Access = public)        
        function obj = RefineMaskProcess(owner,funParams)
            
            if nargin == 0
                super_args = {};
            else
                nChan = numel(owner.channelPath_);
                
                super_args{1} = owner;
                super_args{2} = @refineMovieMasksNEW;                               
                
                if nargin < 2 || isempty(funParams)                                       
                    
                    %----Defaults----%                                            
                    funParams.ChannelIndex = 1:nChan; %Default is to attempt to refine masks for all channels
                    funParams.MaskCleanUp = true;
                    funParams.MinimumSize = 10;
                    funParams.ClosureRadius = 3;
                    funParams.ObjectNumber = 1; %Default is to keep only 1 object per mask
                    funParams.FillHoles = true;
                    funParams.EdgeRefinement = false; %This off by default because it sort of sucks, and is slow.
                    funParams.MaxEdgeAdjust = []; %Use refineMaskEdges.m function defaults for these settings
                    funParams.MaxEdgeGap = [];    
                    funParams.PreEdgeGrow = [];
                    funParams.BatchMode = false;                                              
                end
                %Make sure the input parameters are legit??
                super_args{3} = funParams;                    
            end
            
            obj = obj@SegmentationProcess(super_args{:});
            obj.name_ = mfilename('class');
            
        end        
    end
end
    