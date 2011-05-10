classdef MaskObjectTrackingProcess < ImageAnalysisProcess
% A process class for detecting branches in movies with detectMovieBranches     
%     
% Hunter Elliott
% 7/2010    
%
    properties (SetAccess = protected, GetAccess = public) 
        
        
    end
    
    methods (Access = public)
        
        function obj = MaskObjectTrackingProcess(owner,funParams)
                                              
            if nargin == 0
                super_args = {};
            else                
                
                super_args{1} = owner;
                super_args{2} = 'ObjectTracking';
                super_args{3} = @track3DMaskObjects;                               
                
                if nargin < 2 || isempty(funParams)                                       
                    
                    %----Defaults----%      
                                        
                    funParams.ChannelIndex = 1;
                    funParams.OutputDirectory = ...
                        [owner.outputDirectory_  filesep 'object_tracking'];
                    funParams.BatchMode = false;                                                      
                                    
                end
                
                super_args{4} = funParams;    
                                
            end
            
            obj = obj@ImageAnalysisProcess(super_args{:});
        end
        
        function OK = checkChannelOutput(obj,iChan)
            
           %Checks if the selected channels have valid output files
           nChanTot = numel(obj.owner_.channels_);
           if nargin < 2 || isempty(iChan)
               iChan = 1:nChanTot;
           end
           %Makes sure there's at least one .mat file in the speified
           %directory
           OK =  arrayfun(@(x)(x <= nChanTot && ...
                             x > 0 && isequal(round(x),x) && ...
                             exist(obj.outFilePaths_{x},'file')),iChan);
                             
        end
                        
    end
    
end
        
    
    