classdef SkeletonizationProcess < ImageProcessingProcess
% A process class for performing skeletonization on a movie
%     
% Hunter Elliott
% 12/2010    
%
    properties (SetAccess = protected, GetAccess = public)
        
        
    end
    
    methods (Access = public)
        
        function obj = SkeletonizationProcess(owner,funParams)
                                              
            if nargin == 0
                super_args = {};
            else                
                
                super_args{1} = owner;
                super_args{2} = 'Skeletonization';
                super_args{3} = @skeletonizeMovie;                               
                
                if nargin < 2 || isempty(funParams)                                       
                    
                    %----Defaults----%      
                                        
                    funParams.ChannelIndex = 1; 
                    funParams.GetGraph = false;
                    funParams.ClearBoundary = true;
                    funParams.OutputDirectory = ...
                        [owner.outputDirectory_  filesep 'skeletonization'];                    
                    funParams.BatchMode = false;                                                      
                                    
                end
                
                super_args{4} = funParams;    
                                
            end
            
            obj = obj@ImageProcessingProcess(super_args{:});
        end
                        
    end
    
end
        
    
    