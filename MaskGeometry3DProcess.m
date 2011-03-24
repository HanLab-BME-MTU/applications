classdef MaskGeometry3DProcess < ImageAnalysisProcess
%MASKGEOMETRY3DPROCESS process for analyzing 3D mask geometry with analyze3DMovieMaskGeometry.m
%
%
% Hunter Elliott
% 3/2011
%

    properties (SetAccess = protected, GetAccess = public)
        
        
    end
    
    methods (Access = public)
        
        function obj = MaskGeometry3DProcess(owner,funParams)
                                              
            if nargin == 0
                super_args = {};
            else                                                
                
                super_args{1} = owner;
                super_args{2} = 'MaskGeometry';
                super_args{3} = @analyze3DMovieMaskGeometry;                               
                
                if nargin < 2 || isempty(funParams)                                       
                    
                    %----Defaults----%      
                                        
                    funParams.ChannelIndex = 1;
                    funParams.SmoothSigma =[];%Use the analyze3DMaskGeometry.m defaults.
                    funParams.IsoValue =[];%Use the analyze3DMaskGeometry.m defaults.                    
                    funParams.OutputDirectory = ...
                        [owner.outputDirectory_  filesep 'mask_geometry'];
                    funParams.BatchMode = false;                                                      
                                    
                end
                
                super_args{4} = funParams;    
                                
            end
            
            obj = obj@ImageAnalysisProcess(super_args{:});
        end
        
    end
    
end