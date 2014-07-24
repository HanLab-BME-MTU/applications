classdef MaskedIntensity3DProcess < ImageAnalysisProcess
%MASKEDINTENSITY3DPROCESS process for analyzing 3D masked intensity with analyze3DMovieMaskedIntensities.m
%
%
% Hunter Elliott
% 1/2012
%

%NOT FINISHED YET!! FINISH THIS SHIT!!

    properties (SetAccess = protected, GetAccess = public)
        
        
    end
    
    methods (Access = public)
        
        function obj = MaskedIntensity3DProcess(owner,funParams)
                                              
            if nargin == 0
                super_args = {};
            else                                                
                
                super_args{1} = owner;
                super_args{2} = 'MaskedIntensities';
                super_args{3} = @analyze3DMovieMaskedIntensities;                               
                
                if nargin < 2 || isempty(funParams)                                       
                    
                    funParams = MaskedIntensity3DProcess.getDefaultParams(owner);                 
                                    
                end
                
                super_args{4} = funParams;    
                                
            end
            
            obj = obj@ImageAnalysisProcess(super_args{:});
        end
        function OK = checkChannelOutput(obj,iChan)
            
           %Checks if the selected channel has valid output files           
           if nargin < 2 || isempty(iChan) || numel(iChan) > 1 || ... 
                   ~all(obj.checkChanNum(iChan));
               error('You must specify a valid channel number!')
           end
           OK =  arrayfun(@(x)(exist(obj.outFilePaths_{x},'file')),iChan);
        end
        
        
        function mg = loadChannelOutput(obj,iChan,iFrame)
                                   
            mg = load(obj.outFilePaths_{1});
                        
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

                    funParams.ChannelIndex = 1:numel(owner.channels_);                                        
                    funParams.OutputDirectory = ...
                        [owner.outputDirectory_  filesep 'masked_intensity_analysis'];
                    funParams.CurvSampRad = 1000;% In nanometers.
                    funParams.PhotoBleachMeth = 'Self';%Default is to use the PB correct from the movie being processed                    
                    funParams.TrendRemoval = 'Linear';
                    funParams.BatchMode = false;                                                   
        end
    end            
    
end