classdef ExternalRatioProcess < DoubleProcessingProcess
    
    %A class for inputing / handling ratio .mat files created externally.
    %
    %Hunter Elliott,
    %11/2014
    
    methods (Access = public)
        
        function obj = ExternalRatioProcess(owner,outputDir)
            
                                                        
                
%             if nargin > 3
%                 super_args{5} = inImagePaths;
%             end
%             if nargin > 4
%                 super_args{6} = outImagePaths;
%             end
                                                        
            obj = obj@DoubleProcessingProcess(owner,ExternalRatioProcess.getName,...
                @ExternalRatioProcess.dummyRun);%,ExternalRatioProcess.getDefaultParams(owner,outputDir))
        end

                
    end
    
    
    methods(Static)
        
        function dummyRun(obj,owner)
            %Do nothing!
            %TEMP - maybe check the images?            
        end
        
        function name =getName()
            name = 'Ratios';
        end
%         function h =GUI()
%             h = @ratioProcessGUI;
%         end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            funParams.OutputDirectory =  [outputDir  filesep 'ratio_images'];
%             % Set default parameters
%             
%             funParams.ChannelIndex = [];
%             funParams.ApplyMasks = true;
%             funParams.SegProcessIndex = []; %No default
%             funParams.MaskChannelIndex = [];
%             funParams.BatchMode = false;
        end
    end

end
        
        