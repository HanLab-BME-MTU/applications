classdef RatioProcess < DoubleProcessingProcess
    
    %A class for creating ratios by dividing one channel by another using
    %ratioMovie.m
    %
    %Hunter Elliott,
    %6/2010
    
    methods (Access = public)
        
        function obj = RatioProcess(owner,outputDir,funParams,...                                              
                                    inImagePaths,outImagePaths)
            
                                
            
            super_args{1} = owner;
            super_args{2} = RatioProcess.getName;
            super_args{3} = @ratioMovie;                
            
            if nargin < 3 || isempty(funParams)
                
                %----Defaults----%      
                funParams.OutputDirectory = ...
                    [outputDir  filesep 'ratio_images'];
                funParams.ChannelIndex = [];                
                funParams.ApplyMasks = true;
                funParams.SegProcessIndex = []; %No default
                funParams.MaskChannelIndex = [];
                funParams.BatchMode = false;                
                
            end
            
            super_args{4} = funParams;    
                
            if nargin > 3
                super_args{5} = inImagePaths;
            end
            if nargin > 4
                super_args{6} = outImagePaths;
            end
                                                        
            obj = obj@DoubleProcessingProcess(super_args{:});
        end

                
    end
    methods(Static)
        function name =getName()
            name = 'Ratioing';
        end
        function h =GUI()
            h = @ratioProcessGUI;
        end
    end

end
        
        