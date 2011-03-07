classdef OutputRatioProcess < DoubleProcessingProcess
 
    %A class for creating ratios by dividing one channel by another using
    %ratioMovie.m
    %
    %Hunter Elliott,
    %6/2010
    
    methods (Access = public)
        
        function obj = OutputRatioProcess(owner,outputDir,funParams,...                                              
                                          inImagePaths,outImagePaths)
            
                                
            
            super_args{1} = owner;
            super_args{2} = 'RatioOutput';
            super_args{3} = @outputMovieRatios;                
            
            if nargin < 3 || isempty(funParams)
                
                %----Defaults----%      
                funParams.OutputDirectory = ...
                    [outputDir  filesep 'ratio_tiffs'];
                funParams.ChannelIndex = [];
                funParams.ScaleFactor = 1000;
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
        
        function sanityCheck(obj)
        end
        
        function figHan = resultDisplay(obj)
                        
            figHan = msgbox(['The ratio images have been multiplied by a scale factor and saved as .tif images to the folder "' ...
                obj.funParams_.OutputDirectory '". They can be viewed with ImageJ or a comparable image viewing program.']); 
            
        end
        
    end
    
end    
    
    