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
                funParams.MakeMovie=0;
                funParams.MovieOptions.Saturate=0;
                funParams.MovieOptions.ConstantScale=0;
                funParams.MovieOptions.ColorBar=1;
                funParams.MovieOptions.MakeAvi=0;
                funParams.MovieOptions.MakeMov=1;
            end
            
            super_args{4} = funParams;    
                
            if nargin > 3
                super_args{5} = inImagePaths;
            end
            if nargin > 4
                super_args{6} = outImagePaths;
            end
                                                        
            obj = obj@DoubleProcessingProcess(super_args{:});
            obj.setFunc_ = @outputRatioProcessGUI; % FOr analyzability/ to be implemented
        end
        
        function sanityCheck(obj)
        end
        function OK = checkChannelOutput(obj,iChan)
            
           %Checks if the selected channels have valid output images          
           nChanTot = numel(obj.owner_.channels_);
           if nargin < 2 || isempty(iChan)
               iChan = 1:nChanTot;
           end
           
           OK =  arrayfun(@(x)(x <= nChanTot && ...
                             x > 0 && isequal(round(x),x) && ...
                             (length(dir([obj.outFilePaths_{1,x} filesep '*.tif']))...
                             == obj.owner_.nFrames_)),iChan);
        end   
        function fileNames = getOutImageFileNames(obj,iChan)
            if obj.checkChannelOutput(iChan)
                fileNames = cellfun(@(x)(imDir(x)),obj.outFilePaths_(1,iChan),'UniformOutput',false);
                fileNames = cellfun(@(x)(arrayfun(@(x)(x.name),x,'UniformOutput',false)),fileNames,'UniformOutput',false);
                nIm = cellfun(@(x)(length(x)),fileNames);
                if ~all(nIm == obj.owner_.nFrames_)                    
                    error('Incorrect number of images found in one or more channels!')
                end                
            else
                error('Invalid channel numbers! Must be positive integers less than the number of image channels!')
            end    
            
            
        end
        function figHan = resultDisplay(obj)
                        
            figHan = msgbox(['The ratio images have been multiplied by a scale factor and saved as .tif images to the folder "' ...
                obj.funParams_.OutputDirectory '". They can be viewed with ImageJ or a comparable image viewing program.']); 
            
        end
    
    end
end
   