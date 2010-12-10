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
    methods (Static)
        function text = getHelp(all)
            %Note: This help is designed for the GUI, and is a simplified
            %and shortened version of the help which can be found in the
            %function.
            if nargin < 1  % Static method does not have object as input
                all = false;
            end
            description = 'This process allows you to export the ratio images you have created in step 9 (or step 10, if photobleach correction has been applied) as .tif images. This requires multiplying the images by a large number, called a "scale factor". This is necessary because the ratio images themselves are small, non-integer values (usually between .5 and 5), while .tif images can only store positive integer values. If these ratios are converted to integers without multiplying by a scale factor, there would be severe rounding error in the images. Alternatively, you can use the ratio images themselves, which are stored as floating-point matlab .mat files in the movie''s output directory.';
            paramList = {'Ratio Channel',...
                         'Scale Factor',...
                         'Select Path'};
                         
                         
            paramDesc = {'This box allows you to select the channel which is the NUMERATOR of the ratio images you want to export (usually the FRET channel).',...
                         'This is the number that the ratio images will be multiplied by before being saved as .tif images. It is recommended that this number be fairly large (1000 is a good starting point) to minimize rounding error. It is also important that this number be kept constant among different experiments, if the resulting ratio images are to be compared.',...
                         'This allows you to specify the directory to save the .tif ratio images to. They will be saved to a sub-directory of this folder, called "ratio_tiffs", with one .tif file per ratio image.'};
            if all
                text = makeHelpText(description,paramList,paramDesc);
            else
                text = makeHelpText(description);
            end
             
        end
    end    
    
    
end    
    
    