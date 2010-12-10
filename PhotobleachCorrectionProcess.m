classdef PhotobleachCorrectionProcess < DoubleProcessingProcess
    
    %A class for performing photobleach correction on ratio images.
    %
    %Hunter Elliott, 5/2010
    %

    methods (Access = public)
        
        function obj = PhotobleachCorrectionProcess(owner,outputDir,funParams)
                                              
            if nargin == 0
                super_args = {};
            else                
                
                super_args{1} = owner;
                super_args{2} = 'Photobleach Correction';
                super_args{3} = @photobleachCorrectMovieRatios;                               
                
                if nargin < 3 || isempty(funParams)                                       
                    
                    %----Defaults----%      
                    funParams.OutputDirectory = ...
                        [outputDir  filesep 'photobleach_corrected_images'];                      
                    funParams.ChannelIndex = [];%No default
                    funParams.CorrectionType = 'RatioOfAverages';
                    funParams.BatchMode = false;                                                                                
                                    
                end
                
                super_args{4} = funParams;    
                                
            end
            
            obj = obj@DoubleProcessingProcess(super_args{:});
        end   
        
        function figHan = resultDisplay(obj)
            
            %Open the figure
            figHan = open([obj.funParams_.OutputDirectory ...
                            filesep obj.funParams_.figName]);
            
            
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
            description = 'Photobleaching causes the fluorophores in each channel to lose intensity over time as they are imaged. The rate of photobleaching is almost always different in the different channels, and therefore can introduce artifacts in the ratio image if not corrected: The apparent average ratio will vary as the two fluorophores bleach at different rates. This process allows the photobleaching to be corrected by fitting a double-exponential to the intensities in each channel. After the fitting, a figure will be displayed showing the average ratio over time, and the fitted value along with confidence intervals. It is important to inspect this fit to ensure that the correction is valid, and that the confidence intervals are not excessively large. If photobleach correction is not applied, it is generally not valid to measure changes in the ratio images over time.';
            paramList = {'Ratio Channels',...
                         'Ratio of Averages',...
                         'Average of Ratios',...
                         'Ratio of Totals'};
                         
            paramDesc = {'This channel should be set to be the same as the numerator channel used in ratioing (usually FRET).',...
                         'This method fits a double-exponential to the ratio of the average intensity in each of the channels used for ratioing, and then uses this fit to correct the ratio images. This is standard method.',...
                         'This method fits a double-exponential to the average of the ratio image, and then uses this fit to correct the ratio images. This will suppress any long-term changes in the ratio, and should only be used if the other methods fail.',...
                         'This method fits a double-exponential to the ratio of the total intensity in each of the channels used for ratioing, and then uses this fit to correct the ratio images. This is generally only used when there is a drastic change in area of the object being imaged during the movie. (e.g. a cell rounds up, spreads etc.). This sort of change will have a smaller effect on the total intensity than the average intensity, making the total a better measure of bleaching in these cases.'};
            if all
                text = makeHelpText(description,paramList,paramDesc);
            else
                text = makeHelpText(description);
            end
             
        end
    end    
    
    
end                                   
            