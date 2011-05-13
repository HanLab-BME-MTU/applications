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
            obj.setFunc_ = @photobleachCorrectionProcessGUI; % FOr analyzability/ to be implemented

        end
        
        function figHan = resultDisplay(obj)
            
            %Open the figure
            figHan = open([obj.funParams_.OutputDirectory ...
                            filesep obj.funParams_.figName]);
            
        end
        
    end

    
end                                   
            