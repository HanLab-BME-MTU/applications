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
            super_args{2} = 'Ratioing';
            super_args{3} = @ratioMovie;                
            
            if nargin < 3 || isempty(funParams)
                
                %----Defaults----%      
                funParams.OutputDirectory = ...
                    [outputDir  filesep 'ratio_images'];
                funParams.ChannelIndex = [];                
                funParams.ApplyMasks = true;
                funParams.SegProcessIndex = []; %No default
                funParams.CreateMasks = false;
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
    methods (Static)
        function text = getHelp(all)
            %Note: This help is designed for the GUI, and is a simplified
            %and shortened version of the help which can be found in the
            %function.
            if nargin < 1  % Static method does not have object as input
                all = false;
            end
            description = 'This is the most important step, where the fully-corrected activation channel (usually FRET) is divided by the fully-corrected localization channel (usually the FRET donor). This allows changes in the activation of the protein of interest to be separated from changes in its localization, giving a measure of the protein''s signaling state.';
            paramList = {'Apply Masks to Ratio',...
                         'Create New Masks',...
                         'Numerator',...
                         'Denominator',...
                         'Numerator Mask',...
                         'Denominator Mask'};
                         
                         
            paramDesc = {'If this box is checked, any areas in the ratio image which are outside of the masks from either the numerator or denominator channels will be set to zero. This is normally desireable because areas outside the mask are background areas, and the ratio values there are very noisy and generally meaningless.',...
                         'If this box is checked, new masks will be created and saved for the ratio channel, which are combinations of the masks from the two input channels: Only pixels which are included in both masks will be included in these masks. This is generally only useful if the ratio images are NOT being masked.',...
                         'This allows you to select which channel is the numerator, or activity channel, in the ratio. This is usually the FRET channel.',...
                         'This allows you to select which channel is the denominator, or localization channel, in the ratio. This is usually the FRET donor channel.',...
                         'This allows you to select which channel to use masks from for the numerator. Normally this is the numerator channel itself, but in the case that the numerator has poor-quality masks, a different channel can be used. As long as the channels are well-aligned, this is an acceptable alternative.',...
                         'This allows you to select which channel to use masks from for the denominator. Normally this is the denominator channel itself, but in the case that the denominator has poor-quality masks, a different channel can be used. As long as the channels are well-aligned, this is an acceptable alternative.'};
            if all
                text = makeHelpText(description,paramList,paramDesc);
            else
                text = makeHelpText(description);
            end
             
        end
                
    end    
    
    
    
    
end
        
        