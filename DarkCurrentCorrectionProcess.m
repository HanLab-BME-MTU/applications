classdef DarkCurrentCorrectionProcess < ImageCorrectionProcess
    
    %A class for performing dark current correction on images.
    %
    %Hunter Elliott, 5/2010
    %
        
    methods (Access = public)
        
        function obj = DarkCurrentCorrectionProcess(owner,outputDir,funParams,darkImagePaths,...
                                              inImagePaths,outImagePaths)
            
            if nargin == 0
                super_args = {};
            else
                nChan = numel(owner.channels_);
                
                super_args{1} = owner;
                super_args{2} = 'Dark Current Correction';
                super_args{3} = @darkCurrentCorrectMovie;                               
                
                if nargin < 3 || isempty(funParams)                                       
                    
                    %----Defaults----%      
                    funParams.OutputDirectory = ...
                        [outputDir  filesep 'dark_current_corrected_images'];  
                    funParams.DarkImageDirectories = []; %No default for this! It will be handled differently...
                    funParams.ChannelIndex = 1:nChan;
                    funParams.MedianFilter = true;
                    funParams.GaussFilterSigma = 0;                    
                    funParams.BatchMode = false;                                                                                
                                       
                end
                
                super_args{4} = funParams;    
                
                if nargin > 3           
                    %Set the correction image paths to the darkImage paths
                    %input.
                    super_args{7} = darkImagePaths;                
                end
                
                if nargin > 4
                    super_args{5} = inImagePaths;
                end
                
                if nargin > 5
                    super_args{6} = outImagePaths;
                end                
                
            end
            
            obj = obj@ImageCorrectionProcess(super_args{:});
        end   
        
        function sanityCheck(obj)
        % Sanity check will check the correction images
            for i = obj.funParams_.ChannelIndex
                if ~isempty(obj.correctionImagePaths_{i})
                    
                    if ~exist(obj.correctionImagePaths_{i}, 'dir')
                        error('lccb:set:fatal', ...
                            ['The specified shade image channel:\n\n ',obj.correctionImagePaths_{i}, ...
                            '\n\ndoes not exist. Please double check your channel path.'])
                    end
                    fileNames = imDir(obj.correctionImagePaths_{i},true);
                    
                    if isempty(fileNames)
                        error('lccb:set:fatal', ...
                            ['No proper image files are detected in:\n\n ',obj.correctionImagePaths_{i}, ...
                            '\n\nPlease double check your channel path.'])                        
                    end
                    
                    for j = 1:length(fileNames)
                        imInfo = imfinfo([obj.correctionImagePaths_{i} filesep fileNames(j).name]);
                        if imInfo.Width ~= obj.owner_.imSize_(1) || imInfo.Height ~= obj.owner_.imSize_(2)
                            error('lccb:set:fatal', ...
                                ['Dark current correction image - \n\n',...
                                obj.correctionImagePaths_{i},filesep,fileNames(j).name,...
                                '\n\nmust have the same size as input images. Please double check your correction image data.'])
                        end
                    end
                end
            end
        end        
           
       
    end
    methods (Static)
        function text = getHelp(all)
            %Note: This help is designed for the GUI, and is a simplified
            %and shortened version of the help which can be found in the
            %function.
            if nargin < 1  % Static method does not have object as input
                all = false;
            end                                                                                                                                                                                             %double percent is to escape % formating
            description = '"Dark Current" is background noise in the camera used for imaging. This noise is generally a few hundred counts, and is not uniform throughout the image, sometimes varying by > 10%% within the imaging area of a single camera, and can vary even more between cameras.m. This noise is measured by taking several "dark-current images" - images where there is no light incident on the camera CCD. These images are then averaged together, and used to correct for the inhomogeneity of the dark-current by subtracting the dark-current value from each pixel. If a single camera is being used for all image acquisition, a single set of dark-current images can be used to correct all channels. If different cameras are used, it is best to collect and use dark images for each camera. The final, averaged (and possibly filtered) correction images, in addition to the corrected images, can be seen by clicking the "Result" buton.';
            paramList = {'Input Channels',...
                         'Dark-Current Image Channels',...
                         '3x3 Median Filter',...
                         'Gaussian Filter'};
                         
            paramDesc = {'This allows you to select which channels you want to perform dark-current correction on. This should be applied to all channels that are going to be used for ratioing or bleedthrough correction. Select the channels by clicking on them in the "Available Input Channels" box and then clicking "Select->" to move them to the "Selected Channels" box. You can un-select a channel by clicking the "Delete" button',...
                         'This box allows you to specify a directory containing the dark-current images corresponding to each channel to be corrected. You must specify a directory for each channel to be dark-current corrected, but the same directory may be specified multiple times. The directories specified should contain one or more "dark-current images". It is recommended to take 5 or more, as these will be averaged together to improve the quality of the correction. Using only a single correction image can introduce error into the final ratio images.',...
                         'If this box is checked, a median filter will be applied to the dark-current images prior to their use as a correction. This is usefull because it minizes the contribution of noise in the dark images, and removes "hot pixels" - pixels which have a much higher-than-normal background value (several hundred counts).',...
                         'If checked, the dark-current images will also be filtered (smoothed) using a gaussian filter whose sigma (in pixels) is specified by the value in the "Sigma" box. Larger sigmas will give smoother images, but if the sigma is too large, spatial information will be lost. A sigma of 1 is a good starting point. This can be used to further reduce noise in the dark-current images and is especially important if only 1 or a few dark-current image(s) are taken.'};

            if all
                text = makeHelpText(description,paramList,paramDesc);
            else
                text = makeHelpText(description);
            end
             
        end
    end    
end