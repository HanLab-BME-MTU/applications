classdef BackgroundMasksProcess < MaskProcessingProcess
  
    methods(Access = public)
        
        function obj = BackgroundMasksProcess(owner,outputDir,funParams)
            if nargin == 0
                super_args = {};
            else
                nChan = numel(owner.channels_);
                
                super_args{1} = owner;
                super_args{2} = 'Background Mask';
                super_args{3} = @createMovieBackgroundMasks;                               
                
                if nargin < 3 || isempty(funParams)                                       
                    
                    %----Defaults----%
                    funParams.OutputDirectory = ...
                        [outputDir  filesep 'BackgroundMasks'];                       
                    funParams.ChannelIndex = 1:nChan; %Default is to attempt to creat background masks for all channels
                    funParams.SegProcessIndex = []; %No default
                    funParams.GrowthRadius = 20;
                    funParams.BatchMode = false;                                              
                end
                %Make sure the input parameters are legit??
                super_args{4} = funParams;                    
            end
            
            obj = obj@MaskProcessingProcess(super_args{:});
            
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
            description = 'This process uses the (previously created) masks, which label the foreground (i.e. fluorescent objects of interest), to create masks which label the background. To ensure that no foreground is included in these masks, the original masks are "grown" by a user-specified amount. This makes sure that the background mask is far away from any foreground objects.';
            paramList = {'Input Channels',...
                         'Growth Radius'};                     
                         
            paramDesc = {'This allows you to select which channels you want to create background masks for. By default, all channels will have background masks created. These masks are generally used to perform background subtraction. Select the channels by clicking on them in the "Available Input Channels" box and then clicking "Select->" to move them to the "Selected Channels" box. You can un-select a channel by clicking the "Delete" button',...
                         'This parameter determines how much the masks will be "grown" before being used to create background masks. That is, any pixel which is further than this distance away from any mask objects will be considered background. If you have very little background area in your images, you can decrease this value to make sure that some background is included in the background masks.'};
            if all
                text = makeHelpText(description,paramList,paramDesc);
            else
                text = makeHelpText(description);
            end
             
        end
    end    
end