classdef BackgroundMasksProcess < MaskProcessingProcess
    % A concrete class for creating background masks
    %
    
    methods(Access = public)
        
        function obj = BackgroundMasksProcess(owner,outputDir,funParams)
            if nargin == 0
                super_args = {};
            else
                nChan = numel(owner.channels_);
                
                super_args{1} = owner;
                super_args{2} = BackgroundMasksProcess.getName;
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
    methods(Static)
        function name =getName()
            name = 'Background Mask';
        end
        function h = GUI()
            h= @backgroundMasksProcessGUI;
        end
    end
end