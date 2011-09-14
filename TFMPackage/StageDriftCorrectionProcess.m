classdef StageDriftCorrectionProcess < ImageProcessingProcess
    % Concrete class for a stage drift correction process
    %
    % Sebastien Besson, Sep 2011
    
    methods
        function obj = StageDriftCorrectionProcess(owner,outputDir,funParams)
            
            if nargin == 0
                super_args = {};
            else
                super_args{1} = owner;
                super_args{2} = StageDriftCorrectionProcess.getName;
                super_args{3} = @correctMovieStageDrift;
    
                if nargin < 3 || isempty(funParams)
                    
                    %----Defaults----%
                    funParams.ChannelIndex = 1 : numel(owner.channels_);
                    funParams.OutputDirectory = [outputDir  filesep 'stageDriftCorrection'];
                    funParams.referenceFramePath = '';
                    funParams.minCorLength = 51;
                    funParams.maxFlowSpeed =5;
                    funParams.I0=[];
                    funParams.sDN=[];
                    funParams.GaussRatio=[];
                    funParams.alpha=.05;
                    funParams.cropROI=[1 1 owner.imSize_(end:-1:1)];
                    funParams.doPreReg=1;
                end
                super_args{4} = funParams;
            end
             obj = obj@ImageProcessingProcess(super_args{:});
        end
        function sanityCheck(obj)
            
        end

    end
    methods (Static)
        
        function name =getName()
            name = 'Stage Drift Correction';
        end
        function h = GUI()
            h= @stageDriftCorrectionProcessGUI;
        end
        function output = getDrawableOutput()
            output(1).name='Registered images';
            output(1).var='';
            output(1).formatData=@mat2gray;
            output(1).type='image';
            output(1).defaultDisplayMethod=@ImageDisplay;
        end

    end
end

