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
        
        function h=draw(obj,varargin)
            % Function to draw process output
            
            outputList = obj.getDrawableOutput();
            drawFlow = any(strcmpi(outputList(2).var,varargin) |...
                strcmpi(outputList(3).var,varargin) );
            
            if drawFlow
                % Use dedicated draw method for plotting flow histograms
                ip = inputParser;
                ip.addRequired('obj');
                ip.addParamValue('output',outputList(2).var,@(x) any(cellfun(@(y) isequal(x,y),{outputList.var})));
                ip.KeepUnmatched = true;
                ip.parse(obj,varargin{:})
                
                iOutput= find(cellfun(@(y) isequal(ip.Results.output,y),{outputList.var}));
                data=obj.outFilePaths_{3+iOutput,:};

                try
                    assert(~isempty(obj.displayMethod_{iOutput,1}));
                catch ME
                    obj.displayMethod_{iOutput,1}=...
                        outputList(iOutput).defaultDisplayMethod();
                end
                
                % Delegate to the corresponding method
                tag = [obj.getName '_output' num2str(iOutput)];
                drawArgs=reshape([fieldnames(ip.Unmatched) struct2cell(ip.Unmatched)]',...
                    2*numel(fieldnames(ip.Unmatched)),1);
                h=obj.displayMethod_{iOutput}.draw(data,tag,drawArgs{:});
            else
                h=draw@ImageProcessingProcess(obj,varargin{1},varargin{2},...
                    varargin{3:end});
            end
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
            output(2).name='Flow along x-axis';
            output(2).var='x-flow';
            output(2).formatData=[];
            output(2).type='movieGraph';
            output(2).defaultDisplayMethod=@FigFileDisplay;
            output(3).name='Flow along y-axis';
            output(3).var='y-flow';
            output(3).formatData=[];
            output(3).type='movieGraph';
            output(3).defaultDisplayMethod=@FigFileDisplay;
        end

    end
end

