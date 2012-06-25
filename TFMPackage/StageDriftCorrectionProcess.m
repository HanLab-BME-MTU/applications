classdef StageDriftCorrectionProcess < ImageProcessingProcess
    % Concrete class for a stage drift correction process
    %
    % Sebastien Besson, Sep 2011
    
    methods
        function obj = StageDriftCorrectionProcess(owner,varargin)
            
            if nargin == 0
                super_args = {};
            else
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieData'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                % Define arguments for superclass constructor
                super_args{1} = owner;
                super_args{2} = StageDriftCorrectionProcess.getName;
                super_args{3} = @correctMovieStageDrift;
                if isempty(funParams)
                    funParams=StageDriftCorrectionProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            obj = obj@ImageProcessingProcess(super_args{:});
        end
        
        function h=draw(obj,varargin)
            % Function to draw process output
            
            outputList = obj.getDrawableOutput();
            drawGraph = any(strcmp(varargin,'x-flow') | strcmp(varargin,'y-flow') |...
                strcmp(varargin,'refFrame'));
            
            
            if drawGraph
                % Use dedicated draw method for plotting flow histograms
                ip = inputParser;
                ip.addRequired('obj');
                ip.addParamValue('output',outputList(2).var,@(x) all(ismember(x,{outputList.var})));
                ip.KeepUnmatched = true;
                ip.parse(obj,varargin{:})
                
                [~,iOutput] =ismember(ip.Results.output,{outputList.var});
                if regexp(outputList(iOutput).var,'(.+)flow','once')
                    s=load(obj.outFilePaths_{3,1},'flow');
                    data=s.flow;
                else
                    data=imread(obj.funParams_.referenceFramePath);
                end
                
                if ~isempty(outputList(iOutput).formatData),
                    data=outputList(iOutput).formatData(data);
                end
                
                try
                    assert(~isempty(obj.displayMethod_{iOutput,1}));
                catch ME
                    obj.displayMethod_{iOutput,1}=...
                        outputList(iOutput).defaultDisplayMethod();
                end
                
                % Delegate to the corresponding method
                tag = ['process' num2str(obj.getIndex) '_output' num2str(iOutput)];
                drawArgs=reshape([fieldnames(ip.Unmatched) struct2cell(ip.Unmatched)]',...
                    2*numel(fieldnames(ip.Unmatched)),1);
                h=obj.displayMethod_{iOutput}.draw(data,tag,drawArgs{:});
            else
                % Call superclass method
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
            output(2).name='Reference frame';
            output(2).var='refFrame';
            output(2).formatData=@mat2gray;
            output(2).type='movieGraph';
            output(2).defaultDisplayMethod=@ImageDisplay;
            output(3).name='Flow along x-axis';
            output(3).var='x-flow';
            output(3).formatData=@(x)getFlow(x,1);
            output(3).type='movieGraph';
            output(3).defaultDisplayMethod=@FlowHistogramDisplay;
            output(4).name='Flow along y-axis';
            output(4).var='y-flow';
            output(4).formatData=@(x)getFlow(x,2);
            output(4).type='movieGraph';
            output(4).defaultDisplayMethod=@FlowHistogramDisplay;
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
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
    end
end

function data=getFlow(data,i)

data=cellfun(@(x) x(:,i+2)-x(:,i),data,'UniformOutput',false);
end