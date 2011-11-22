classdef FlowTrackingProcess < ImageAnalysisProcess
    % Concrete class for a flow tracking process
    %
    % Sebastien Besson, 5/2011
    
    methods
        function obj = FlowTrackingProcess(owner,varargin)
            
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
                super_args{2} = FlowTrackingProcess.getName;
                super_args{3} = @trackMovieFlow;
                if isempty(funParams)
                   funParams=FlowTrackingProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            obj = obj@ImageAnalysisProcess(super_args{:});
        end
        
        function OK = checkChannelOutput(obj,varargin)
           % Input check
           ip =inputParser;
           ip.addOptional('iChan',1:numel(obj.owner_.channels_),...
               @(x) all(ismember(x,1:numel(obj.owner_.channels_))));
           ip.parse(varargin{:});
           iChan=ip.Results.iChan;

           %Makes sure there's at least one output file per channel
           OK =  arrayfun(@(x) exist(obj.outFilePaths_{1,x},'file'),iChan);
        end
        
        function varargout = loadChannelOutput(obj,iChan,varargin)
            
            % Input check
            outputList = {'flow','corLen'};
            ip =inputParser;
            ip.addRequired('obj',@(x) isa(x,'FlowTrackingProcess'));
            ip.addRequired('iChan',@(x) ismember(x,1:numel(obj.owner_.channels_)));
            ip.addOptional('iFrame',1:obj.owner_.nFrames_,...
                @(x) all(ismember(x,1:obj.owner_.nFrames_)));
            ip.addParamValue('output',outputList,@(x) all(ismember(x,outputList)));
            ip.parse(obj,iChan,varargin{:})
            iFrame = ip.Results.iFrame;
            output = ip.Results.output;
            if ischar(output), output={output}; end
            
            % Initialize output
            for j=1:numel(output)
                if numel(iFrame)==1
                    varargout{j}=[];
                else
                    varargout{j} = cell(size(iFrame));
                end
            end
            
            % Fill output
            for i=1:numel(iFrame)
                flowFile= [obj.outFilePaths_{1,iChan}...
                    filesep 'flow' num2str(iFrame(i)) '.mat'];
                % Check existence of output file
                if ~exist(flowFile,'file'), continue; end
                
                % Read variables and dispatch them
                s = load(flowFile,output{:});
                for j=1:numel(output)
                    if numel(iFrame)==1,
                        varargout{j}=s.(output{j});
                    else
                        varargout{j}{i}=s.(output{j});
                    end
                end
            end    
        end
        
        function checkValue=checkValue(obj,property,value)
            % Test the validity of a property value
            
            switch property
                case {'lastImage','timeWindow','timeStepSize'}
                    checkTest=@(x) isnumeric(x) && x>0 && x<=obj.owner_.nFrames_;
                case {'firstImage','minCorLength', 'maxCorLength',...
                        'minFeatureSize','edgeErodeWidth',...
                        'numStBgForAvg','maxFlowSpeed'}
                    checkTest=@(x) isnumeric(x) && x>0;
            end
            checkValue = checkTest(value);
        end
        
        function output = getDrawableOutput(obj)
            colors = hsv(numel(obj.owner_.channels_));
            output(1).name='Flow vectors';
            output(1).var='flow';
            output(1).formatData=@formatFlow;
            output(1).type='overlay';
            output(1).defaultDisplayMethod=@(x) VectorFieldDisplay('Color',colors(x,:));
            output(2).name='Template sizes';
            output(2).var='corLen';
            output(2).formatData=@formatBoxes;
            output(2).type='overlay';
            output(2).defaultDisplayMethod=@(x) RectangleDisplay('Color',colors(x,:));
        end
        
    end
    methods (Static)
        function name =getName()
            name = 'Flow Tracking';
        end
        function h = GUI()
            h= @flowTrackingProcessGUI;
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex = 1:numel(owner.channels_);
            funParams.OutputDirectory = [outputDir  filesep 'flow'];
            funParams.firstImage = 1;
            funParams.lastImage = owner.nFrames_;
            funParams.timeWindow = 2;
            funParams.timeStepSize = 1;
            funParams.minCorLength = 11;
            funParams.maxCorLength = 11;
            funParams.numStBgForAvg = 0;
            funParams.minFeatureSize = 11;
            funParams.edgeErodeWidth =5;
            funParams.maxFlowSpeed =10;
            funParams.outlierThreshold = 2;
            funParams.ROI = [1 1;owner.imSize_];
        end
    end
end

function flow = formatFlow(initFlow)
if isempty(initFlow),
    flow=zeros(1,4);
else
    flow=[initFlow(:,[2 1]) initFlow(:,[4 3])-initFlow(:,[2 1])];
end
end

function data = formatBoxes(initData)
if isempty(initData),
    data=[];
else
    finiteLength = isfinite(initData(:,3));
    data=[initData(finiteLength,1:2) repmat(initData(finiteLength,3)/2,1,2)];

end
end

