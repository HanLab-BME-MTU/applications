classdef FlowTrackingProcess < DataProcessingProcess
    % Concrete class for a flow tracking process
    %
    % Sebastien Besson, 5/2011
    
    methods
        function obj = FlowTrackingProcess(owner,outputDir, funParams)
            
            if nargin == 0
                super_args = {};
            else
                super_args{1} = owner;
                super_args{2} = FlowTrackingProcess.getName;
                super_args{3} = @trackMovieFlow;
                if nargin < 3 || isempty(funParams)
                    
                    %----Defaults----%
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
                super_args{4} = funParams;
            end
            
            obj = obj@DataProcessingProcess(super_args{:});
        end
        function sanityCheck(obj)
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
            %
            % INPUT:
            %    property - a property name (string)
            %    value - the property value to be checked
            %
            % OUTPUT:
            %    checkValue - a boolean containing the result of the test
            %
            % Sebastien Besson, 5/2011
            
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
        end
        
    end
    methods (Static)
        function name =getName()
            name = 'Flow Tracking';
        end
        function h = GUI()
            h= @flowTrackingProcessGUI;
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