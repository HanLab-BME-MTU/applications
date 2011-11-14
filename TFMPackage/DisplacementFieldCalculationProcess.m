classdef DisplacementFieldCalculationProcess < ImageAnalysisProcess
    % Concrete class for a displacement field calculation process
    %
    % Sebastien Besson, Aug 2011
    
    methods
        function obj = DisplacementFieldCalculationProcess(owner,varargin)
            
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
                super_args{2} = DisplacementFieldCalculationProcess.getName;
                super_args{3} = @calculateMovieDisplacementField;
                if isempty(funParams)
                    funParams=DisplacementFieldCalculationProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4}=funParams;
            end
            
            obj = obj@ImageAnalysisProcess(super_args{:});
            
        end
        function sanityCheck(obj)
            
        end
        
        function status = checkChannelOutput(obj,varargin)
            status = logical(exist(obj.outFilePaths_{1},'file'));
        end
        
        function varargout = loadChannelOutput(obj,varargin)
            
            outputList = {'displField'};
            ip =inputParser;
            ip.addRequired('obj',@(x) isa(x,'DisplacementFieldCalculationProcess'));
            ip.addOptional('iFrame',1:obj.owner_.nFrames_,...
                @(x) ismember(x,1:obj.owner_.nFrames_));
            ip.addParamValue('output',outputList{1},@(x) all(ismember(x,outputList)));
            ip.parse(obj,varargin{:})
            iFrame = ip.Results.iFrame;
            
            % Data loading
            output = ip.Results.output;
            if ischar(output), output = {output}; end
            s = load(obj.outFilePaths_{1},output{:});
            
            if numel(iFrame)>1,
                for i=1:numel(output),
                    varargout{i}=s.(output{i});
                end
            else
                for i=1:numel(output),
                    varargout{i}=s.(output{i})(iFrame);
                end
            end
        end
        
        function h=draw(obj,iFrame,varargin)
            % Function to draw process output (template method)
            
            if ~ismember('getDrawableOutput',methods(obj)), h=[]; return; end
            outputList = obj.getDrawableOutput();
            ip = inputParser;
            ip.addRequired('obj',@(x) isa(x,'Process'));
            ip.addRequired('iFrame',@isnumeric);
            ip.addParamValue('output',outputList(1).var,@(x) any(cellfun(@(y) isequal(x,y),{outputList.var})));
            ip.KeepUnmatched = true;
            ip.parse(obj,iFrame,varargin{:})
            
            data=obj.loadChannelOutput(iFrame,'output',ip.Results.output);
            iOutput= find(cellfun(@(y) isequal(ip.Results.output,y),{outputList.var}));
            if ~isempty(outputList(iOutput).formatData),
                data=outputList(iOutput).formatData(data);
            end
            try
                assert(~isempty(obj.displayMethod_{iOutput}));
            catch ME
                obj.displayMethod_{iOutput}=...
                    outputList(iOutput).defaultDisplayMethod();
            end
            
            % Delegate to the corresponding method
            tag = [obj.getName '_output' num2str(iOutput)];
            drawArgs=reshape([fieldnames(ip.Unmatched) struct2cell(ip.Unmatched)]',...
                2*numel(fieldnames(ip.Unmatched)),1);
            h=obj.displayMethod_{iOutput}.draw(data,tag,drawArgs{:});
        end
        
    end
    methods (Static)
        function name =getName()
            name = 'Displacement Field Calculation';
        end
        function h = GUI()
            h= @displacementFieldCalculationProcessGUI;
        end
        
        function output = getDrawableOutput()
            output(1).name='Displacement field';
            output(1).var='displField';
            output(1).formatData=@(x) [x.pos x.vec];
            output(1).type='movieOverlay';
            output(1).defaultDisplayMethod=@(x) VectorFieldDisplay('Color','r');
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex = 1;
            funParams.OutputDirectory = [outputDir  filesep 'displacementField'];
            funParams.referenceFramePath='';
            funParams.I0=[];
            funParams.sDN=[];
            funParams.GaussRatio=[];
            funParams.alpha=.05;
            funParams.minCorLength = 21;
            funParams.maxFlowSpeed =20;
        end
    end
end