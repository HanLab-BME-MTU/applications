classdef ForceFieldCalculationProcess < DataProcessingProcess
    % Concrete process for calculating a force field
    %
    % Sebastien Besson, Aug 2011
    
    methods
        function obj = ForceFieldCalculationProcess(owner,varargin)
            
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
                super_args{2} = ForceFieldCalculationProcess.getName;
                super_args{3} = @calculateMovieForceField;
                if isempty(funParams)
                    funParams=ForceFieldCalculationProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            
            obj = obj@DataProcessingProcess(super_args{:});
            
        end
        
        function status = checkChannelOutput(obj,varargin)
            
            status = logical(exist(obj.outFilePaths_{1},'file'));
            
        end
        
        function varargout = loadChannelOutput(obj,varargin)
            
            outputList = {'forceField'};
            ip =inputParser;
            ip.addRequired('obj',@(x) isa(x,'ForceFieldCalculationProcess'));
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
        
        function h=draw(obj,varargin)
            % Function to draw process output
            
            outputList = obj.getDrawableOutput();
            drawLcurve = any(strcmpi('lcurve',varargin));
            
            if drawLcurve
                ip = inputParser;
                ip.addRequired('obj',@(x) isa(x,'Process'));
                ip.addParamValue('output',outputList(1).var,...
                    @(x) any(cellfun(@(y) isequal(x,y),{outputList.var})));
                ip.KeepUnmatched = true;
                ip.parse(obj,varargin{:})
                data=obj.outFilePaths_{3,1};
            else
                % Input parser
                ip = inputParser;
                ip.addRequired('obj',@(x) isa(x,'Process'));
                ip.addRequired('iFrame',@isnumeric);
                ip.addParamValue('output',outputList(1).var,...
                    @(x) any(cellfun(@(y) isequal(x,y),{outputList.var})));
                ip.KeepUnmatched = true;
                ip.parse(obj,varargin{1},varargin{2:end})
                iFrame=ip.Results.iFrame;
                
                data=obj.loadChannelOutput(iFrame,'output',ip.Results.output);
            end
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
            tag = ['process' num2str(obj.getIndex) '_output' num2str(iOutput)];
            drawArgs=reshape([fieldnames(ip.Unmatched) struct2cell(ip.Unmatched)]',...
                2*numel(fieldnames(ip.Unmatched)),1);
            h=obj.displayMethod_{iOutput}.draw(data,tag,drawArgs{:});
        end
        
        function output = getDrawableOutput(obj)
            output(1).name='Force  field';
            output(1).var='forceField';
            output(1).formatData=@(x) [x.pos x.vec];
            output(1).type='movieOverlay';
            output(1).defaultDisplayMethod=@(x) VectorFieldDisplay('Color','r');
            if ~strcmp(obj.funParams_.solMethodBEM,'QR')
                output(2).name='Lcurve';
                output(2).var='lcurve';
                output(2).formatData=[];
                output(2).type='movieGraph';
                output(2).defaultDisplayMethod=@FigFileDisplay;
            end
        end
        
        
    end
    methods (Static)
        function name =getName()
            name = 'Force Field Calculation';
        end
        function h = GUI()
            h= @forceFieldCalculationProcessGUI;
        end
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.OutputDirectory = [outputDir  filesep 'forceField'];
            funParams.YoungModulus = 10000;
            funParams.PoissonRatio = .5;
            funParams.method = 'FastBEM';
            funParams.meshPtsFwdSol = 4096;
            funParams.regParam=1e-7;
            funParams.solMethodBEM='QR';
            funParams.basisClassTblPath='';
            funParams.LcurveFactor=10;
        end
    end
end