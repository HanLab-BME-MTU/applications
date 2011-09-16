classdef ForceFieldCalculationProcess < Process
    % Concrete class for a force field calculation process
    %
    % Sebastien Besson, Aug 2011
    
    methods
        function obj = ForceFieldCalculationProcess(owner,outputDir,funParams)
            
            if nargin == 0
                super_args = {};
            else
                super_args{1} = owner;
                super_args{2} = ForceFieldCalculationProcess.getName;
            end
            
            obj = obj@Process(super_args{:});
            obj.funName_ = @calculateMovieForceField;
            
            %----Defaults----%
            defaultParams.OutputDirectory = [outputDir  filesep 'forceField'];
            defaultParams.YoungModulus = 10000;
            defaultParams.PoissonRatio = .5;
            defaultParams.method = 'FastBEM';
            defaultParams.meshPtsFwdSol = 4096;
            defaultParams.regParam=1e-7;
            defaultParams.solMethodBEM='QR';
            defaultParams.basisClassTblPath='';
            
            if nargin < 3 || isempty(funParams)
                obj.funParams_=defaultParams;
            else
                obj.funParams_=parseProcessParams(defaultParams,funParams);
            end

        end
        function sanityCheck(obj)
            
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
            name = 'Force Field Calculation';
        end
        function h = GUI()
            h= @forceFieldCalculationProcessGUI;
        end
        function output = getDrawableOutput()
            output(1).name='Force  field';
            output(1).var='forceField';
            output(1).formatData=@(x) [x.pos x.vec];
            output(1).type='movieOverlay';
            output(1).defaultDisplayMethod=@(x) VectorFieldDisplay('Color','r');
        end
        
    end
end

