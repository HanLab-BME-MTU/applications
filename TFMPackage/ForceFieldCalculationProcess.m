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
            
            if nargin < 3 || isempty(funParams)
                obj.funParams_=defaultParams;
            else
                obj.funParams_=parseProcessParams(defaultParams,funParams);
            end

        end
        function sanityCheck(obj)
            
        end
        
        function status = checkChannelOutput(obj,varargin)
            
           %Checks if the selected channels have valid output files
           ip =inputParser;
           ip.addRequired('obj',@(x) isa(x,'KineticAnalysisProcess'));
           ip.addOptional('iChan',1:numel(obj.owner_.channels_),...
               @(x) ismember(x,1:numel(obj.owner_.channels_)));
           ip.parse(obj,varargin{:});
           iChan=ip.Results.iChan;
           
           %Makes sure there's at least one .mat file in the speified
           %directory
           status = all(arrayfun(@(x) exist(obj.outFilePaths_{x},'file'),...
               iChan));           
        end
        
        function varargout = loadChannelOutput(obj,iChan,varargin)
            
            outputList = {'speckleArray','kinScore','polyMap',...
                'depolyMap','kinMap2C'};
            ip =inputParser;
            ip.addRequired('obj',@(x) isa(x,'KineticAnalysisProcess'));
            ip.addRequired('iChan',@(x) isscalar(x) && ...
                ismember(x,1:numel(obj.owner_.channels_)));
            ip.addOptional('iFrame',1:obj.owner_.nFrames_,...
                @(x) ismember(x,1:obj.owner_.nFrames_));
            ip.addParamValue('output',outputList{1},@(x) all(ismember(x,outputList)));
            ip.parse(obj,iChan,varargin{:})
            iFrame = ip.Results.iFrame;
                                  
            % Data loading
            output = ip.Results.output;
            if ischar(output), output = {output}; end
            s = load(obj.outFilePaths_{1,iChan},output{:});
            
            if numel(iFrame)>1,
                for i=1:numel(output),
                    varargout{i}=s.(output{i});
                end
            else
                for i=1:numel(output),
                    varargout{i}=s.(output{i}){iFrame};
                end
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
        
    end
end

