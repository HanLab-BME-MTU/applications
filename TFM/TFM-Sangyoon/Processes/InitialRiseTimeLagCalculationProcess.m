classdef InitialRiseTimeLagCalculationProcess < DataProcessingProcess
    methods (Access = public)
        function obj = InitialRiseTimeLagCalculationProcess(owner,varargin)
%             obj = obj@DataProcessingProcess(owner, InitialRiseTimeLagCalculationProcess.getName);
%             obj.funName_ = @calculateInitialRiseTimeLagFromTracks; % This should be variation from colocalizationAdhesionWithTFM
%             obj.funParams_ = InitialRiseTimeLagCalculationProcess.getDefaultParams(owner,varargin{1});
            if nargin == 0
                super_args = {};
            else
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieData'));
                ip.addOptional('outputDir', owner.outputDirectory_,@ischar);
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                % Define arguments for superclass constructor
                super_args{1} = owner;
                super_args{2} = InitialRiseTimeLagCalculationProcess.getName;
                super_args{3} = @calculateInitialRiseTimeLagFromTracks;
                
                if isempty(funParams)
                    funParams = InitialRiseTimeLagCalculationProcess.getDefaultParams(owner,outputDir);
                end
                
                super_args{4} = funParams;
            end
            
            obj = obj@DataProcessingProcess(super_args{:});
        end
        
        function output = loadChannelOutput(obj, iChan, varargin)
            outputList = {};
            nOutput = length(outputList);

            ip.addRequired('iChan',@(x) obj.checkChanNum(x));
            ip.addOptional('iOutput',1,@(x) ismember(x,1:nOutput));
            ip.addParamValue('output','',@(x) all(ismember(x,outputList)));
            ip.addParamValue('useCache',false,@islogical);
            ip.parse(iChan,varargin{:})
    
            s = cached.load(obj.outFilePaths_{iChan},'-useCache',ip.Results.useCache);

            output = s.Imean;          
        end
    end
    methods (Static)
        function name = getName()
            name = 'Initial-Rise Time Lag Calculation';
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner', @(x) isa(x, 'MovieObject'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            adhAnalProc = owner.getProcess(owner.getProcessIndex('AdhesionAnalysisProcess'));
            pAnal=adhAnalProc.funParams_;
            
            ip.addOptional('ChannelIndex',pAnal.ChannelIndex,...
               @(x) all(owner.checkChanNum(x)));
            ip.parse(owner,varargin{:})
            
            % Set default parameters
            funParams.OutputDirectory = [ip.Results.outputDir filesep 'InitialRiseTimeLagCalculation'];
            funParams.ChannelIndex = ip.Results.ChannelIndex;
            funParams.mainSlave = 1;
        end
        
        function h = GUI()
            h = @initialRiseTimeLagCalculationProcessGUI;
        end
    end
end
