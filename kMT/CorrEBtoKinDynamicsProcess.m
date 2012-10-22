classdef CorrEBtoKinDynamicsProcess < DataProcessingProcess
    % A class for correlating kEB signal with kinetochore dynamics
    %
    % Khuloud Jaqaman, October 2012
    
    methods (Access = public)
        
        function obj = CorrEBtoKinDynamicsProcess(owner, varargin)
            
            if nargin == 0
                super_args = {};
            else
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) (isa(x,'MovieData')||isa(x,'MovieList')));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                super_args{1} = owner;
                super_args{2} = CorrEBtoKinDynamicsProcess.getName;
                super_args{3} = @corrEBtoKinDynamics;
                if isempty(funParams)  % Default funParams
                    funParams = CorrEBtoKinDynamicsProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            obj = obj@DataProcessingProcess(super_args{:});
            
        end
        
        function varargout = loadChannelOutput(obj,iChan,varargin)
            
            % Input check
            outputList = {'measurementsEB'};
            ip =inputParser;
            ip.addRequired('iChan',@(x) isscalar(x) && obj.checkChanNum(x));
            %             ip.addOptional('iFrame',1:obj.owner_.nFrames_,@(x) all(obj.checkFrameNum(x)));
            ip.addParamValue('output',outputList,@(x) all(ismember(x,outputList)));
            ip.parse(iChan,varargin{:})
            %             iFrame = ip.Results.iFrame;
            output = ip.Results.output;
            if ischar(output),output={output}; end
            
            % Data loading
            s = load(obj.outFilePaths_{1,iChan},output{:});
            
            for i=1:numel(output),varargout{i}=s.(output{i}); end
            
        end
        
        %         function output = getDrawableOutput(obj)
        %
        %             colors = hsv(numel(obj.owner_.channels_));
        %             output(1).name='Sister pairs';
        %             output(1).var='sisterList';
        %             output(1).formatData=[];
        %             output(1).type='overlay';
        %             output(1).defaultDisplayMethod=@(x)PairsDisplay('Color',colors(x,:));
        %
        %         end
        
    end
    
    methods (Static)
        
        function name = getName()
            name = 'kEB - kin dynamics correlation';
        end
        
        function funParams = getDefaultParams(owner,varargin)
            
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) (isa(x,'MovieData')||isa(x,'MovieList')));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex = 1; % : numel(owner.channels_);
            funParams.OutputDirectory = [outputDir  filesep 'EBtoKinDynamicsCorr'];
            funParams.minDisp = 0.5;
            
        end
        
    end
    
end
