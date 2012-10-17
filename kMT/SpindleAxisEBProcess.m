classdef SpindleAxisEBProcess < DataProcessingProcess
    % A concrete class associated to the estimation of spindle axis from EB images
    %
    % Khuloud Jaqaman, October 2012
    
    methods (Access = public)
        
        function obj = SpindleAxisEBProcess(owner, varargin)
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
                
                super_args{1} = owner;
                super_args{2} = SpindleAxisEBProcess.getName;
                super_args{3} = @getMovieSpindleAxisEB;
                if isempty(funParams)  % Default funParams
                    funParams = SpindleAxisEBProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            obj = obj@DataProcessingProcess(super_args{:});
        end
        
        function varargout = loadChannelOutput(obj,iChan,varargin)
            
            % Input check
            outputList = {'spindleAxisVec','poleInfo'};
            ip =inputParser;
            ip.addRequired('iChan',@(x) isscalar(x) && obj.checkChanNum(x));
            ip.addOptional('iFrame',1:obj.owner_.nFrames_,@(x) all(obj.checkFrameNum(x)));
            ip.addParamValue('output',outputList,@(x) all(ismember(x,outputList)));
            ip.parse(iChan,varargin{:})
            iFrame = ip.Results.iFrame;
            output = ip.Results.output;
            if ischar(output),output={output}; end
            
            % Data loading
            s = load(obj.outFilePaths_{1,iChan},output{:});
            if ~isempty(iFrame) && isfield(s,'spindleAxisVec'),
                s.spindleAxisVec = s.spindleAxisVec(iFrame,:);
            end
            if ~isempty(iFrame) && isfield(s,'poleInfo'),
                s.poleInfo = s.poleInfo(iFrame);
            end
            
            for i=1:numel(output),varargout{i}=s.(output{i}); end
            
        end
        
        function output = getDrawableOutput(obj)
            
            colors = hsv(numel(obj.owner_.channels_));
            output(1).name='Spindle Poles';
            output(1).var='poleInfo';
            output(1).formatData=@SpindleAxisEBProcess.formatOutput;
            output(1).type='overlay';
            output(1).defaultDisplayMethod=@(x) LineDisplay('Marker','o',...
                'LineStyle','none','Color',colors(x,:));
            
        end
        
    end
    
    methods (Static)
        
        function name = getName()
            name = 'Spindle Axis';
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
            funParams.OutputDirectory = [outputDir  filesep 'SpindleAxisEB'];
            funParams.doPlot = 1;
            
        end
        
        function y = formatOutput(x)
            
            % Format output in xy coordinate system
            if isempty(x.xCoord)
                y = NaN(1,2);
            else
                y = horzcat(x.xCoord(:,1),x.yCoord(:,1));
            end
            
        end
        
    end
    
end
