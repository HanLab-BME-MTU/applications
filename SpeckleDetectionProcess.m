classdef SpeckleDetectionProcess < ImageAnalysisProcess
    % Concrete class for the speckle detection process
    %
    % Sebastien Besson, May 2011
    
    methods
        function obj = SpeckleDetectionProcess(owner,varargin)
            
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
                super_args{2} = SpeckleDetectionProcess.getName;
                super_args{3} = @detectMovieSpeckles;
                if isempty(funParams)
                    funParams = SpeckleDetectionProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            obj = obj@ImageAnalysisProcess(super_args{:});
        end
        
        
        function varargout = loadChannelOutput(obj,iChan,varargin)
            
            % Input check
            outputList = {'cands','locMax'};
            ip =inputParser;
            ip.addRequired('iChan',@(x) isscalar(x) && obj.checkChanNum(x));
            ip.addOptional('iFrame',1:obj.owner_.nFrames_,@(x) all(obj.checkFrameNum(x)));
            ip.addParamValue('output',outputList,@(x) all(ismember(x,outputList)));
            ip.parse(iChan,varargin{:})
            iFrame = ip.Results.iFrame;
            output = ip.Results.output;
            if ischar(output),output={output}; end
            
            % Data loading
            outFileNames = arrayfun(@(x) x.name,...
                dir([obj.outFilePaths_{1,iChan} filesep '*.mat']),'Unif',false);
            for j=1:numel(output)
                varargout{j} = cell(size(iFrame));
            end
            
            for i=1:numel(iFrame)
                inSpeckle= [obj.outFilePaths_{1,iChan}...
                    filesep outFileNames{iFrame(i)}(1:end-4) '.mat'];
                s = load(inSpeckle,output{:});
                for j=1:numel(output)
                    varargout{j}{i} = s.(output{j});
                end
            end
            if numel(iFrame)==1,
                for j=1:numel(output)
                    varargout{j} = varargout{j}{1};
                end
            end
        end
        function output = getDrawableOutput(obj)

            output(1).name='Speckle candidates';
            output(1).var='cands';
            output(1).formatData=@(x) x([x.status]==1);
            output(1).type='overlay';
            colors = hsv(numel(obj.owner_.channels_));
            output(1).defaultDisplayMethod=@(x) SpeckleDisplay('Color',colors(x,:));
        end
        
    end
    methods (Static)
        function name =getName()
            name = 'Speckle Detection';
        end
        function h = GUI()
            h= @speckleDetectionProcessGUI;
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Define default process parameters
            funParams.ChannelIndex = 1 : numel(owner.channels_);
            funParams.MaskChannelIndex = 1:numel(owner.channels_);
            funParams.OutputDirectory = [outputDir  filesep 'speckles'];
            funParams.paramSpeckles = [1 0];
            funParams.alpha = .01;
            funParams.I0 = [];
            funParams.sDN = [];
            funParams.GaussRatio = [];
            psfSigmaCheck =~arrayfun(@(x)isempty(x.psfSigma_),owner.channels_);
            if all(psfSigmaCheck)
                funParams.filterSigma = [owner.channels_.psfSigma_];
            else
                funParams.filterSigma = zeros(size(owner.channels_));
            end
        end
        
    end
end