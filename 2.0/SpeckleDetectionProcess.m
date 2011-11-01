classdef SpeckleDetectionProcess < ImageAnalysisProcess
    % Concrete class for the speckle detection process
    %
    % Sebastien Besson, May 2011
    
    methods
        function obj = SpeckleDetectionProcess(owner,outputDir, funParams)
            
            if nargin == 0
                super_args = {};
            else
                super_args{1} = owner;
                super_args{2} = SpeckleDetectionProcess.getName;
                super_args{3} = @detectMovieSpeckles;
                if nargin < 3 || isempty(funParams)
                    
                    %----Defaults----%
                    funParams.ChannelIndex = 1 : numel(owner.channels_);
                    funParams.MaskChannelIndex = 1:numel(owner.channels_);
                    funParams.OutputDirectory = [outputDir  filesep 'speckles'];
                    funParams.paramSpeckles = [1 0];                
                    funParams.alpha = .01;
                    funParams.I0 = [];
                    funParams.sDN = [];
                    funParams.GaussRatio = [];
                    psfSigmaCheck =arrayfun(@(x)isempty(x.psfSigma_),owner.channels_);
                    if any(psfSigmaCheck)
                        funParams.filterSigma = zeros(size(owner.channels_));
                    else
                        funParams.filterSigma = [owner.channels_.psfSigma_];
                    end
                end
                super_args{4} = funParams;
            end
            
            obj = obj@ImageAnalysisProcess(super_args{:});
        end

            
        function varargout = loadChannelOutput(obj,iChan,varargin)
            
            % Input check
            outputList = {'cands','locMax'};
            ip =inputParser;
            ip.addRequired('obj',@(x) isa(x,'SpeckleDetectionProcess'));
            ip.addRequired('iChan',@(x) ismember(x,1:numel(obj.owner_.channels_)));
            ip.addOptional('iFrame',1:obj.owner_.nFrames_,...
                @(x) ismember(x,1:obj.owner_.nFrames_));
            ip.addParamValue('output','cands',@(x) all(ismember(x,outputList)));
            ip.parse(obj,iChan,varargin{:})
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
            markers={'o','^','s','d'};
            output(1).name='Speckle candidates';
            output(1).var='cands';
            output(1).formatData=@(x) x([x.status]==1);
            output(1).type='overlay';
            output(1).defaultDisplayMethod=@(x) SpeckleDisplay('Marker',markers{x});
        end
        
    end
    methods (Static)
        function name =getName()
            name = 'Speckle Detection';
        end
        function h = GUI()
            h= @speckleDetectionProcessGUI;
        end
        
    end
        
end

