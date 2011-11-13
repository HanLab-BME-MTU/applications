classdef FlowAnalysisProcess < DataProcessingProcess
    % Concrete class for a process analyzing flow
    %
    % Sebastien Besson, 5/2011
    
    methods
        function obj = FlowAnalysisProcess(owner,varargin)
            
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
                super_args{2} = FlowAnalysisProcess.getName;
                super_args{3} = @analyzeMovieFlow;
                if isempty(funParams)
                    funParams=FlowAnalysisProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            obj = obj@DataProcessingProcess(super_args{:});
            
        end
        function sanityCheck(obj)
            sanityCheck@DataProcessingProcess(obj)
        end
        
        function varargout = loadChannelOutput(obj,iChan,varargin)
            
            % Input check
            outputList = {'speedMap','Md','Mv','Ms','E','S','img3C_map','img3C_SNR'};
            ip =inputParser;
            ip.addRequired('iChan',@(x) isscalar(x) && ...
                ismember(x,1:numel(obj.owner_.channels_)));
            ip.addOptional('iFrame',1:obj.owner_.nFrames_,...
                @(x) ismember(x,1:obj.owner_.nFrames_));
            ip.addParamValue('output',outputList{1},@(x) all(ismember(x,outputList)));
            ip.parse(iChan,varargin{:})
            iFrame = ip.Results.iFrame;
            
            % Data loading
            output = ip.Results.output;
            if ischar(output), output = {output}; end
            
            % Read file name
            outFileNames = arrayfun(@(x) x.name,...
                dir([obj.outFilePaths_{1,iChan} filesep '*.mat']),'Unif',false);
            for j=1:numel(output)
                varargout{j} = cell(size(iFrame));
            end
            
            % Load output
            for i=1:numel(iFrame)
                kineticMapFile= [obj.outFilePaths_{1,iChan}...
                    filesep outFileNames{iFrame(i)}(1:end-4) '.mat'];
                s = load(kineticMapFile,output{:});
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
            colors = hsv(numel(obj.owner_.channels_));
            output(1).name='Speed maps';
            output(1).var='speedMap';
            output(1).formatData=[];
            output(1).type='image';
            output(1).defaultDisplayMethod=@(x)ImageDisplay('Colormap','jet',...
                'Colorbar','on','Units','nm/min','CLim',obj.getSpeedLimits(x));
            output(2).name='Interpolated vectors';
            output(2).var='Md';
            output(2).formatData=@(x)[x(:,[2 1]) x(:,[4 3])-x(:,[2 1])];
            output(2).type='overlay';
            output(2).defaultDisplayMethod=@(x) VectorFieldDisplay('Color',colors(x,:));
            output(3).name='Noise vectors';
            output(3).var='Ms';
            output(3).formatData=@(x)[x(:,[2 1]) x(:,[4 3])-x(:,[2 1])];
            output(3).type='overlay';
            output(3).defaultDisplayMethod=@(x) VectorFieldDisplay('Color',colors(x,:));
            output(4).name='Error circles';
            output(4).var='E';
            output(4).formatData=[];
            output(4).type='overlay';
            output(4).defaultDisplayMethod=@(x) CircleDisplay('Color',colors(x,:));
            output(5).name='SNR circles';
            output(5).var='S';
            output(5).formatData=[];
            output(5).type='overlay';
            output(5).defaultDisplayMethod=@(x) CircleDisplay('Color',colors(x,:));
            output(6).name='Error maps';
            output(6).var='img3C_map';
            output(6).formatData=[];
            output(6).type='image';
            output(6).defaultDisplayMethod=@ImageDisplay;
            output(7).name='SNR maps';
            output(7).var='img3C_SNR';
            output(7).formatData=[];
            output(7).type='image';
            output(7).defaultDisplayMethod=@ImageDisplay;
        end
    end
    methods (Access=protected)
        function output = getSpeedLimits(obj,iChan)
            speedMap=loadChannelOutput(obj,iChan,'output','speedMap');
            allMaps = vertcat(speedMap{:});
            output=[min(allMaps(:)) max(allMaps(:))];
        end        
    end
    
    methods (Static)
        function name =getName()
            name = 'Flow Analysis';
        end
        function h = GUI()
            h= @flowAnalysisProcessGUI;
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
            funParams.OutputDirectory = [outputDir  filesep 'flowAnalysis'];
            funParams.timeWindow = 1;
            funParams.corrLength = 33;
            funParams.gridSize = 11;
            funParams.interpolate = 1;
            funParams.noise = 1;
            funParams.error = 1;
        end        
    end
end