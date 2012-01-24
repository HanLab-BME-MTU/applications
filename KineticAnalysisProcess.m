classdef KineticAnalysisProcess < DataProcessingProcess
    % Concrete class for a kinetic analysis process
    %
    % Sebastien Besson, 5/2011
    
    properties (SetAccess = protected)  
        kineticLimits_
    end
    
    methods
        function obj = KineticAnalysisProcess(owner,varargin)
            
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
                super_args{2} = KineticAnalysisProcess.getName;
                super_args{3} = @analyzeMovieSpeckles;
                if isempty(funParams)
                    funParams=KineticAnalysisProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            obj = obj@DataProcessingProcess(super_args{:});            
        end
        
        function varargout = loadChannelOutput(obj,iChan,varargin)
            
            outputList = {'kinMap2C','polyMap','depolyMap','netMap'};
            ip =inputParser;
            ip.addRequired('iChan',@(x) isscalar(x) && obj.checkChanNum(x));
            ip.addOptional('iFrame',1:obj.owner_.nFrames_,@(x) all(obj.checkFrameNum(x)));
            ip.addParamValue('output',outputList,@(x) all(ismember(x,outputList)));
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
                s = load(kineticMapFile);
                for j=1:numel(output)
                    if strcmp(output{j},'netMap')
                        varargout{j}{i} = s.polyMap+s.depolyMap;
                    else
                        varargout{j}{i} = s.(output{j});
                    end
                end
            end
            if numel(iFrame)==1,
                 for j=1:numel(output)
                    varargout{j} = varargout{j}{1};
                 end
            end
        end
        
        function setKineticLimits(obj,kineticLimits)
            obj.kineticLimits_=kineticLimits;
        end
        
        
        function output =getDrawableOutput(obj)
            output(1).name='Net polymerization map';
            output(1).formatData=[];
            output(1).var='netMap';
            output(1).type='image';
            output(1).defaultDisplayMethod=@(x)ImageDisplay('Colormap',obj.createNetColormap(x),...
                'Colorbar','on','CLim',obj.kineticLimits_{x},'Units','Kinetic score');
            
			output(2).name='Polymerization map';
            output(2).var='polyMap';
            output(2).formatData=[];
            output(2).type='image';
            output(2).defaultDisplayMethod=@(x)ImageDisplay('Colormap',(0:1/64:1)'*[1 0 0],...
            'Colorbar','on','CLim',[0 obj.kineticLimits_{x}(2)],'Units','Kinetic score');
            output(3).name='Depolymerization map';
            output(3).formatData=[];
            output(3).var='depolyMap';
            output(3).type='image';
            output(3).defaultDisplayMethod=@(x)ImageDisplay('Colormap',(1:-1/64:0)'*[0 1 0],...
                'Colorbar','on','CLim',[obj.kineticLimits_{x}(1) 0],'Units','Kinetic score');
%             output(4).name='Combined map';
%             output(4).formatData=[];
%             output(4).var='kinMap2C';
%             output(4).type='image';
%             output(4).defaultDisplayMethod=@(x)ImageDisplay('Colormap',...
%                 [(1:-1/32:0)'*[0 1 0]; (0:1/32:1)'*[1 0 0]],'Colorbar','on');
        end
        
    end
    
    methods (Access = protected)
        function cMap =createNetColormap(obj,x)
            % Get scores and normalize them
            minScore= obj.kineticLimits_{x}(1);
            maxScore = obj.kineticLimits_{x}(2);        
            dScore=maxScore-minScore;
            nBins = 256;
            nr= round(nBins*maxScore/dScore);
            ng= nBins-nr;
            g=[0 1 0];
            r=[1 0 0];
            cMap= vertcat((1:-1/ng:0)'*g,(0:1/nr:1)'*r) ;
        end
    end
    
    methods (Static)
        function name =getName()
            name = 'Kinetic Analysis';
        end
        function h = GUI()
            h= @kineticAnalysisProcessGUI;
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
            funParams.OutputDirectory = [outputDir  filesep 'kineticAnalysis'];
            funParams.bleachRed = 0;
            funParams.timeWindow = 5;
            funParams.sigma = 5;
        end        
    end
end
