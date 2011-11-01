classdef KineticAnalysisProcess < Process
    % Concrete class for a kinetic analysis process
    %
    % Sebastien Besson, 5/2011
    
    methods
        function obj = KineticAnalysisProcess(owner,outputDir,funParams)
            
            if nargin == 0
                super_args = {};
            else
                super_args{1} = owner;
                super_args{2} = KineticAnalysisProcess.getName;
            end
            
            obj = obj@Process(super_args{:});
            
            obj.funName_ = @analyzeMovieSpeckles;
            if nargin < 3 || isempty(funParams)
                
                %----Defaults----%
                funParams.ChannelIndex = 1 : numel(owner.channels_);
                funParams.OutputDirectory = [outputDir  filesep 'kineticAnalysis'];
                funParams.bleachRed = 0;
                funParams.timeWindow = 1;
                funParams.sigma = 5;
            end
            %Make sure the input parameters are legit??
            obj.funParams_ = funParams;
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
            
            outputList = {'kinScore','polyMap',...
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
    end
    methods (Static)
        function name =getName()
            name = 'Kinetic Analysis';
        end
        function h = GUI()
            h= @kineticAnalysisProcessGUI;
        end

        function output =getDrawableOutput()
            output(1).name='Combined map';
            output(1).formatData=[];
            output(1).var='kinMap2C';
            output(1).type='image';
            output(1).defaultDisplayMethod=@ImageDisplay;
			output(2).name='Polymerization map';
            output(2).var='polyMap';
            output(2).formatData=@mat2gray;
            output(2).type='image';
            output(2).defaultDisplayMethod=@ImageDisplay;
			output(3).name='Depolymerization map';
            output(3).formatData=@(x)mat2gray(-x);
            output(3).var='depolyMap';
            output(3).type='image';
            output(3).defaultDisplayMethod=@ImageDisplay;
        end
    end
end

