classdef StrainEnergyCalculationProcess < DataProcessingProcess
    % Concrete process for calculating a force field
    %
    % Sangyoon Han & Andrew Jamieson, May 2017
    properties (SetAccess = protected)  
        tMapLimits_
        dELimits_
        distBeadMapLimits_
    end
    
    methods
        function obj = StrainEnergyCalculationProcess(owner,varargin)
            
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
                super_args{2} = StrainEnergyCalculationProcess.getName;
                super_args{3} = @calculateMovieStrainEnergy;
                if isempty(funParams)
                    funParams=StrainEnergyCalculationProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            obj = obj@DataProcessingProcess(super_args{:});
        end
        
        function varargout = loadChannelOutput(obj, iFrame, varargin)
            % Input check
            outputList = {'SE_Blobs','SE_FOV','SE_Cell'};
            ip = inputParser;
            ip.addRequired('obj');
            ip.addRequired('iFrame', @(x) obj.checkFrameNum(x));
            ip.addParameter('useCache',true, @islogical);
            ip.addParameter('output', outputList{3}, @(x) all(ismember(x,outputList)));
            ip.parse(obj,iFrame, varargin{:})
            output = ip.Results.output;
            varargout = cell(numel(output), 1);
            iFrame = ip.Results.iFrame;
            if ischar(output),output={output}; end

            % Load all the output data
            s = cached.load(obj.outFilePaths_{1,10}, '-useCache', ip.Results.useCache);
            
            for iout = 1:numel(output)
                
                switch output{iout} 
                    case 'SE_FOV' 
                        val = s.(output{iout}).SE(iFrame);                   
                    case 'SE_Cell'
                        val = s.(output{iout}).SE(iFrame);
                    case 'SE_Blobs'
                        val = s.(output{iout}).SE(iFrame);
                    otherwise
                        error('Incorrect Output Var type');
                end
                varargout{iout} = val;
            end 
        % outputFile{1,1} = [p.OutputDirectory filesep 'strainEnergyInFOV.mat'];
        % outputFile{1,2} = [p.OutputDirectory filesep 'strainEnergyInCell.mat'];
        % outputFile{1,3} = [p.OutputDirectory filesep 'forceBlobs.mat'];

        % @(proc,iChan,iFrame,varargin) proc.getParameters().text;
        % frame = arrayfun(@(x) ['Frame ' num2str(x)], 1:3, 'Uniform' ,0)';
        % t = struct2table(SE_Blobs);
        % t.Properties.RowNames = t.Frame
        % t.Frame = []
        end

        function output = getDrawableOutput(obj)
            i = 1;
            output(i).name = 'SE FOV';
            output(i).var = 'SE_FOV';
            output(i).formatData = @(val) struct('String',{{['SE FOV: ' num2str(val)]}},'Color',[1 0 0],'Position',[10 10]);
            output(i).type = 'movieOverlay';
            output(i).defaultDisplayMethod = @TextDisplay;
            
            i = 2;
            output(i).name = 'SE Cell';
            output(i).var = 'SE_Cell';
            output(i).formatData = @(val) struct('String', {{['SE Cell: ' num2str(val)]}},'Color',[0 0 1],'Position',[10 20]);
            output(i).type = 'movieOverlay';
            output(i).defaultDisplayMethod = @TextDisplay;            
            
            i = 3;
            output(i).name = 'SE Cell';
            output(i).var = 'SE_Blobs';
            output(i).formatData = @(val) struct('String', {{['SE Blobs: ' num2str(val)]}},'Color',[0 0 1],'Position',[10 20]);
            output(i).type = 'movieOverlay';
            output(i).defaultDisplayMethod = @TextDisplay;
            
        end
    end
%         function status = checkChannelOutput(obj,varargin)
%             
%             status = logical(exist(obj.outFilePaths_{1},'file'));
%             
%         end
    methods (Static)
        function name =getName()
            name = 'Strain Energy Calculation';
        end
        function h = GUI()
            h= @strainEnergyCalculationProcessGUI;
        end
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.OutputDirectory = [outputDir  filesep 'strainEnergy'];
            funParams.exportCSV=true;
            funParams.performForceBlobAnalysis=true;
            funParams.useFOV=true;
            funParams.useCellMask=true;
        end
%         function units = getUnits(varargin)
%             units = 'Traction (Pa)';
%         end
    end
end