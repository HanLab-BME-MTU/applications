classdef FocalAdhesionPackage < Package
    % The main class of the Integrator package
    
    % Sebastien Besson, May 2011
    
    methods
        function obj = FocalAdhesionPackage(owner,varargin)
            % Constructor of class QFSMPackage
            
            if nargin == 0
                super_args = {};
            else
                % Check input
                ip =inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieData'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;

                % Owner: MovieData object
                super_args{1} = owner;
                super_args{2} = [outputDir  filesep 'FocalAdhesionPackage'];
            end
            
            % Call the superclass constructor
            obj = obj@Package(super_args{:});
        end
        
        function [status processExceptions] = sanityCheck(obj,varargin) % throws Exception Cell Array
            nProcesses = length(obj.getProcessClassNames);
            
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.addRequired('obj');
            ip.addOptional('full',true, @(x) islogical(x));
            ip.addOptional('procID',1:nProcesses,@(x) (isvector(x) && ~any(x>nProcesses)) || strcmp(x,'all'));
            ip.parse(obj,varargin{:});
            full = ip.Results.full;
            procID = ip.Results.procID;
            if strcmp(procID,'all'), procID = 1:nProcesses;end
            
            [status processExceptions] = sanityCheck@Package(obj,full,procID);
            
            if ~full, return; end
            
            validProc = procID(~cellfun(@isempty,obj.processes_(procID)));

            % Set the MaskProcessIndex of the detection
            if ismember(3,validProc) &&  ~isempty(obj.processes_{2})
                maskProcIndex = find(cellfun(@(x)isequal(x, obj.processes_{2}), obj.owner_.processes_));
                assert(numel(maskProcIndex)== 1,'User-defined: More than one identical process exists in movie data''s process list.');
                funParams.MaskProcessIndex = maskProcIndex;
                parseProcessParams(obj.processes_{3},funParams);
            end
            
            % Set detection process index
            if ~isempty(obj.processes_{3})
                detProcIndex = find(cellfun(@(x)isequal(x, obj.processes_{3}), obj.owner_.processes_));
                assert(numel(detProcIndex)== 1,'User-defined: More than one identical process exists in movie data''s process list.');
                funParams.DetProcessIndex = detProcIndex;
                if ismember(4,validProc)
                    parseProcessParams(obj.processes_{4},funParams);
                end  
                if ismember(5,validProc)
                    parseProcessParams(obj.processes_{5},funParams);
                end  
            end
            
            % Set the MaskProcessIndex of the detection
            if ismember(5,validProc) &&  ~isempty(obj.processes_{4})
                trackProcIndex = find(cellfun(@(x)isequal(x, obj.processes_{4}), obj.owner_.processes_));
                assert(numel(trackProcIndex)== 1,'User-defined: More than one identical process exists in movie data''s process list.');
                funParams.TrackProcessIndex = trackProcIndex;
                parseProcessParams(obj.processes_{5},funParams);
            end
        end
        
    end
    
    methods (Static)
        
        function m = getDependencyMatrix(i,j)
            
            m = [0 0 0 0 0;  %1 SignalPreprocessingProcess
                1 0 0 0 0;
                0 1 0 0 0;
                0 0 1 0 0;
                0 0 0 1 0;]; %2 SignalProcessingProcess

            if nargin<2, j=1:size(m,2); end
            if nargin<1, i=1:size(m,1); end
            m=m(i,j);
        end
        
        function name = getName()
            name='Focal Adhesion';
        end
        function varargout = GUI(varargin)
            % Start the package GUI
            varargout{1} = focalAdhesionPackageGUI(varargin{:});
        end
        function procConstr = getDefaultProcessConstructors(index)
            integratorProcConstr = {
                @ThresholdProcess,...
                @MaskRefinementProcess,...                
                @AnisoGaussianDetectionProcess,...
                @(x,y)TrackingProcess(x,y,FocalAdhesionPackage.getDefaultTrackingParams(x,y)),...
                @TrackGroupingProcess};
            
            if nargin==0, index=1:numel(integratorProcConstr); end
            procConstr=integratorProcConstr(index);
        end
        function classes = getProcessClassNames(index)
            integratorClasses = {
                'ThresholdProcess',...
                'MaskRefinementProcess',...
                'AnisoGaussianDetectionProcess',...
                'TrackingProcess',...
                'TrackGroupingProcess'};
            if nargin==0, index=1:numel(integratorClasses); end
            classes=integratorClasses(index);
        end
        function funParams = getDefaultTrackingParams(owner,outputDir)
            funParams = TrackingProcess.getDefaultParams(owner,outputDir);
            
            % Set default gap closing parameters
            funParams.gapCloseParam.timeWindow = 3;
            funParams.gapCloseParam.minTrackLen = 1;
            funParams.gapCloseParam.diagnostics = 0;

            % Set default kalman functions
            funParams.kalmanFunctions = TrackingProcess.getKalmanFunctions(1);           

            % Set default cost matrices
            funParams.costMatrices(1) = TrackingProcess.getDefaultLinkingCostMatrices(owner, funParams.gapCloseParam.timeWindow,1);
            funParams.costMatrices(1).parameters.linearMotion = 2;
            funParams.costMatrices(1).parameters.minSearchRadius = 5;
            funParams.costMatrices(1).parameters.maxSearchRadius = 5;
            funParams.costMatrices(1).parameters.diagnostics = [];
            
            funParams.costMatrices(2) = TrackingProcess.getDefaultGapClosingCostMatrices(owner, funParams.gapCloseParam.timeWindow,1);
            funParams.costMatrices(2).parameters.linearMotion = 2;
            funParams.costMatrices(2).parameters.minSearchRadius = 5;
            funParams.costMatrices(2).parameters.maxSearchRadius = 5;
            funParams.costMatrices(2).parameters.maxAngleVV = 45;

        end
    end
end
