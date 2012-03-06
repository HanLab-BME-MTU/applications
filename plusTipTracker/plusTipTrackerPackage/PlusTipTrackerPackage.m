classdef PlusTipTrackerPackage < TrackingPackage
    % A concrete package for tracking microtubules (with pre-set constructors)
    
    methods (Access = public)
        function obj = PlusTipTrackerPackage (owner,varargin)
            if nargin == 0
                super_args = {};
            else
                % Check input
                ip =inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieObject'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                
                super_args{1} = owner;
                super_args{2} = [outputDir  filesep 'PlusTipTrackerPackage'];
            end
            % Call the superclass constructor
            obj = obj@TrackingPackage(super_args{:});        
        end
    end
    methods (Static)
        
        function name = getName()
            name = 'plusTipTracker';
        end
        
        function varargout = GUI(varargin)
            % Start the package GUI
            varargout{1} = plusTipTrackerPackageGUI(varargin{:});
        end
        
        function classes = getProcessClassNames(index)
            classes = TrackingPackage.getProcessClassNames;
            classes{3}='CometPostTrackingProcess';
            
            if nargin==0, index=1:numel(classes); end
            classes=classes(index);
        end
        
        
        function procConstr = getDefaultProcessConstructors(index)
            procConstr = {
                @CometDetectionProcess,...
                @(x,y)TrackingProcess(x,y,PlusTipTrackerPackage.getDefaultTrackingParams(x,y)),...
                @CometPostTrackingProcess};
            if nargin==0, index=1:numel(procConstr); end
            procConstr=procConstr(index);
        end
        
        function funParams = getDefaultTrackingParams(owner,outputDir)
            funParams = TrackingProcess.getDefaultParams(owner,outputDir);
            % Set default minimum track length
            funParams.gapCloseParam.minTrackLen = 3;
            % Set default kalman functions
            kalmanFunctions = TrackingProcess.getKalmanFunctions(2);
            fields = fieldnames(kalmanFunctions);
            validFields = {'reserveMem','initialize','calcGain','timeReverse'};
            kalmanFunctions = rmfield(kalmanFunctions,fields(~ismember(fields,validFields)));
            funParams.kalmanFunctions = kalmanFunctions;            
            % Set default cost matrices
            funParams.costMatrices(1) = TrackingProcess.getDefaultLinkingCostMatrices(owner, funParams.gapCloseParam.timeWindow,2);
            funParams.costMatrices(2) = TrackingProcess.getDefaultGapClosingCostMatrices(owner, funParams.gapCloseParam.timeWindow,2);
        end
        
    end
    
end