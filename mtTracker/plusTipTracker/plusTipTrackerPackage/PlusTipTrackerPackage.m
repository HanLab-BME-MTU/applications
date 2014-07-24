classdef PlusTipTrackerPackage < TrackingPackage
    % A concrete package for tracking microtubules
    
    methods (Access = public)
        function obj = PlusTipTrackerPackage (varargin)
            
            % Call the superclass constructor
            obj = obj@TrackingPackage(varargin{:});
        end
        
        function [status, processExceptions] = sanityCheck(obj,varargin)
            
            % Check that the time interval is correctly setup
            missingMetadataMsg = ['Missing %s! The %s is necessary to analyze '...
                'microtubule plus-end tracks. Please edit the movie and fill the %s.'];
            errorMsg = @(x) sprintf(missingMetadataMsg, x, x, x);
            
            assert(~isempty(obj.owner_.pixelSize_), errorMsg('pixel size'));
            assert(~isempty(obj.owner_.timeInterval_), errorMsg('time interval'));
            assert(~isempty(obj.owner_.camBitdepth_), errorMsg('camera bit depth'));
            
            [status, processExceptions] = sanityCheck@Package(obj, varargin{:});
        end
    end
    methods (Static)
        
        function varargout = GUI(varargin)
            % Start the package GUI
            varargout{1} = plusTipTrackerPackageGUI(varargin{:});
        end
        
        function procConstr = getDefaultProcessConstructors(index)
            procConstr = {
                @CometDetectionProcess,...
                @(x,y)TrackingProcess(x,y,PlusTipTrackerPackage.getDefaultTrackingParams(x,y)),...
                @CometPostTrackingProcess};
            if nargin==0, index=1:numel(procConstr); end
            procConstr=procConstr(index);
        end
        
        function procConstr = getAlternateDetector()
            procConstr = @AnisoGaussianDetectionProcess;
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
        
        function tools = getTools(index)
            plusTipTools(1).name = '+TIP Group Analysis';
            plusTipTools(1).funHandle = @plusTipGroupAnalysisGUI;
            if nargin==0, index=1:numel(plusTipTools); end
            tools = plusTipTools(index);
        end
    end
end