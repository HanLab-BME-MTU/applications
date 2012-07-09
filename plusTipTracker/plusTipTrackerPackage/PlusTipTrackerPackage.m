classdef PlusTipTrackerPackage < TrackingPackage
    % A concrete package for tracking microtubules
    
    methods (Access = public)
        function obj = PlusTipTrackerPackage (varargin)

            % Call the superclass constructor
            obj = obj@TrackingPackage(varargin{:});        
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