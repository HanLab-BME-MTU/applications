classdef FocalAdhesionPackage < Package
    % The main class of the focal adhesion package
    %
    % NOTE: This is the package that is based on anisotropic spot detection,
    %   tracking and track grouping.
    %   This package is therefore better suited to detection of small
    %   (though not truly nascent, as the detection assumes anisotropy) and
    %   primarily fibrilar (linear) adhesions. It does not support sampling
    %   of image intensities within focal adhesions. 
    %   **For sampling of image intensities within adhesions, and for
    %   support which includes larger non-fibrilar adhesions (but still
    %   misses very small/dim nascent adhesions) see
    %   FocalAdhesionSegmentationPackage 
    %               - Hunter
    % Sebastien Besson, May 2011
    % Updated by Andrew R. Jamieson Feb 2017 
    
    methods
        function obj = FocalAdhesionPackage(owner, varargin)
            % Constructor of class FocalAdhesionPackage
            
            if nargin == 0
                super_args = {};
            else
                % Check input
                ip =inputParser;
                ip.addRequired('owner', @(x) isa(x,'MovieData'));
                ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
                ip.parse(owner, varargin{:});
                outputDir = ip.Results.outputDir;

                % Owner: MovieData object
                super_args{1} = owner;
                super_args{2} = [outputDir  filesep 'FocalAdhesionPackage'];
            end
            
            % Call the superclass constructor
            obj = obj@Package(super_args{:});
        end
        
        function [status processExceptions] = sanityCheck(obj, varargin) % throws Exception Cell Array
        %     nProcesses = length(obj.getProcessClassNames);
            
        %     ip = inputParser;
        %     ip.CaseSensitive = false;
        %     ip.addRequired('obj');
        %     ip.addOptional('full',true, @(x) islogical(x));
        %     ip.addOptional('procID',1:nProcesses,@(x) (isvector(x) && ~any(x>nProcesses)) || strcmp(x,'all'));
        %     ip.parse(obj,varargin{:});
        %     full = ip.Results.full;
        %     procID = ip.Results.procID;
        %     if strcmp(procID,'all'), procID = 1:nProcesses;end
            
            %% TODO - Add Additional appropriate checks
            [status, processExceptions] = sanityCheck@Package(obj, varargin{:});
            
        %     if ~full, return; end
            
        %     validProc = procID(~cellfun(@isempty,obj.processes_(procID)));
        %     maskProc = obj.processes_{2};
        %     detProc = obj.processes_{3};
        %     trackProc = obj.processes_{4};
        %     groupProc = obj.processes_{5};
            
        %     % Set the MaskProcessIndex of the detection
        %     if all(ismember([2 3], validProc))
        %         maskProcIndex = maskProc.getIndex();
        %         funParams.MaskProcessIndex = maskProcIndex;
        %         funParams.MaskChannelIndex = maskProc.funParams_.ChannelIndex;
        %         parseProcessParams(detProc, funParams);
        %     end
            
        %     % Set detection process index for tracking and track grouping
        %     if ~isempty(detProc),
        %         detProcIndex = detProc.getIndex();
        %         funParams.DetProcessIndex = detProcIndex;
        %         if ismember(4,validProc)
        %             parseProcessParams(trackProc, funParams);
        %         end  
        %         if ismember(5,validProc)
        %             parseProcessParams(groupProc, funParams);
        %         end  
        %     end
            
        %     % Set the tracking process index of the track grouping
        %     if all(ismember([4 5], validProc))
        %         trackProcIndex = trackProc.getIndex();
        %         funParams.TrackProcessIndex = trackProcIndex;
        %         parseProcessParams(groupProc, funParams);
        %     end
        % end
        
    end
    
    methods (Static)
        
        function m = getDependencyMatrix(i,j)
            %    1 2 3 4 5 6          
            m = [0 0 0 0 0 0;  %1 Thresholding [optional]
                 2 0 0 0 0 0;  %2 Mask Refinement [optional]
                 0 2 0 0 0 0;  %3 PointSourceProcess
                 0 0 1 0 0 0;  %4 TrackingProcess
                 0 0 1 1 0 0;  %5 FocalAdhesionSegmentationProcess
                 0 0 0 1 1 0;];%6 AnalyzeAdhesionMaturationProcess

            if nargin<2, j=1:size(m,2); end
            if nargin<1, i=1:size(m,1); end
            m=m(i,j);
        end
        
        function name = getName()
            name='New Focal Adhesion';
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
                @FocalAdhesionSegmentationProcess, ...
                @TrackGroupingProcess};
            
            if nargin==0, index=1:numel(integratorProcConstr); end
            procConstr=integratorProcConstr(index);
        end
        function classes = getProcessClassNames(index)
            procContrs = {
                'ThresholdProcess',...
                'MaskRefinementProcess',...
                'AnisoGaussianDetectionProcess',... % Make a generic detectionprocess?
                'TrackingProcess',...
                'FocalAdhesionSegmentationProcess'
                'AdhesionAnalysisProcess'};
            if nargin==0, index=1:numel(procContrs); end
            classes=procContrs(index);
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
