classdef FocalAdhesionPackage < Package
    % The main class of the focal adhesion package
    %
    % NOTE: This is the package that is NOT based on anisotropic spot detection,
    %   anymore.  - Sangyoon
    % Sebastien Besson, May 2011
    % Updated by Andrew R. Jamieson Feb 2017 
    % Updated by Sangyoon J. Han Oct 2017 
    
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
            % Check that the time interval is correctly setup
            % Check if there is at least one channel have a value for the psf sigma
            
            psfSigmaCheck =arrayfun(@(x) ~isempty(x.psfSigma_),obj.owner_.channels_);
            assert(any(psfSigmaCheck),...
                ['Missing standard deviation of the theoretical point-spread function! '...
                'Please fill the numerical aperture, pixel size and'...
                ' emission wavelengths of all channels!']);

            missingMetadataMsg = ['Missing %s! The %s is necessary to analyze '...
                'Focal Adhesions. Please edit the movie and fill the %s.'];
            errorMsg = @(x) sprintf(missingMetadataMsg, x, x, x);
            
            assert(~isempty(obj.owner_.pixelSize_), errorMsg('pixel size'));
            assert(~isempty(obj.owner_.timeInterval_), errorMsg('time interval'));
            assert(~isempty(obj.owner_.camBitdepth_), errorMsg('camera bit depth'));

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
        end
        
    end
    
    methods (Static)

        function name = getName()
            name = 'New Focal Adhesion Package';
        end

        function varargout = GUI(varargin)
            %% TODO - Update GUI
            varargout{1} = FocalAdhesionPackageGUI(varargin{:});
        end        

        function classes = getProcessClassNames(index)
            procContrs = {
                'EfficientSubpixelRegistrationProcess',...
                'ThresholdProcess',...
                'MaskRefinementProcess',...
                'DetectionProcess',... % default should be pointsourcedetection
                'TrackingProcess',...
                'FocalAdhesionSegmentationProcess',...
                'AdhesionAnalysisProcess',...
                'AdhesionClassificationProcess',...
                'TheOtherChannelReadingProcess',...
                'TractionForceReadingProcess'...
                'InitialRiseTimeLagCalculationProcess'};
            if nargin==0, index=1:numel(procContrs); end
            classes=procContrs(index);
        end        

        function m = getDependencyMatrix(i,j)
            %    1 2 3 4 5 6 7 8 9 10 11 {Processes}    
            m = [0 0 0 0 0 0 0 0 0 0 0;  %1 Stage Drift Correction [optional]
                 0 0 0 0 0 0 0 0 0 0 0;  %2 Thresholding [optional]
                 0 1 0 0 0 0 0 0 0 0 0;  %3 Mask Refinement [optional]
                 0 2 2 0 0 0 0 0 0 0 0;  %4 DetectionProcess
                 0 0 0 1 0 0 0 0 0 0 0;  %5 TrackingProcess
                 0 2 2 0 0 0 0 0 0 0 0;  %6 FocalAdhesionSegmentationProcess
                 2 0 1 1 1 1 0 0 0 0 0;  %7 AnalyzeAdhesionMaturationProcess
                 0 0 2 1 1 1 1 0 0 0 0;  %8 AdhesionClassificationProcess
                 0 0 0 0 0 0 1 0 0 0 0;  %9 TheOtherChannelReadingProcess
                 0 0 0 0 0 0 1 0 0 0 0;  %10 TractionForceReadingProcess
                 0 0 0 0 0 0 1 2 2 2 0]; %11 InitialRiseTimeLagCalculationProcess

            if nargin<2, j=1:size(m,2); end
            if nargin<1, i=1:size(m,1); end
            m=m(i,j);
        end

        function procConstr = getDefaultProcessConstructors(index)
            procContrs = {
                @EfficientSubpixelRegistrationProcess,...
                @ThresholdProcess,...
                @MaskRefinementProcess,...                
                @(x,y)PointSourceDetectionProcess(x,y,FocalAdhesionPackage.getDefaultDetectionParams(x,y)),...
                @(x,y)TrackingProcess(x,y,FocalAdhesionPackage.getDefaultTrackingParams(x,y)),...
                @(x,y)FocalAdhesionSegmentationProcess(x,y,FocalAdhesionPackage.getDefaultFASegParams(x,y)), ...
                @(x,y)AdhesionAnalysisProcess(x,y,FocalAdhesionPackage.getDefaultAnalysisParams(x,y)),...
                @(x,y)AdhesionClassificationProcess(x,y,FocalAdhesionPackage.getDefaultClassificationParams(x,y)), ...
                @(x,y)TheOtherChannelReadingProcess(x,y,FocalAdhesionPackage.getDefaultTheOtherChannelReadingParams(x,y)), ...
                @(x,y)TractionForceReadingProcess(x,y,FocalAdhesionPackage.getDefaultForceReadingParams(x,y)), ...
                @(x,y)InitialRiseTimeLagCalculationProcess(x,y,FocalAdhesionPackage.getDefaultInitialRiseTimeLagCalculationParams(x,y))};
            
            if nargin==0, index=1:numel(procContrs); end
            procConstr=procContrs(index);
        end

        function funParams = getDefaultDetectionParams(owner, outputDir)

            funParams = PointSourceDetectionProcess.getDefaultParams(owner, outputDir);
            
            %% TODO - Verify ideal default settings here.
            % Set default parameters
            funParams.ChannelIndex = 1;

            % Check if segmentation occured.
            %% TODO -- Update 
            iProc = owner.getProcessIndex('MaskRefinementProcess', 'askUser', false);
            if isempty(iProc)
                disp('Note: No Cell Segmentation Mask found');
                funParams.MaskChannelIndex = []; %1:numel(owner.channels_);
                funParams.MaskProcessIndex = [];            
            else
                funParams.MaskProcessIndex = iProc; % Specify Process Index with cell mask output
                funParams.MaskChannelIndex = 1:numel(owner.channels_);
            end

            funParams.OutputDirectory = [outputDir  filesep 'point_sources'];
            funParams.alpha=.05;
            funParams.maskRadius=40;
            funParams.Mode = {'xyAc'};
            funParams.FitMixtures = false;
            funParams.MaxMixtures = 5;
            funParams.RedundancyRadius = .25;
            funParams.UseIntersection = true;            
            funParams.PreFilter = true;
            %list of parameters which can be specified at a per-channel
            %level. If specified as scalar these will  be replicated
            funParams.PerChannelParams = {'alpha','Mode','FitMixtures','MaxMixtures','RedundancyRadius','filterSigma','PreFilter','ConfRadius','WindowSize'};
            
            nChan = numel(owner.channels_);
            funParams.filterSigma = 1.2*ones(1,nChan); %Minimum numerically stable sigma is ~1.2 pixels.
            hasPSFSigma = arrayfun(@(x) ~isempty(x.psfSigma_), owner.channels_);
            funParams.filterSigma(hasPSFSigma) = [owner.channels_(hasPSFSigma).psfSigma_];            
            funParams.filterSigma(funParams.filterSigma<1.2) = 1.2; %Make sure default isn't set to too small.
            
            funParams.ConfRadius = arrayfun(@(x)(2*x),funParams.filterSigma);
            funParams.WindowSize = arrayfun(@(x)(ceil(4*x)),funParams.filterSigma);
            
            funParams = prepPerChannelParams(funParams,nChan);
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

        function funParams = getDefaultFASegParams(owner, outputDir)

            funParams = FocalAdhesionSegmentationProcess.getDefaultParams(owner,outputDir);

            % TODO FIX Sanity checks (be careful of optional processes)
%             try owner.getProcessIndex('TrackingProcess')
%                 trackNAProc = owner.getProcess(owner.getProcessIndex('TrackingProcess'));
%             catch
%                 trackNAProc = -1
%             end
%             detectedNAProc = owner.getProcess(owner.getProcessIndex('DetectionProcess'));
%             
%             %%TODO - Move to sanity check section
%             if length(detectedNAProc.funParams_.ChannelIndex) > 1
%                 numChan = length(trackNAProc.funParams_.ChannelIndex);
%                 for i = 1:numChan
%                     assert(detectedNAProc.funParams_.ChannelIndex(i) == trackNAProc.funParams_.ChannelIndex(i), 'ChannelInex should match');        
%                 end
%             elseif length(detectedNAProc.funParams_.ChannelIndex) == 1
%                 assert(detectedNAProc.funParams_.ChannelIndex == trackNAProc.funParams_.ChannelIndex, 'ChannelInex should match');    
%             end
            
            
            %%TODO - Sangyoon suggests an automated selection of parameters
            funParams.ChannelIndex = 1;%trackNAProc.funParams_.ChannelIndex;
            funParams.SteerableFilterSigma = 72; % %Sigma in nm of steerable filter to use in splitting adjacent adhesions
            funParams.OpeningRadiusXY = 0; % %Spatial radius in nm of structuring element used in opening.
            funParams.MinVolTime = 1; %um2*s Minimum spatiotemporal "Volume" in micron^2 * seconds of segmented adhesions to retain.
            funParams.OpeningHeightT = 10; % Temporal "height" in seconds of structuring element used in opening            
        end
    
        function funParams = getDefaultAnalysisParams(owner, outputDir)

            funParams = AdhesionAnalysisProcess.getDefaultParams(owner, outputDir);
            
            iProcDetect = owner.getProcessIndex('DetectionProcess');
            iProcTrack = owner.getProcessIndex('TrackingProcess');
            trackNAProc = owner.getProcess(iProcTrack);
            iProcFASeg = owner.getProcessIndex('FocalAdhesionSegmentationProcess');
            FASegProc = owner.getProcess(iProcFASeg);
            
            % Specify Channels where adhesions are segmented (% Assume Pax Channel?_)
            %%TODO - Move to sanity check section
            assert(FASegProc.funParams_.ChannelIndex == trackNAProc.funParams_.ChannelIndex, 'ChannelInex should match');
            funParams.ChannelIndex = FASegProc.funParams_.ChannelIndex;
            
            funParams.detectedNAProc = iProcDetect;
            funParams.trackFAProc = iProcTrack;
            funParams.FAsegProc = iProcFASeg;
            
            % Check if segmentation occured.
            % TODO -- Update 
            iProc = owner.getProcessIndex('MaskRefinementProcess', 'askUser', false);
            if isempty(iProc)
                disp('Note: No Cell Segmentation Mask found');
                funParams.ApplyCellSegMask = false;
            else
                funParams.SegCellMaskProc = iProc; % Specify Process Index with cell mask output
            end
            
            funParams.showAllTracks = false;
            funParams.plotEachTrack = false;
                        
            funParams.reTrack = true;  %%TODO - Check with sangyoon
            
            funParams.onlyEdge = false; 
            funParams.matchWithFA = true;  
            funParams.getEdgeRelatedFeatures = true;
            funParams.bandwidthNA = 7;
            funParams.minLifetime = 3; % Checked: changed to 3 from 5 -Sangyoon10/25/2017
            
        end
        %% getDefaultClassificationParams
        function funParams = getDefaultClassificationParams(owner, outputDir)

            funParams = AdhesionClassificationProcess.getDefaultParams(owner, outputDir);
            
            iAdhAnal = owner.getProcessIndex('AdhesionAnalysisProcess');
            adhAnalProc = owner.getProcess(iAdhAnal);
            
            % Specify Channels where adhesions are segmented (% Assume Pax Channel?_)
            %%TODO - Move to sanity check section
            funParams.ChannelIndex = adhAnalProc.funParams_.ChannelIndex;
            
            % Check if segmentation occured.            
            funParams.labeledData=[];
            funParams.useAutomaticallySelectedData = true;
            funParams.manualLabeling=false;
        end

        %% getDefaultTheOtherChannelReading
        function funParams = getDefaultTheOtherChannelReadingParams(owner, outputDir)
            funParams = TheOtherChannelReadingProcess.getDefaultParams(owner, outputDir);
        end
        
        %% getDefaultForceReading
        function funParams = getDefaultForceReadingParams(owner, outputDir)
            funParams = TractionForceReadingProcess.getDefaultParams(owner, outputDir);
        end
        %% getDefaultInitialRiseTimeLagCalculation
        function funParams = getDefaultInitialRiseTimeLagCalculationParams(owner, outputDir)
            funParams = InitialRiseTimeLagCalculationProcess.getDefaultParams(owner, outputDir);
        end
    end
end
