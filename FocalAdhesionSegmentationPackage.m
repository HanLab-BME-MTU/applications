classdef FocalAdhesionSegmentationPackage < Package
    % The main class of the focal adhesion segmentation package This is the
    % package that is based on segmentation of adhesions. This package
    % handles both small and large, globular adhesions, and allows sampling
    % of image intensities (and processed images) within segmented focal
    % adhesions. For detection and tracking of small dim adhesions see
    % FocalAdhesionPackage
    
    % Hunter Elliott, April 2013
    
    methods
        function obj = FocalAdhesionSegmentationPackage(owner,varargin)            
            
            if nargin == 0
                super_args = {};
            else
                % Check input
                ip =inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieObject'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;

                % Owner: MovieData object
                super_args{1} = owner;
                super_args{2} = [outputDir  filesep 'FocalAdhesionSegmentationPackage'];
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
            
            [status, processExceptions] = sanityCheck@Package(obj,full,procID);
            
            if ~full, return; end
            
            validProc = procID(~cellfun(@isempty,obj.processes_(procID)));

            %I guess this should not be hard-coded... oh well
            dcProc = obj.processes_{1};
            scProc = obj.processes_{2};            
            threshProc = obj.processes_{3};
            bmProc = obj.processes_{4};
            bsProc = obj.processes_{5};
            asProc = obj.processes_{6};
            ssProc = obj.processes_{7};
               
            %If corrections have been performed, use these images to
            %threshold
            if ~isempty(threshProc) && ~isempty(scProc)                
                thParams.ProcessIndex = scProc.getIndex;
                parseProcessParams(threshProc,thParams);
            end
            
            
            %Set process index for background mask creation so it always
            %uses the thresholding.
            if ~isempty(threshProc) && ~isempty(bmProc)
               threshInd = threshProc.getIndex;
               bmParams.SegProcessIndex = threshInd;
               parseProcessParams(bmProc,bmParams);
            end                
 
% Removed - This is now handled by the GUI.
%             %If corrections are performed, make sure the seg sampling uses
%             %corrected images
%             if ~isempty(bsProc) && ~isempty(ssProc)
%                 ssParams.ProcessIndex = bsProc.getIndex;
%                 parseProcessParams(ssProc,ssParams);
%             end
            
            %Make sure the adhesion uses the adhesion masks 
            if ~isempty(asProc) && ~isempty(ssProc)
                ssParams.SegProcessIndex = asProc.getIndex;
                parseProcessParams(ssProc,ssParams);
            end
                
                
            
%             maskProc = obj.processes_{2};
%             detProc = obj.processes_{3};
%             trackProc = obj.processes_{4};
%             groupProc = obj.processes_{5};            
%             % Set the MaskProcessIndex of the detection
%             if all(ismember([2 3], validProc))
%                 maskProcIndex = maskProc.getIndex();
%                 funParams.MaskProcessIndex = maskProcIndex;
%                 funParams.MaskChannelIndex = maskProc.funParams_.ChannelIndex;
%                 parseProcessParams(detProc, funParams);
%             end
%             
%             % Set detection process index for tracking and track grouping
%             if ~isempty(detProc),
%                 detProcIndex = detProc.getIndex();
%                 funParams.DetProcessIndex = detProcIndex;
%                 if ismember(4,validProc)
%                     parseProcessParams(trackProc, funParams);
%                 end  
%                 if ismember(5,validProc)
%                     parseProcessParams(groupProc, funParams);
%                 end  
%             end
%             
%             % Set the tracking process index of the track grouping
%             if all(ismember([4 5], validProc))
%                 trackProcIndex = trackProc.getIndex();
%                 funParams.TrackProcessIndex = trackProcIndex;
%                 parseProcessParams(groupProc, funParams);
%             end
        end
        
    end
    
    methods (Static)
        
        function m = getDependencyMatrix(i,j)
               % 1 2 3 4 5 6 7
            m = [0 0 0 0 0 0 0 ; %1 - DarkCurrentCorrection 
                 2 0 0 0 0 0 0 ; %2 - Shade correction process
                 0 2 0 0 0 0 0 ; %3 - Thresholding        
                 0 0 2 0 0 0 0 ; %4 - background masks
                 0 2 0 2 0 0 0 ; %5 - background subtraction
                 0 0 0 0 0 0 0 ; %6 - adhesion seg
                 0 0 0 0 2 0 0 ];%7 - seg sampling

            if nargin<2, j=1:size(m,2); end
            if nargin<1, i=1:size(m,1); end
            m=m(i,j);
        end
        
        function name = getName()
            name='Focal Adhesion Segmentation';
        end
        function varargout = GUI(varargin)
            % Start the package GUI
            varargout{1} = focalAdhesionSegmentationPackageGUI(varargin{:});
        end
        function procConstr = getDefaultProcessConstructors(index)
            procContrs = {
                @DarkCurrentCorrectionProcess,...
                @ShadeCorrectionProcess,...
                @ThresholdProcess,...
                @BackgroundMasksProcess,...                
                @BackgroundSubtractionProcess,...
                @FocalAdhesionSegmentationProcess,...
                @SegmentationSamplingProcess};
            
            if nargin==0, index=1:numel(procContrs); end
            procConstr=procContrs(index);
        end
        function classes = getProcessClassNames(index)
            classes = {
                'DarkCurrentCorrectionProcess',...
                'ShadeCorrectionProcess',...
                'ThresholdProcess',...                
                'BackgroundMasksProcess',...
                'BackgroundSubtractionProcess',...
                'FocalAdhesionSegmentationProcess',...
                'SegmentationSamplingProcess'};
            
            if nargin==0 || isempty(index), index=1:numel(classes); end
            classes=classes(index);
        end
        
    end
end
