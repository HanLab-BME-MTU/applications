classdef QFSMPackage < Package
    % The main class of the QFSM package
    %
    % Sebastien Besson, 5/2011
    
    methods
        function obj = QFSMPackage(owner,outputDir)
            % Constructor of class QFSMPackage
            
            if nargin == 0
                super_args = {};
            else
                % Owner: MovieData object
                super_args{1} = owner;
                super_args{2} = [outputDir  filesep 'QFSMPackage'];
            end
            
            % Call the superclass constructor
            obj = obj@Package(super_args{:});
        end

        function [status processExceptions] = sanityCheck(obj,varargin) 
            
            % Check that the channels have a value for the spsf sigma
            psfSigmaCheck =arrayfun(@(x)isempty(x.psfSigma_),obj.owner_.channels_);
            if any(psfSigmaCheck)
                error(['Missing standard deviation of the theoretical point-spread function! '...
                    'Please fill the numerical aperture, pixel size and'...
                    ' emission wavelengths!']);            
            end
            % Check that the time interval is correctly setup
            if isempty(obj.owner_.timeInterval_)
                error('Missing frame rate! Please fill the time interval!');            
            end
            [status processExceptions] = sanityCheck@Package(obj,varargin{:});
        end

    end
    
    methods (Static)
        
        function m = getDependencyMatrix(i,j)
            
            %    1 2 3 4 5 6 7 8 9
            m = [0 0 0 0 0 0 0 0;  %1 NoiseEstimationProcess
                0 0 0 0 0 0 0 0;   %2 ThresholdProcess
                0 1 0 0 0 0 0 0;   %3 MaskRefinementProcess
                2 0 1 0 0 0 0 0;   %4 SpeckleDetectionProcess
                0 0 1 1 0 0 0 0;   %5 FlowTrackingProcess
                0 0 0 1 2 0 0 0;   %6 SpeckleTrackingProcess
                2 0 1 1 0 1 0 0;   %7 KineticAnalysisProcess
                0 0 1 0 0 1 0 0;]; %8 FlowAnalysisProcess
            if nargin<2, j=1:size(m,2); end
            if nargin<1, i=1:size(m,1); end
            m=m(i,j);
        end
        
        function name = getName()
            name = 'Quantitative Fluorescent Speckle Microscopy';
        end
        
        function varargout = GUI(varargin)
            % Start the package GUI
            varargout{1} = qfsmPackageGUI(varargin{:});
        end
        
        function funParams = getDefaultThresholdParams(owner,outputDir)
            funParams = ThresholdProcess.getDefaultParams(owner,outputDir);
            funParams.GaussFilterSigma = 1;
        end
        
        function funParams = getDefaultMaskRefinementParams(owner,outputDir)
            funParams = MaskRefinementProcess.getDefaultParams(owner,outputDir);
            funParams.ClosureRadius = 5;
        end
        
        function procConstr = getDefaultProcessConstructors(index)
            QFSMProcConstr = {
                @NoiseEstimationProcess,...
                @(x,y)ThresholdProcess(x,y,QFSMPackage.getDefaultThresholdParams(x,y)),...
                @(x,y)MaskRefinementProcess(x,y,QFSMPackage.getDefaultMaskRefinementParams(x,y)),... 
                @SpeckleDetectionProcess,...
                @FlowTrackingProcess,...
                @SpeckleTrackingProcess,...
                @KineticAnalysisProcess,...
                @FlowAnalysisProcess};
              
            if nargin==0, index=1:numel(QFSMProcConstr); end
            procConstr=QFSMProcConstr(index);
        end
        function classes = getProcessClassNames(index)
            QFSMClasses = {
                'NoiseEstimationProcess',...
                'ThresholdProcess',...
                'MaskRefinementProcess',...
                'SpeckleDetectionProcess',...
                'FlowTrackingProcess',...
                'SpeckleTrackingProcess',...
                'KineticAnalysisProcess',...
                'FlowAnalysisProcess'};         
            if nargin==0, index=1:numel(QFSMClasses); end
            classes=QFSMClasses(index);
        end
    end

    
end

