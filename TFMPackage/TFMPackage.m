classdef TFMPackage < Package
    % Main class of the TFM package
    
    % Sebastien Besson, Aug 2011
    
    methods
        function obj = TFMPackage(owner,outputDir)
            % Constructor of class TFMPackage
            
            if nargin == 0
                super_args = {};
            else
                % Owner: MovieData object
                super_args{1} = owner;
                super_args{2} = TFMPackage.getDependencyMatrix;               
                super_args{3} = [outputDir  filesep 'TFMPackage'];
            end
            
            % Call the superclass constructor
            obj = obj@Package(super_args{:});
        end
        
        
        function processExceptions = sanityCheck(obj,varargin) 
            
            % Check that the first channel has a psf function
            psfSigmaCheck = isempty(obj.owner_.channels_(1).psfSigma_);
            if psfSigmaCheck
                error(['Missing standard deviation of the theoretical point-spread function! '...
                    'Please fill the numerical aperture, pixel size and'...
                    ' emission wavelengths for the beads channel!']);            
            end
            
            % Check that the movie has a frame rate
            if isempty(obj.owner_.timeInterval_)
                error('Missing frame rate! Please fill the time interval!');            
            end
            processExceptions = sanityCheck@Package(obj,varargin{:});
        end

    end
    
    methods (Static)
        
        function m = getDependencyMatrix()
            % Get dependency matrix
            %    1 2 3 4
            m = [0 0 0 0;
                2 0 0 0;
                0 1 0 0
                0 1 2 0;];
        end
        
        function name = getName()
            name = 'Traction Force Microscopy';
        end
        
        function id = getOptionalProcessId()
            M=TFMPackage.getDependencyMatrix;
            % Get the optional process id
            id=find(sum(M==2,1));
        end
        
        function varargout = GUI(varargin)
            % Start the package GUI
            varargout{1} = tfmPackageGUI(varargin{:});
        end
        
        function procConstr = getDefaultProcessConstructors(index)
            TFMProcConstr = {
                @StageDriftCorrectionProcess,...
                @DisplacementFieldCalculationProcess,...
                @DisplacementFieldCorrectionProcess,...
                @ForceFieldCalculationProcess};
            
            if nargin==0, index=1:numel(TFMProcConstr); end
            procConstr=TFMProcConstr(index);
        end
        function classes = getProcessClassNames(index)
            TFMClasses = {
                'StageDriftCorrectionProcess',...
                'DisplacementFieldCalculationProcess',...
                'DisplacementFieldCorrectionProcess',...
                'ForceFieldCalculationProcess'};
            if nargin==0, index=1:numel(TFMClasses); end
            classes=TFMClasses(index);
        end
    end
   
end