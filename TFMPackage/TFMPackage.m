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
                super_args{2} = TFMPackage.getName;
                % Dependency Matrix (same length as process class name
                % string)
                super_args{3} = TFMPackage.getDependencyMatrix;
                
                % Process CLASS NAME string (same length as dependency matrix)
                % Must be accurate process class name
                TFMProcConstr = {
                    @StageDriftCorrectionProcess,...
                    @DisplacementFieldCalculationProcess,...
                    @DisplacementFieldCorrectionProcess,...
                    @ForceFieldCalculationProcess};
                super_args{4} = cellfun(@func2str,TFMProcConstr,...
                    'UniformOutput',false);
                super_args{5} = [outputDir  filesep 'TFMPackage'];
            end
            
            % Call the superclass constructor
            obj = obj@Package(super_args{:},...
                'processClassHandles_',TFMProcConstr);
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
        
        function varargout = start(varargin)
            % Start the package GUI
            varargout{1} = tfmPackageGUI(varargin{:});
        end
    end

    
end

