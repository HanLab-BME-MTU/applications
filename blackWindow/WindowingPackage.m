classdef WindowingPackage < Package
    % The main class of the Windowing package
    %
    % Sebastien Besson, July 2011
    
    methods
        function obj = WindowingPackage(owner,outputDir)
            % Constructor of class QFSMPackage
            
            if nargin == 0
                super_args = {};
            else
                % Owner: MovieData object
                super_args{1} = owner;
                super_args{2} = 'Windowing';
                % Dependency Matrix (same length as process class name
                % string)
                super_args{3} = WindowingPackage.getDependencyMatrix;
                
                % Process CLASS NAME string (same length as dependency matrix)
                % Must be accurate process class name
                WindowingProcConstr = {
                    @ProtrusionProcess,...
                    @WindowingProcess,...
                    @ProtrusionSamplingProcess,...
                    @WindowSamplingProcess};
                super_args{4} = cellfun(@func2str,WindowingProcConstr,...
                    'UniformOutput',false);
                super_args{5} = [outputDir  filesep 'WindowingPackage'];
            end
            
            % Call the superclass constructor
            obj = obj@Package(super_args{:},...
                'processClassHandles_',WindowingProcConstr);
        end
    end
    
    methods (Static)
        
        function m = getDependencyMatrix()
            % Get dependency matrix
            
            %    1 2 3 4
            m = [0 0 0 0; %1 ProtrusionProcess
                 2 0 0 0; %2 WindowingProcess
                 1 1 0 0;  %3 ProtrusionSamplingProcess
                 0 1 0 0;]; %4 WindowSamplingProcess
        end
        
        function id = getOptionalProcessId()
            M=WindowingPackage.getDependencyMatrix;
            % Get the optional process id
            id=find(sum(M==2,1));
        end
        
        function varargout = start(varargin)
            % Start the package GUI
            varargout{1} = windowingPackageGUI(varargin{:});
        end
        
    end

    
end

