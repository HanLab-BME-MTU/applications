classdef WindowingPackage < Package
    % The main class of the Windowing package
    
    % Sebastien Besson, July 2011
    
    methods
        function obj = WindowingPackage(owner,varargin)
            % Constructor of class QFSMPackage
            
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
                super_args{2} = [outputDir  filesep 'WindowingPackage'];
            end
            
            % Call the superclass constructor
            obj = obj@Package(super_args{:});
        end
        
        function parentID = getParent(obj,procID)
            % Use default getParent method
            parentID=getParent@Package(obj,procID);
            
            % Refine dependency relationship between protrusion and
            % windowing processes
            if procID==2
                if isempty(obj.processes_{procID})
                    parentID(parentID==1)=[];
                elseif ~strcmp(obj.processes_{procID}.funParams_.MethodName,...
                        'ProtrusionBased')
                    parentID(parentID==1)=[];
                end
            end
        end
    end
    
    methods (Static)
        
        function m = getDependencyMatrix(i,j)
            
            m = [0 0 0 0; %1 ProtrusionProcess
                2 0 0 0; %2 WindowingProcess
                1 1 0 0;  %3 ProtrusionSamplingProcess
                0 1 0 0;]; %4 WindowSamplingProcess
            if nargin<2, j=1:size(m,2); end
            if nargin<1, i=1:size(m,1); end
            m=m(i,j);
        end
        
        function name = getName()
            name = 'Windowing';
        end
        
        function varargout = GUI(varargin)
            % Start the package GUI
            varargout{1} = windowingPackageGUI(varargin{:});
        end
        
        function procConstr = getDefaultProcessConstructors(index)
            windowingProcConstr = {
                @ProtrusionProcess,...
                @WindowingProcess,...
                @ProtrusionSamplingProcess,...
                @WindowSamplingProcess};
            
            if nargin==0, index=1:numel(windowingProcConstr); end
            procConstr=windowingProcConstr(index);
        end
        function classes = getProcessClassNames(index)
            windowingClasses = {
                'ProtrusionProcess',...
                'WindowingProcess',...
                'ProtrusionSamplingProcess',...
                'WindowSamplingProcess'};
            if nargin==0, index=1:numel(windowingClasses); end
            classes=windowingClasses(index);
        end
        
    end
    
    
end

