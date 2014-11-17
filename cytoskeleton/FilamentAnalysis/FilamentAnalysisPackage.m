classdef FilamentAnalysisPackage < Package
    % A concrete package for segmentation of filament and the orientation
    % Liya Ding 06. 2012.
    % Liya is still working on it. It is not ready yet.
    
    methods (Access = public)
        function obj = FilamentAnalysisPackage(owner,varargin)
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
                super_args{2} = [outputDir  filesep 'FilamentAnalysisPackage'];
            end
            % Call the superclass constructor
            obj = obj@Package(super_args{:});
        end
        
        function [status processExceptions] = sanityCheck(obj,varargin)
            [status processExceptions] = sanityCheck@Package(obj,varargin{:});
            
        end
    end
    
    methods (Static)
        
        function name = getName()
            name = 'FilamentAnalysis';
        end
        
        function m = getDependencyMatrix(i,j)
            m = [0 0 0 0 0 ;
                 1 0 0 0 0 ;
                 0 0 0 0 0 ;
                 0 0 0 0 0 ;
                 0 0 0 0 0 ];

        if nargin<2, j=1:size(m,2); end
            if nargin<1, i=1:size(m,1); end
            m=m(i,j);
        end
        
        function varargout = GUI(varargin)
            % Start the package GUI
            varargout{1} = FilamentAnalysisPackageGUI(varargin{:});
        end
        
        function procConstr = getDefaultProcessConstructors(index)
            FilemantProcConstr = {
                @ThresholdProcess, ...                
                @MaskRefinementProcess, ...
                @ImageFlattenProcess, ... 
                @SteerableFilteringProcess, ...
                @FilamentSegmentationProcess ...                
                };
            if nargin==0, index=1:numel(FilemantProcConstr); end
            procConstr=FilemantProcConstr(index);
        end
        
        function classes = getProcessClassNames(index)
            FilamentClasses = {'SegmentationProcess', ...
                'MaskRefinementProcess', ...
                'ImageFlattenProcess', ...
                'SteerableFilteringProcess', ...
                'FilamentSegmentationProcess' ...
                };
            
            if nargin==0, index=1:numel(FilamentClasses); end
            classes=FilamentClasses(index);
        end
        
    end
    
end