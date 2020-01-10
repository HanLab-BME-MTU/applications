classdef TFMPackage < Package
    % Main class of the TFM package
    
    % Sebastien Besson, Aug 2011
    
    methods
        function obj = TFMPackage(owner,varargin)
            % Constructor of class TFMPackage
            
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
                super_args{2} = [outputDir  filesep 'TFMPackage'];
            end
            
            % Call the superclass constructor
            obj = obj@Package(super_args{:});
        end
        
        
        function [status processExceptions] = sanityCheck(obj,varargin) 
            
            % Check that the at least one channel has a valid psf sigma
%             hasPsfSigma = any(~cellfun(@isempty, {obj.owner_.channels_.psfSigma_}));
%             if ~hasPsfSigma
%                 error(['Missing standard deviation of the theoretical point-spread function! '...
%                     'Please fill the numerical aperture, pixel size and'...
%                     ' emission wavelengths for the beads channel!']);            
%             end
            
            % Check that the movie has a frame rate
            if isempty(obj.owner_.timeInterval_) && obj.owner_.nFrames_>1
                error('Missing frame rate! Please fill the time interval!');            
            end
            [status, processExceptions] = sanityCheck@Package(obj,varargin{:});
        end

    end
    
    methods (Static)
        function name = getName()
            name = 'Traction Force Microscopy';
        end
        
        function m = getDependencyMatrix(i,j)
  
            m = [0 0 0 0 0 0;   %1 StageDriftCorrectionProcess [optional]
                 2 0 0 0 0 0;   %2 DisplacementFieldCalculationProcess
                 0 1 0 0 0 0;   %3 DisplacementFieldCorrectionProcess [optional]
                 0 1 2 0 0 0;   %4 ForceFieldCalculationProcess
                 0 1 2 1 0 0;   %5 StrainEnergyCalculationProcess [optional]
                 0 1 2 1 0 0];  %6 OutputTFMProcess [optional]
            if nargin<2, j=1:size(m,2); end
            if nargin<1, i=1:size(m,1); end
            m=m(i,j);
        end
        
        function varargout = GUI(varargin)
            % Start the package GUI
            varargout{1} = tfmPackageGUI(varargin{:});
        end
        
        function procConstr = getDefaultProcessConstructors(index)
            TFMProcConstr = {
                @EfficientSubpixelRegistrationProcess,...
                @DisplacementFieldCalculationProcess,...
                @DisplacementFieldCorrectionProcess,...
                @ForceFieldCalculationProcess,...
                @StrainEnergyCalculationProcess,...
                @OutputTFMProcess};
            
            if nargin==0, index=1:numel(TFMProcConstr); end
            procConstr=TFMProcConstr(index);
        end
        function classes = getProcessClassNames(index)
            TFMClasses = {
                'StageDriftCorrectionProcess',...
                'DisplacementFieldCalculationProcess',...
                'DisplacementFieldCorrectionProcess',...
                'ForceFieldCalculationProcess',...
                'StrainEnergyCalculationProcess'...
                'OutputTFMProcess'};
            if nargin==0, index=1:numel(TFMClasses); end
            classes=TFMClasses(index);
        end
    end
   
end