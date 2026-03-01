classdef TFMPackage < Package
    % Main class of the TFM package
    %
    % Patched to:
    %   1) Fix a long-standing typo in getProcessClassNames (missing comma)
    %   2) Add an optional process: OtherChannelSamplingProcess (new step 7)
    %   3) Provide backward compatibility with legacy MovieData that stored
    %      TFMPackage with only 6 process slots.
    %
    % Sebastien Besson, Aug 2011
    % Patch: Sangyoon Han + ChatGPT (2026)

    methods
        function obj = TFMPackage(owner,varargin)
            % Constructor of class TFMPackage

            if nargin == 0
                super_args = {};
            else
                % Check input
                ip = inputParser;
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

        function [status, processExceptions] = sanityCheck(obj,varargin)
            % Backward compatibility: pad legacy process cell array
            % (Older saved MovieData may have TFMPackage with 6 slots)
            nExpected = numel(TFMPackage.getProcessClassNames());
            if numel(obj.processes_) < nExpected
                obj.processes_{nExpected} = [];
            end

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
            % 1 StageDriftCorrectionProcess [optional]
            % 2 DisplacementFieldCalculationProcess
            % 3 DisplacementFieldCorrectionProcess [optional]
            % 4 ForceFieldCalculationProcess
            % 5 StrainEnergyCalculationProcess [optional]
            % 6 OutputTFMProcess [optional]
            % 7 OtherChannelSamplingProcess [optional]
            m = [0 0 0 0 0 0 0;   %1
                 2 0 0 0 0 0 0;   %2 depends on 1
                 0 1 0 0 0 0 0;   %3 depends on 2
                 0 1 2 0 0 0 0;   %4 depends on 2 and 3
                 0 1 2 1 0 0 0;   %5 depends on 2,3,4
                 0 1 2 1 0 0 0;   %6 depends on 2,3,4
                 0 0 0 0 0 0 0];  %7 independent (sampling only)
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
                @OutputTFMProcess,...
                @OtherChannelSamplingProcess};

            if nargin==0, index=1:numel(TFMProcConstr); end
            procConstr=TFMProcConstr(index);
        end

        function classes = getProcessClassNames(index)
            % IMPORTANT: commas between entries (no accidental string concat)
            TFMClasses = {
                'StageDriftCorrectionProcess',...
                'DisplacementFieldCalculationProcess',...
                'DisplacementFieldCorrectionProcess',...
                'ForceFieldCalculationProcess',...
                'StrainEnergyCalculationProcess',...
                'OutputTFMProcess',...
                'OtherChannelSamplingProcess'};
            if nargin==0, index=1:numel(TFMClasses); end
            classes=TFMClasses(index);
        end
    end
end
