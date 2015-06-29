classdef SubcellMaskProcess < DataProcessingProcess
    % untitled3 Add summary here
    %
    % This template includes the minimum set of functions required
    % to define a System object with discrete state.
    
    properties
        imSize_;% Image size 1x2 array[height width]
        nFrames_;% Number of frames
    end
    
    methods (Access = public)
        %The constructor
        function obj = SubcellMaskProcess(owner, varargin)
            if nargin == 0
                super_args = {};
            else
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieData'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                % Define arguments for superclass constructor
                super_args{1} = owner;
                super_args{2} = SubcellMaskProcess.getName;
                super_args{3} = @maskDetectedStructure;
                if isempty(funParams)
                    funParams = SubcellMaskProcess.getDefaultParams(outputDir);
                end
                super_args{4} = funParams;
            end
            obj = obj@DataProcessingProcess(super_args{:});
        end
        %loads the output
        function mask = loadChannelOutput(obj,iChan,varargin)
            load(obj.outFilePaths_{1,iChan});
        end
    end
    methods
        function set.imSize_(obj, values)
            obj.imSize_ = values;
        end
        function set.nFrames_(obj, value)
            obj.nFrames_ = value;
        end
    end
    methods (Static)
        function funParams = getDefaultParams(outputDir)
            funParams.channel = 1;
            funParams.psfSigmaMult = 3;
            funParams.outputDirectory = [outputDir filesep 'TrackingPackage' filesep 'StructureMasking'];
            %funParams.psfSigma = 1.191;
        end
        function name = getName()
            name = 'SubcellMaskProcess';
        end
    end
end
