classdef PartitionAnalysisProcess < DataProcessingProcess
    % untitled2 Add summary here
    %
    % This template includes the minimum set of functions required
    % to define a System object with discrete state.
    
    properties
        % Public, tunable properties.
    end
   
    methods (Access = public)
        function obj = PartitionAnalysisProcess(owner, varargin)
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
                super_args{2} = PartitionAnalysisProcess.getName;
                super_args{3} = @trackPartitioning;
                if isempty(funParams)
                    funParams = PartitionAnalysisProcess.getDefaultParams(outputDir);
                end
                super_args{4} = funParams;
            end
            obj = obj@DataProcessingProcess(super_args{:});
        end
    end
    methods
      %for set.  
    end
    methods (Static)
        function funParams = getDefaultParams(outputDir)
            %funParams.
            funParams.channel_tracks = 1;
            funParams.channel_struct = 1;
            funParams.outputDirectory = [outputDir filesep 'TrackingPackage' filesep 'PartitionAnalysis'];
        end
        function name = getName()
            name = 'PartitionAnalysisProcess';
        end
    end
end
