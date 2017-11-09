classdef PartitionAnalysisProcess < DataProcessingProcess
    % MovieData process that stores information about the paritioning fraction of tracks
    %
    % funParams : function parameters for maskDetectedStrucutre similar
    %                   to other DataProcessingProcesses
    %   .channel_tracks     : channel index of the MD that contains the track
    %                       information
    %   .channel_struct     : channel index of the MD that contains the mask
    %                       information
    %   .outputDirectory    : directory where partitioning information will
    %                       be stored
    %   .scrambleTracks     : randomize track location to create control
    %                       data set <Unused>
    %   .nControl           : number of times the randomized control
    %                         (scrambling track location) is repeated.
    
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
            %funParams.scrambleTracks = false;
            funParams.nControl = 5;
        end
        function name = getName()
            name = 'PartitionAnalysisProcess';
        end
    end
end
