classdef TrackPartitionProcess < DataProcessingProcess & NonSingularProcess
    %TRACKPARTITIONPROCESS Process for trackPartition
    properties
        trackChannel_ % Channel in owner MD with tracks to partition
        maskChannel_ % Channel of maskMD with particles to use for mask
        % MD with particles used to mask image
        maskMovieDataPath_
        maskMovieDataFileName_
        mask_ 
    end
    
    methods (Access = public)
        % Constructor
        function obj = TrackPartitionProcess(owner,varargin)
            nChannels = numel(owner.channels_);
            channels = 1:nChannels;
                       
            ip = inputParser;
            % Movie whose tracks to partition
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            % Movie containing structure/particle for masking image
            ip.addRequired('maskMD',@(x) isa(x,'MovieData'));
            ip.addRequired('trackChannel',@(x) isa(x,'numeric'));
            ip.addRequired('maskChannel',@(x) isa(x,'numeric'));
            ip.addOptional('funParams',[],@isstruct);
            ip.parse(owner,varargin{:});
            funParams = ip.Results.funParams;
            
            obj = obj@DataProcessingProcess(owner,'TrackPartitionProcess',@trackPartitionWrapper);
            obj.trackChannel_ = ip.Results.trackChannel;
            obj.maskChannel_ = ip.Results.maskChannel;
            
            maskMD = ip.Results.maskMD;
            obj.maskMovieDataPath_ = maskMD.movieDataPath_;
            obj.maskMovieDataFileName_ = maskMD.movieDataFileName_;
            
            outFileDir = [obj.owner_.outputDirectory_,...
                filesep,'TrackPartition'];
            obj.outFilePaths_{obj.trackChannel_} = [outFileDir,filesep,...
                'Channel-',num2str(obj.trackChannel_),...
                '-partition-result.mat'];
            if ~exist(outFileDir,'dir')
                mkdir(outFileDir)
            end
            
            % Use default parameters if none are specified
            if isempty(funParams)
                obj.funParams_ = obj.getDefaultParams(owner);
            else
                obj.funParams_ = funParams;
            end
        end       
    end
    
    methods (Static)
        function params = getDefaultParams(owner)
            params.minTrackLength = 3; % 3 frames
            params.minMaskDiam = 100; % 100 nm
            params.gaussianThresh = 0; % default to constant size masking
            params.upscale = 1;
            params.analysisStart = 1; % start at first frame
            params.analysisEnd = owner.nFrames_; % end at last frame
            params.runDiffAnalysis = true; % run diffusion analysis afterwards
        end
        
        function name = getName()
            name = 'Partitioning';
        end
    end
    
end

