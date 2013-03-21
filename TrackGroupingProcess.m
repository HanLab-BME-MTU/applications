classdef TrackGroupingProcess < DataProcessingProcess
    % A concrete class associated to the grouping of sister tracks
    %
    % Sebastien Besson, March 2012

    methods (Access = public)
        function obj = TrackGroupingProcess(owner, varargin)            
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
                
                super_args{1} = owner;
                super_args{2} = TrackGroupingProcess.getName;
                super_args{3} = @groupMovieTracks;
                if isempty(funParams)  % Default funParams
                    funParams = TrackGroupingProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;  
            end
            
            obj = obj@DataProcessingProcess(super_args{:});
        end
        function varargout = loadChannelOutput(obj,iChan,varargin)
            
            % Input check
            outputList = {'segments','tracks'};
            ip =inputParser;
            ip.addRequired('iChan',@(x) isscalar(x) && obj.checkChanNum(x));
            ip.addOptional('iFrame',[],@(x) all(obj.checkFrameNum(x)));
            ip.addParamValue('output',outputList,@(x) all(ismember(x,outputList)));
            ip.parse(iChan,varargin{:})
            iFrame = ip.Results.iFrame;
            output = ip.Results.output;
            if ischar(output),output={output}; end
            
            % Data loading
            for j=1:numel(output)
                if ismember(output{j},'tracks')
                    s = load(obj.outFilePaths_{1,iChan});
                    tracksFinal = s.tracksFinal;
                    
                    if ~isempty(ip.Results.iFrame),
                        % Filter tracks existing in input frame
                        trackSEL=getTrackSEL(tracksFinal);
                        validTracks = (iFrame>=trackSEL(:,1) &iFrame<=trackSEL(:,2));
                        [tracksFinal(~validTracks).tracksCoordAmpCG]=deal([]);
                        
                        for i=find(validTracks)'
                            tracksFinal(i).tracksCoordAmpCG = tracksFinal(i).tracksCoordAmpCG(1:8*(iFrame-trackSEL(i,1)+1));
                            tracksFinal(i).label = s.trackLabels(i);
                        end
                    else
                        for i=1:numel(tracksFinal), 
                            tracksFinal(i).label = s.trackLabels(i);
                        end
                    end
                    varargout{j} = tracksFinal;
                elseif ismember(output{j},'segments')
                    s = load(obj.outFilePaths_{2,iChan},'segments');
                    if ~isempty(ip.Results.iFrame),
                        varargout{j} = s.segments{iFrame};
                    else
                        varargout{j} = s.segments;
                    end
                end
            end


        end       
        
         function output = getDrawableOutput(obj)
            colors = hsv(numel(obj.owner_.channels_));
            output(1).name='Classified tracks';
            output(1).var='tracks';
            output(1).formatData=@TrackingProcess.formatTracks;
            output(1).type='overlay';
            output(1).defaultDisplayMethod=@(x)TracksDisplay('Color',hsv(64));
            output(2).name='Classified segments';
            output(2).var='segments';
            output(2).formatData=@(x) reshape([x(:,1:4) NaN(size(x,1),2)]',2,3*size(x,1))';
            output(2).type='overlay';
            output(2).defaultDisplayMethod=@(x)LineDisplay('Color',colors(x,:));

        end
    end
    methods (Static)
        
        function name = getName()
            name = 'Track Grouping';
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex = 1 : numel(owner.channels_);
            funParams.OutputDirectory = [outputDir  filesep 'groupedTracks'];
            funParams.DetProcessIndex = [];
            funParams.TrackProcessIndex = [];
            funParams.MaskProcessIndex = [];

            funParams.minLifetime = 1;
            funParams.maxDistance = 2000;
            funParams.minOverlap = 1;
            funParams.bandWidth = 1000;
            funParams.minDistance = 350;
            funParams.alpha = .05;

        end
    end    
end