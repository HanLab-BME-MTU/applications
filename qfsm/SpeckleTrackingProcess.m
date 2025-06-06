classdef SpeckleTrackingProcess < DataProcessingProcess
    % Concrete class for a speckle tracking process
    %
    % Sebastien Besson, May 2011
    
    methods
        function obj = SpeckleTrackingProcess(owner,varargin)
            
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
                super_args{2} = SpeckleTrackingProcess.getName;
                super_args{3} = @trackMovieSpeckles;
                if isempty(funParams)
                    funParams=SpeckleTrackingProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            obj = obj@DataProcessingProcess(super_args{:});
        end
        
         function varargout = loadChannelOutput(obj,iChan,varargin)
             
             % Input check
             outputList={'MPM','M','gapList','flow'};
             ip =inputParser;
             ip.addRequired('iChan',@(x) isscalar(x) && obj.checkChanNum(x));
             ip.addOptional('iFrame',1:obj.owner_.nFrames_,@(x) all(obj.checkFrameNum(x)));
             ip.addParamValue('output',outputList,@(x) all(ismember(x,outputList)));
             ip.parse(iChan,varargin{:})
             iFrame=ip.Results.iFrame;
             
             % Data loading
             output = ip.Results.output;
             if ischar(output), output = {output}; end
             s = load(obj.outFilePaths_{1,iChan},output{:});
             
             for i=1:numel(output), varargout{i}=s.(output{i}); end
             
             % If single frame is selected
             if numel(iFrame)==1
                 for i=1:numel(output)
                     switch output{i}
                         case 'M'
                             if size(varargout{i}, 3) < iFrame,
                                 varargout{i} = zeros(1,4);
                             else
                                 varargout{i}=varargout{i}(:,:,iFrame);
                             end
                         case 'MPM'
                             index = all(varargout{i}(:,2*iFrame-1:2*iFrame),2)~=0;
                             varargout{i}=varargout{i}(index,1:2*iFrame);    
                             
                         case {'flow','gapList'}
                             varargout{i}=varargout{i}{iFrame};
                     end
                 end
             end
         end
                     
         function status = checkChannelOutput(obj,varargin)

           %Checks if the selected channels have valid output files
           ip =inputParser;
           ip.addRequired('obj',@(x) isa(x,'SpeckleTrackingProcess'));
           ip.addOptional('iChan',1:numel(obj.owner_.channels_),...
               @(x) all(obj.checkChanNum(x)));
           ip.parse(obj,varargin{:});
           iChan=ip.Results.iChan;
           
           %Makes sure there's at least one .mat file in the speified
           %directory
           status = arrayfun(@(x) exist(obj.outFilePaths_{x},'file'),iChan);           
         end
        
         function output = getDrawableOutput(obj)
            colors = hsv(numel(obj.owner_.channels_));
            output(1).name='Tracks';
            output(1).var='MPM';
            output(1).formatData=@formatTracks;
            output(1).type='overlay';
            output(1).defaultDisplayMethod=@(x)...
                TracksDisplay('Color',colors(x,:),'showLabel',false);
            output(2).name='Frame to frame displacement';
            output(2).var='M';
            output(2).formatData=@(x) [x(all(x(:,[1 3])~=0,2),[2 1])...
                x(all(x(:,[1 3])~=0,2),[4 3])-x(all(x(:,[1 3])~=0,2),[2 1])];
            output(2).type='overlay';
            output(2).defaultDisplayMethod=@(x)...
                 VectorFieldDisplay('Color',colors(x,:));
            output(3).name='Tracking initializers';
            output(3).var='flow';
            output(3).formatData=@formatFlow;
            output(3).type='overlay';
            output(3).defaultDisplayMethod=@(x)...
                VectorFieldDisplay('Color',[0 1 0]);

        end

    end
    methods (Static)
        function name =getName()
            name = 'Speckle Tracking';
        end
        function h = GUI()
            h= @speckleTrackingProcessGUI;
        end
        function methods = getInterpolationMethods(varargin)
            intMethods(1).name = 'none';
            intMethods(1).description = ' Do not interpolate';
            intMethods(2).name = 'nearest-neighbor';
            intMethods(2).description = ' Use nearest-neighbor flow';
            intMethods(3).name = 'gaussian';
            intMethods(3).description = ' Use Gaussian weighted neighbors interpolation';
            
            ip=inputParser;
            ip.addOptional('index',1:length(intMethods),@isvector);
            ip.parse(varargin{:});
            methods=intMethods(ip.Results.index);
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
            funParams.OutputDirectory = [outputDir  filesep 'speckleTracks'];
            funParams.enhanced = 0;
            funParams.threshold = 3;
            funParams.corrLength = 33;
            funParams.interpolationMethod = 1;
        end
    end
end


function flow = formatFlow(initFlow)
if isempty(initFlow), 
    flow=zeros(1,4);
else
    flow=[initFlow(:,[2 1]) initFlow(:,[4 3])-initFlow(:,[2 1])];
end
end


function tracks=formatTracks(data)
nTracks = size(data,1);

% Find track starts
trackStart=arrayfun(@(x) find(data(x,:)==0,1,'last')+1,1:nTracks,'UniformOutput',false);
trackStart(cellfun(@isempty,trackStart))={1};
trackStart = [trackStart{:}];

% Create a tracks array of structures with 2 fields xCoord and yCoord
tracks(nTracks, 1) =struct('xCoord',[],'yCoord',[]);
for i=1:nTracks
    tracks(i).xCoord = data(i,trackStart(i)+1:2:end);
    tracks(i).yCoord = data(i,trackStart(i):2:end);
end
end