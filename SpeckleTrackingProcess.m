classdef SpeckleTrackingProcess < DataProcessingProcess
    % Concrete class for a speckle tracking process
    %
    % Sebastien Besson, May 2011
    
    methods
        function obj = SpeckleTrackingProcess(owner,outputDir, funParams)
            
            if nargin == 0
                super_args = {};
            else
                super_args{1} = owner;
                super_args{2} = SpeckleTrackingProcess.getName;
                super_args{3} = @trackMovieSpeckles;
                if nargin < 3 || isempty(funParams)
                    
                    %----Defaults----%
                    funParams.ChannelIndex = 1 : numel(owner.channels_);
                    funParams.OutputDirectory = [outputDir  filesep 'speckleTracks'];
                    funParams.enhanced = 0;                  
                    funParams.threshold = 3;
                    funParams.corrLength = 33;
                    funParams.interpolationMethod = 1;
                end
                super_args{4} = funParams;
            end
            
            obj = obj@DataProcessingProcess(super_args{:});
        end
        function sanityCheck(obj)
            
        end
        
         function varargout = loadChannelOutput(obj,iChan,varargin)
             
             % Input check
             outputList={'MPM','M','gapList','flow'};
             ip =inputParser;
             ip.addRequired('obj',@(x) isa(x,'SpeckleTrackingProcess'));
             ip.addRequired('iChan',@(x) ismember(x,1:numel(obj.owner_.channels_)));
             ip.addOptional('iFrame',1:obj.owner_.nFrames_,...
                @(x) ismember(x,1:obj.owner_.nFrames_));
             ip.addParamValue('output','MPM',@(x) all(ismember(x,outputList)));
             ip.parse(obj,iChan,varargin{:})
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
                             varargout{i}=varargout{i}(:,:,iFrame);
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
               @(x) ismember(x,1:numel(obj.owner_.channels_)));
           ip.parse(obj,varargin{:});
           iChan=ip.Results.iChan;
           
           %Makes sure there's at least one .mat file in the speified
           %directory
           status = arrayfun(@(x) exist(obj.outFilePaths_{x},'file'),iChan);           
         end
        
         function output = getDrawableOutput(obj)
            colors = hsv(numel(obj.owner_.channels_));
            output(1).name='Frame to frame displacement';
            output(1).var='M';
            output(1).formatData=@(x) [x(all(x(:,[1 3])~=0,2),[2 1])...
                x(all(x(:,[1 3])~=0,2),[4 3])-x(all(x(:,[1 3])~=0,2),[2 1])];
            output(1).type='overlay';
            output(1).defaultDisplayMethod=@(x)...
                 VectorFieldDisplay('Color',colors(x,:));
            output(2).name='Interpolated flow';
            output(2).var='flow';
            output(2).formatData=@formatFlow;
            output(2).type='overlay';
            output(2).defaultDisplayMethod=@(x)...
                VectorFieldDisplay('Color',colors(x,:));
            output(3).name='Tracks';
            output(3).var='MPM';
            output(3).formatData=@formatTracks;
            output(3).type='overlay';
            output(3).defaultDisplayMethod=@(x)...
                TracksDisplay('Color',colors(x,:),'showLabel',false);
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
    end
end


function flow = formatFlow(initFlow)
if isempty(initFlow), 
    flow=zeros(1,4);
else
    flow=[initFlow(:,[2 1]) initFlow(:,[4 3])-initFlow(:,[2 1])];
end
end


function fdata=formatTracks(data)
trackStart=arrayfun(@(x) find(data(x,:)==0,1,'last')+1,1:size(data,1),'Uniformout',false);
trackStart(cellfun(@isempty,trackStart))={1};
trackStart = [trackStart{:}];
fdata.x=arrayfun(@(x,y)data(x,y+1:2:end),1:size(data,1),trackStart,'UniformOutput',false);
fdata.y=arrayfun(@(x,y)data(x,y:2:end),1:size(data,1),trackStart,'UniformOutput',false);
end