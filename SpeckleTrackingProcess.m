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

             % Data loading
             output = ip.Results.output;
             if ischar(output), output = {output}; end
             s = load(obj.outFilePaths_{1,iChan},output{:});
             
             for i=1:numel(output), varargout{i}=s.(output{i}); end
             
             % If single frame is selected
             if numel(ip.Results.iFrame)==1
                 for i=1:numel(output)
                     switch output{i}
                         case 'M'
                             varargout{i}=varargout{i}(:,:,ip.Results.iFrame);
                         case {'flow','gapList'}
                             varargout{i}=varargout{i}{ip.Results.iFrame};
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
            output(1).name='Tracks';
            output(1).var='M';
            output(1).formatData=@formatTracks;
            output(1).type='overlay';
            output(1).defaultDisplayMethod=@(x)...
                VectorFieldDisplay('Color',colors(x,:));
            output(2).name='Interpolated flow';
            output(2).var='flow';
            output(2).formatData=@formatFlow;
            output(2).type='overlay';
            output(2).defaultDisplayMethod=@(x)...
                VectorFieldDisplay('Color',colors(x,:));
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


function y=formatTracks(x)
ind = (x(:,1)~=0 & x(:,3)~=0);
y=[x(ind,[2 1]) x(ind,[4 3])-x(ind,[2 1])];
end