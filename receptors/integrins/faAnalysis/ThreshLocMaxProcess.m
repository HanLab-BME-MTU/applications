classdef ThreshLocMaxProcess < DetectionProcess
    % A concrete class for detecting large and diffraction-limited objects using thresholding and local maxima detection
    % Tony Vega 12/2014
    % Orginal Khuloud Jaqaman
    
    methods (Access = public)
        function obj = ThreshLocMaxProcess(owner, varargin)
            % Input check
            ip = inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.addOptional('funParams',[],@isstruct);
            ip.parse(owner,varargin{:});
            outputDir = ip.Results.outputDir;
            funParams = ip.Results.funParams;
            
            
            % Constructor of the SubResolutionProcess
            
            super_args{1} = owner;
            super_args{2} = ThreshLocMaxProcess.getName;
            super_args{3} = @detectMovieThreshLocMax;
            if isempty(funParams)  % Default funParams
                funParams = ThreshLocMaxProcess.getDefaultParams(owner,outputDir);
            end
            super_args{4} = funParams;
            if(nargin > 4)
                super_args{5:nargin} = varargin{5:nargin};
            end
            obj = obj@DetectionProcess(super_args{:});            
        end
        function varargout = loadChannelOutput(obj,iChan,varargin)
            
            % Input check
            outputList = {'movieInfo','pixelIndx'};
            ip =inputParser;
            ip.StructExpand = true;
            ip.addRequired('iChan',@(x) isscalar(x) && obj.checkChanNum(x));
            ip.addOptional('iFrame',1:obj.owner_.nFrames_,@(x) all(obj.checkFrameNum(x)));
            ip.addParamValue('useCache',false,@islogical);
            ip.addParamValue('output',outputList,@(x) all(ismember(x,outputList)));
            ip.parse(iChan,varargin{:})
            iFrame = ip.Results.iFrame;
            output = ip.Results.output;
            if ischar(output),output={output}; end
            
            % Data loading
            % load outFilePaths_{1,iChan}
            %
            s = cached.load(obj.outFilePaths_{1,iChan}, '-useCache', ip.Results.useCache, output{:});
           
            if numel(ip.Results.iFrame)>1,
                varargout{1}=s.(output{1});
            else
                switch(output{1})
                    case 'movieInfo'
                        if numel(ip.Results.iFrame)>1,
                            varargout{1}=s.(output{1});
                        else
                            varargout{1}=s.(output{1})(iFrame);
                        end
                    case 'pixelIndx'
                        if(iFrame <= obj.funParams_.lastImageNum)
                            varargout{1}=s.(output{1})(:,:,iFrame);
                        else
                            varargout{1} = zeros(size(s.(output{1})(:,:,1)));
                        end
                end
            end
        end
         function output = getDrawableOutput(obj)
%             nOutput = 2;
% %             colors = parula(numel(obj.owner_.channels_)*nOutput);
            % Rename default detection output
            output(1) = getDrawableOutput@DetectionProcess(obj);
            output(1).name='Sub-resolution objects';
            
            output(2).name='Masks';
            output(2).var='pixelIndx';            
            output(2).type='overlay';
%             if isempty(obj.maxIndex_)            
            output(2).formatData=@MaskProcess.getMaskBoundaries;
            colors = hsv(numel(obj.owner_.channels_));
%             output(2).defaultDisplayMethod=@(x) LineDisplay('Color',colors(x,:));
            output(2).defaultDisplayMethod=@(x) LineDisplay('Color','m');
%             else
%                 cMap = randomColormap(obj.maxIndex_,42);%random colors since index is arbitrary, but constant seed so we are consistent across frames.
%                 output(2).formatData=@(x)(MaskProcess.getMaskOverlayImage(x,cMap));
%                 %If index values for diff channels are to be differentiated
%                 %they must be done at the level of the indexes.
%                 output(2).defaultDisplayMethod=@ImageOverlayDisplay;
%                 
%                 output(3).name='Object Number';
%                 output(3).var = 'number';
%                 output(3).type = 'overlay';
%                 output(3).formatData=@(x)(MaskProcess.getObjectNumberText(x,cMap));
%                 output(3).defaultDisplayMethod=@TextDisplay;
%                 
%                 
%             end
% %             output(1).name='Meshwork';
% %             output(1).var='S4';
% %             output(1).formatData=@(x) x.getEdgeXY;
% %             output(1).type='overlay';
% %             output(1).defaultDisplayMethod=@(x) LineDisplay('Marker','none',...
% %                 'Color',colors((x-1)*nOutput+1,:));
% %             
% %             output(2).name = 'Junctions';
% %             output(2).var='S4v';
% %             output(2).formatData=@(x) x.getVertexXY;
% %             output(2).type='overlay';
% %             output(2).defaultDisplayMethod=@(x) LineDisplay('Marker','x',...
% %                 'LineStyle','none','Color',colors((x-1)*nOutput+2,:));

        end 
%         function output = getDrawableOutput(obj)
%             % Rename default detection output
%             output = getDrawableOutput@DetectionProcess(obj);
%             output(1).name='Sub-resolution objects';
%         end
        
    end
    methods (Static)
        
        function name = getName()
            name = 'ClusterDetection';
        end

        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            % moviePara  
            funParams.ChannelIndex =1:numel(owner.channels_);
            funParams.OutputDirectory = [outputDir  filesep 'ClusterDetections'];
            funParams.firstImageNum = 20;
            funParams.lastImageNum = owner.nFrames_;
            
            % detectionParam
            funParams.detectionParam.alphaLocMax = 0.05;
            funParams.detectionParam.alphaSubRes = 0.05;
            funParams.detectionParam.psfSigma = 1;            
            funParams.detectionParam.thresholdMethod = 'rosin';
            funParams.detectionParam.methodValue = [];
            funParams.detectionParam.filterNoise = 2;
            funParams.detectionParam.filterBackground = 10;
            funParams.detectionParam.minSize = 10;
            funParams.detectionParam.maxSize = 10000;
            
        end

    end
    
end