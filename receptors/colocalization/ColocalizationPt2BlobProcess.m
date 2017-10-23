classdef ColocalizationPt2BlobProcess < ImageAnalysisProcess
   % A concreate class for measuring colocalization between two images
   % Anthony Vega 09/2014
        methods (Access = public)
        function obj = ColocalizationPt2BlobProcess(owner,varargin)
            
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
                super_args{2} = ColocalizationPt2BlobProcess.getName;
                super_args{3} = @colocalizationPt2BlobWrapper;
                if isempty(funParams)
                    funParams = ColocalizationPt2BlobProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
                if(nargin > 4)
                    super_args{5:nargin} = varargin{5:nargin};
                end
            end
            
            obj = obj@ImageAnalysisProcess(super_args{:});
        end
        
    end
    methods (Static)
        function name = getName()
            name = 'ColocalizationPt2Blob';
        end

        function methods = getMethods(varargin)
            colocalizationMethods(1).name = 'Point2Blob';
            colocalizationMethods(1).func = @colocalMeasurePt2Blob;            
            
            ip=inputParser;
            ip.addOptional('index',1:length(colocalizationMethods),@isvector);
            ip.parse(varargin{:});
            index = ip.Results.index;
            methods=colocalizationMethods(index);
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelPt = 1;
            funParams.ChannelBlob = 2;
            funParams.ChannelMask = 2;
            funParams.SearchRadius = 3;
            funParams.OutputDirectory = [outputDir  filesep];
            funParams.ProcessIndex = [];
        end
    end
end