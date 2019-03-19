classdef OutputTFMProcess < DoubleProcessingProcess
    
    %A class for creating traction map in tiff format
    %
    %Sangyoon Han,
    %3/2019
    
    methods (Access = public)
        
        function obj = OutputTFMProcess(owner,varargin)
            % Input check
            ip = inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.addOptional('funParams',[],@isstruct);
            ip.parse(owner,varargin{:});
            outputDir = ip.Results.outputDir;
            funParams = ip.Results.funParams;
            
            super_args{1} = owner;
            super_args{2} = OutputTFMProcess.getName;
            super_args{3} = @createTractionMapTiff;
            
            if isempty(funParams)
                funParams=OutputTFMProcess.getDefaultParams(owner,outputDir);
            end
            
            super_args{4} = funParams;

            obj = obj@DoubleProcessingProcess(super_args{:});
        end
        
        function status = checkChannelOutput(obj,iChan)
            
            %Checks if the selected channels have valid output images
            nChanTot = numel(obj.owner_.channels_);
            if nargin < 2 || isempty(iChan)
                iChan = 1:nChanTot;
            end
            
            status =  arrayfun(@(x)(ismember(x,1:nChanTot) && ...
                (length(dir([obj.outFilePaths_{1,x} filesep '*.tif']))...
                == obj.owner_.nFrames_)),iChan);
        end
        
        function outIm = loadChannelOutput(obj,iChan,iFrame,varargin)
            
            ip =inputParser;
            ip.addRequired('obj');
            ip.addRequired('iChan',@(x) ismember(x,1:numel(obj.owner_.channels_)));
            ip.addRequired('iFrame',@(x) ismember(x,1:obj.owner_.nFrames_));
            ip.addParamValue('output',[],@ischar);
            ip.parse(obj,iChan,iFrame,varargin{:})
            
            %get the image names
            imNames = getOutImageFileNames(obj,iChan);
            outIm = double(imread([obj.outFilePaths_{1,iChan} filesep imNames{1}{iFrame}]));
            
        end
        
        
        function fileNames = getOutImageFileNames(obj,iChan)
            if obj.checkChannelOutput(iChan)
                fileNames = cellfun(@(x)(imDir(x)),obj.outFilePaths_(1,iChan),'UniformOutput',false);
                fileNames = cellfun(@(x)(arrayfun(@(x)(x.name),x,'UniformOutput',false)),fileNames,'UniformOutput',false);
                nIm = cellfun(@(x)(length(x)),fileNames);
                if ~all(nIm == obj.owner_.nFrames_)
                    error('Incorrect number of images found in one or more channels!')
                end
            else
                error('Invalid channel numbers! Must be positive integers less than the number of image channels!')
            end
            
            
        end
        function figHan = resultDisplay(obj)
            
            figHan = msgbox(['The TFM images have been saved as .tif images to the folder "' ...
                obj.funParams_.OutputDirectory '". They can be viewed with ImageJ or a comparable image viewing program.']);
            
        end
        
    end
    methods(Static)
        
        function name = getName()
            name = 'Traction Magnitude Output';
        end
        function h = GUI()
            h= @outputTFMProcessGUI;
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.OutputDirectory =  [outputDir  filesep 'TractionMapTiff'];
            funParams.useCellConfig = true;
            funParams.useRefConfig = false;
            
        end
    end
end