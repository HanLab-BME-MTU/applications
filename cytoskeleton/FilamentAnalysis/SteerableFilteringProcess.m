classdef SteerableFilteringProcess < ImageProcessingProcess
    % A concrete class for steerable filtering
    %
    % Liya Ding, 06. 2012
    
    methods (Access = public)
        
        function obj = SteerableFilteringProcess(owner,varargin)
            
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
                super_args{2} = SteerableFilteringProcess.getName;
                super_args{3} = @steerable_filter_forprocess;
                if isempty(funParams)
                    funParams = SteerableFilteringProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
                
                if nargin > 4
                    super_args{5} = inImagePaths;
                end
                if nargin > 5
                    super_args{6} = outImagePaths;
                end
                
            end
            
            obj = obj@ImageProcessingProcess(super_args{:});
            
        end
        
        function setInImagePath(obj,chanNum,imagePath)
            
            if ~obj.checkChanNum(chanNum)
                error('lccb:set:fatal','Invalid image channel number for image path!\n\n');
            end
            
            if ~iscell(imagePath)
                imagePath = {imagePath};
            end
            nChan = length(chanNum);
            if nChan ~= length(imagePath)
                error('lccb:set:fatal','You must specify a path for every channel!')
            end
            
            for j = 1:nChan
                if ~exist(imagePath{j},'dir')
                    error('lccb:set:fatal',...
                        ['The directory specified for channel ' ...
                        num2str(chanNum(j)) ' is invalid!'])
                else
                    if isempty(imDir(imagePath{j})) && ...
                            isempty(dir([imagePath{j} filesep '*.mat']))
                        error('lccb:set:fatal',...
                            ['The directory specified for channel ' ...
                            num2str(chanNum(j)) ' does not contain any images!!'])
                    else
                        obj.inFilePaths_{1,chanNum(j)} = imagePath{j};
                    end
                end
            end
        end
        
        function fileNames = getOutImageFileNames(obj,iChan)
            if obj.checkChannelOutput(iChan)
                fileNames = cellfun(@(x)(dir([x filesep '*.tif'])),obj.outFilePaths_(1,iChan),'UniformOutput',false);
                fileNames = cellfun(@(x)(arrayfun(@(x)(x.name),x,'UniformOutput',false)),fileNames,'UniformOutput',false);
                nChan = numel(iChan);
                for j = 1:nChan
                    %Sort the files by the trailing numbers
                    fNums = cellfun(@(x)(str2double(...
                        x(max(regexp(x(1:end-4),'\D'))+1:end-4))),fileNames{j});
                    [~,iX] = sort(fNums);
                    fileNames{j} = fileNames{j}(iX);
                end
                nIm = cellfun(@(x)(length(x)),fileNames);
                if ~all(nIm == obj.owner_.nFrames_)
                    error('Incorrect number of images found in one or more channels!')
                end
            else
                error('Invalid channel numbers! Must be positive integers less than the number of image channels!')
            end
        end
        
        function fileNames = getInImageFileNames(obj,iChan)
            if obj.checkChanNum(iChan)
                nChan = numel(iChan);
                fileNames = cell(1,nChan);
                for j = 1:nChan
                    %First check for regular image inputs
                    fileNames{j} = imDir(obj.inFilePaths_{1,iChan(j)});
                    if isempty(fileNames{j})
                        %If none found, check for .mat image inputs
                        fileNames{j} = dir([obj.inFilePaths_{1,inFilePaths_iChan(j)} filesep '*.tif']);
                    end
                    fileNames{j} = arrayfun(@(x)(x.name),fileNames{j},'UniformOutput',false);
                    nIm = length(fileNames{j});
                    if nIm ~= obj.owner_.nFrames_
                        error(['Incorrect number of images found in channel ' num2str(iChan(j)) ' !'])
                    end
                end
            else
                error('Invalid channel numbers! Must be positive integers less than the number of image channels!')
            end
        end
        
        function setOutImagePath(obj,chanNum,imagePath)
            
            if ~obj.checkChanNum(chanNum)
                error('lccb:set:fatal','Invalid image channel number for image path!\n\n');
            end
            
            if ~iscell(imagePath)
                imagePath = {imagePath};
            end
            nChan = length(chanNum);
            if nChan ~= length(imagePath)
                error('lccb:set:fatal','You must specify a path for every channel!')
            end
            
            for j = 1:nChan
                if ~exist(imagePath{j},'dir')
                    error('lccb:set:fatal',...
                        ['The directory specified for channel ' ...
                        num2str(chanNum(j)) ' is invalid!'])
                else
                    obj.outFilePaths_{1,chanNum(j)} = imagePath{j};
                end
            end
            
            
        end
        
        
        
         function out_data = loadChannelOutput(obj,iChan,iFrame,varargin)
            % Input check
            ip =inputParser;
            ip.addRequired('iChan',@obj.checkChanNum);
            ip.addRequired('iFrame',@obj.checkFrameNum);
            
            outputList = {'MAX_st_res','nms','orienation_map','scaleMap'};
            ip.addParamValue('output',{},@(x) all(ismember(x,outputList)));
            
            ip.parse(iChan,iFrame,varargin{:})
            
            % Data loading
            Channel_FilesNames = obj.getInImageFileNames(iChan);
            filename_short_strs = uncommon_str_takeout(Channel_FilesNames{1});
            
            % this line in commandation for shortest version of filename
            filename_shortshort_strs = all_uncommon_str_takeout(Channel_FilesNames{1});
            
            Sub_Sample_Num = obj.funParams_.Sub_Sample_Num;
             Frames_to_Seg = 1:Sub_Sample_Num:obj.owner_.nFrames_;
             
             Frames_results_correspondence = im2col(repmat(Frames_to_Seg, [Sub_Sample_Num,1]),[1 1]);
             Frames_results_correspondence = Frames_results_correspondence(1:obj.owner_.nFrames_);
             
       
             % these loads are for old version of the naming system
             
             try
                out_data_all = load([obj.outFilePaths_{1,iChan},'/steerable_',filename_short_strs{iFrame},'.mat'], ...
                    'orienation_map', 'MAX_st_res','nms','scaleMap');
            catch
                try
                    out_data_all = load([obj.outFilePaths_{1,iChan},'/steerable_',filename_shortshort_strs{iFrame},'.mat'], ...
                                            'orienation_map', 'MAX_st_res','nms','scaleMap');
                catch
                    % when only on channel, one image, it will be last
                    % character of the image, so 'f'
                    out_data_all = load([obj.outFilePaths_{1,iChan},'/steerable_','f','.mat'], ...
                                   'orienation_map', 'MAX_st_res','nms','scaleMap');
                end
            end
            
            
            % if there is no output parameter
            if( isempty(ip.Results.output))
                out_data = out_data_all.current_seg;
            else
                % or according to user's defined output parameter
                switch ip.Results.output
                    case 'orienation_map'
                        out_data = out_data_all.orienation_map;
                    case 'MAX_st_res'
                        out_data = out_data_all.MAX_st_res;
                    case 'nms'
                        out_data = out_data_all.nms;
                    case 'scaleMap'
                        out_data = out_data_all.scaleMap;
                    otherwise
                        out_data = out_data_all.MAX_st_res;
                end
            end
        end
        
        
        function h = draw(obj,iChan,varargin)
            
            outputList = obj.getDrawableOutput();
            drawSteerableFilteringImage = any(strcmpi('SteerableFilteringImage',varargin));
            
            if drawSteerableFilteringImage
                % Input check
                ip =inputParser;
                ip.addRequired('iChan',@(x) ismember(x,1:numel(obj.owner_.channels_)));
                ip.addParamValue('output',[],@ischar);
                ip.KeepUnmatched = true;
                ip.parse(iChan,varargin{:})
                
                % Load average corrected image
                s = load(obj.outFilePaths_{2,iChan});
                tmpFields=fieldnames(s);
                data=s.(tmpFields{1});
                
                iOutput= find(cellfun(@(y) isequal(ip.Results.output,y),{outputList.var}));
                if ~isempty(outputList(iOutput).formatData),
                    data=outputList(iOutput).formatData(data);
                end
                
                try
                    assert(~isempty(obj.displayMethod_{iOutput,iChan}));
                catch ME
                    obj.displayMethod_{iOutput,iChan}=...
                        outputList(iOutput).defaultDisplayMethod(iChan);
                end
                
                % Delegate to the corresponding method
                tag = [obj.getName '_channel' num2str(iChan) '_output' num2str(iOutput)];
                drawArgs=reshape([fieldnames(ip.Unmatched) struct2cell(ip.Unmatched)]',...
                    2*numel(fieldnames(ip.Unmatched)),1);
                h=obj.displayMethod_{iOutput,iChan}.draw(data,tag,drawArgs{:});
            else
                h=draw@ImageProcessingProcess(obj,iChan,varargin{1},varargin{2:end});
            end
        end
    end
    
    
    methods (Static)
        function name =getName()
            name = 'Steerable filtering';
        end
        function h = GUI()
            h= @steerableFilteringProcessGUI;
        end
        
        function output = getDrawableOutput()
            output = ImageProcessingProcess.getDrawableOutput();
            %             output(2).name='Steerable Filtering images';
            %             output(2).var='SteerableFilteringImage';
            %             output(2).formatData=[];
            %             output(2).type='graph';
            %             output(2).defaultDisplayMethod=@(x)LineDisplay('Color',[0 0 0],...
            %                 'LineStyle','-','LineWidth',2,...
            %                 'XLabel','Frame Number','YLabel','SteerableFiltering');
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            
            funParams.ChannelIndex = 1:numel(owner.channels_);
            funParams.Levelsofsteerablefilters = 2;
            funParams.BaseSteerableFilterSigma = 1;
            funParams.ImageFlattenFlag = 2;
            % sub-sample number, since often VIF images are taken at a
            % lower sample rate than the other channel, so use this number
            % to save some time.
            funParams.Sub_Sample_Num = 1;
        end
    end
end