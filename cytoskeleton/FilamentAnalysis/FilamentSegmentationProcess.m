classdef FilamentSegmentationProcess < ImageProcessingProcess
    % A concrete class for steerable filtering
    %
    % Liya Ding, 06. 2012
    
    methods (Access = public)
        
        function obj = FilamentSegmentationProcess(owner,varargin)
            
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
                super_args{2} = FilamentSegmentationProcess.getName;
                super_args{3} = @filament_segmentation;
                if isempty(funParams)
                    funParams = FilamentSegmentationProcess.getDefaultParams(owner,outputDir);
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
            
            for j = 1 : nChan
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
            
            outputList = {'current_seg_orientation','tip_orientation',...
                'tip_int','tip_NMS', 'current_model','RGB_seg_orient_heat_map',''};
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
                out_data_all = load([obj.outFilePaths_{1,iChan},'/DataOutput/steerable_vote_',filename_short_strs{iFrame},'.mat'], ...
                    'current_seg','current_seg_orientation','tip_orientation','tip_int','tip_NMS','current_model','RGB_seg_orient_heat_map');
            catch
                try
                    out_data_all = load([obj.outFilePaths_{1,iChan},'/DataOutput/steerable_vote_',filename_shortshort_strs{iFrame},'.mat'], ...
                        'current_seg','current_seg_orientation','tip_orientation','tip_int','tip_NMS','current_model','RGB_seg_orient_heat_map');
                catch
                    % when only on channel, one image, it will be last
                    % character of the image, so 'f'
                    out_data_all = load([obj.outFilePaths_{1,iChan},'/DataOutput/steerable_vote_','f','.mat'], ...
                        'current_seg','current_seg_orientation','tip_orientation','tip_int','tip_NMS','current_model','RGB_seg_orient_heat_map');
             
                end
            end
            
            
            % if there is no output parameter
            if( isempty(ip.Results.output))
                out_data = out_data_all.current_seg;
            else
                % or according to user's defined output parameter
                switch ip.Results.output
                    case 'current_seg_orientation'
                        out_data = out_data_all.current_seg_orientation;
                    case 'tip_orientation'
                        out_data = out_data_all.tip_orientation;
                    case 'tip_int'
                        out_data = out_data_all.tip_int;
                    case 'tip_NMS'
                        out_data = out_data_all.tip_NMS;
                    case 'current_model'
                        out_data = out_data_all.current_model;
                    case 'RGB_seg_orient_heat_map'
                        out_data = out_data_all.RGB_seg_orient_heat_map;
                    otherwise
                        out_data = out_data_all.current_seg_orientation;
                end
            end
        end
        
        function h = draw(obj,iChan,varargin)
            
            outputList = obj.getDrawableOutput();
            drawFilamentSegmentationImage = any(strcmpi('FilamentSegmentationImage',varargin));
            
            if drawFilamentSegmentationImage
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
            name = 'Filament Segmentation';
        end
        function h = GUI()
            h= @filamentSegmentationProcessGUI;
        end
        
        function output = getDrawableOutput()            
            output(2).name='Segmentation with Orientaion';
            output(2).var='current_seg_orientation';
            output(2).formatData=@mat2gray;
            output(2).type='image';
            output(2).defaultDisplayMethod=@ImageDisplay;
            
            output(5).name='Orientation at Tip Only';
            output(5).var='tip_orientation';
            output(5).formatData=@mat2gray;
            output(5).type='image';
            output(5).defaultDisplayMethod=@ImageDisplay;
            
            output(3).name='Intensity at Tip Only';
            output(3).var='tip_int';
            output(3).formatData=@mat2gray;
            output(3).type='image';
            output(3).defaultDisplayMethod=@ImageDisplay;
            
            output(4).name='SF Results at Tip Only';
            output(4).var='tip_NMS';
            output(4).formatData=@mat2gray;
            output(4).type='image';
            output(4).defaultDisplayMethod=@ImageDisplay;
            
            output(1).name='Heat Display';
            output(1).var='RGB_seg_orient_heat_map';
            output(1).formatData=[];
            output(1).type='image';
            output(1).defaultDisplayMethod=@ImageDisplay;
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            % Set default channels, use all channels
            funParams.ChannelIndex = 1:numel(owner.channels_);
            
            % The parameter to set pace in local segmentation
            funParams.StPace_Size = 3;
            % The parameter to set patch size in local segmentation, for
            % the estimation of local threshold
            funParams.StPatch_Size  = 21;
            
            % The percentage as the lower bound of local thresholding 
            % local threshold has to be larger or equal to this percentage
            % of the global threshold
            funParams.st_lowerbound_localthresholding  = 90; % default 90%
            
                        
            % Same set of parameters for intensity based segmentation
            funParams.IntPace_Size = 3;
            funParams.IntPatch_Size  = 21;
            funParams.int_lowerbound_localthresholding  = 90; % default 90%
            
            
            % The way to combine segmentation results from steerable
            % filtering responce and from intensity, default is : only use
            % steerable filtering result and geobased nms
            funParams.Combine_Way = 'geo_based_GM';
            
            % Flag to set if cell mask is used, if 1, use
            % segmentation(refined) results, if 2, use the user define ROI
            % as in MD_ROI.tif in movieData folder, if 3, a combined
            % version of two channel, if 4 a direction sum of the two
            % channels, if 5 no limit
            funParams.Cell_Mask_ind = 1;
            
            % whole movie constrain index, 1 for completely, 2 for
            % half-half, 3 for none
            
            funParams.Whole_movie_ind = 2;
            
            funParams.Whole_movie_stat_cell = cell(1,max(funParams.ChannelIndex));
            
            % if existing, don't re run whole movie analysis, if set to 1,
            % rerun any way
            funParams.Rerun_WholeMovie = 0;
            
            % Flag to do VIF_outgrowth or not. This is an option made for
            % Gelfand lab
            funParams.VIF_Outgrowth_Flag = 0;
            
            % sub-sample number, since often VIF images are taken at a
            % lower sample rate than the other channel, so use this number
            % to save some time.
            funParams.Sub_Sample_Num = 1;
            
            % the classifiers trained, the mat file's name is saved here.
            funParams.F_classifier = cell(1,max(funParams.ChannelIndex));
 
            % the type of classifier
            funParams.Classifier_Type_ind=1;
            
            % the number of samples in training
            funParams.training_sample_number=30;
            
            % the length threshold
            funParams.LengthThreshold=4;
    
            % the curvature threshold
            funParams.CurvatureThreshold=0.1;
            
            % the default linear plane classifier with offset, Alpha
            funParams.CoefAlpha = 2;
            
            % iteration number
            funParams.IternationNumber=2;            
            
            
             % Canny thresholds in percentage
            funParams.CannyHigherThreshold=80;            
            funParams.CannyLowerThreshold=80;            
            
            
            % for code running preference
            % No jumping out figures disrupt other thing; change to 0 if want to see the figures
            funParams.nofiguredisruption = 1;
            % No saving the step figures or debug ones; change to 1 if want to save these figures
            funParams.savestepfigures = 0;
            % No displaying detailed messages: change to 1 if want to see all debugging messages
            funParams.showdetailmessages = 0;
            
            % the flag for if a channel has been specifically signed
            % setting, without being the same with all other selected of channels
            funParams.channel_specific=0;
        end
    end
end