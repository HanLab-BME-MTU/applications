classdef DisplacementFieldCorrectionProcess < DataProcessingProcess
    % Concrete class for a displacement field correction process
    %
    % Sebastien Besson, Aug 2011
    properties (SetAccess = protected)  
        tMapLimits_
    end
    
    methods
        function obj = DisplacementFieldCorrectionProcess(owner,varargin)
            
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
                super_args{2} = DisplacementFieldCorrectionProcess.getName;
                super_args{3} = @correctMovieDisplacementField;
                if isempty(funParams)
                    funParams=DisplacementFieldCorrectionProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            obj = obj@DataProcessingProcess(super_args{:});
        end
        
        function status = checkChannelOutput(obj,varargin)
            
            status = logical(exist(obj.outFilePaths_{1},'file'));
            
        end
        
        function varargout = loadChannelOutput(obj,varargin)
            outputList = {'displField','dMap','dMapRef','dMapUnshifted'};
            ip =inputParser;
            ip.addRequired('obj',@(x) isa(x,'DisplacementFieldCorrectionProcess'));
            ip.addOptional('iFrame',1:obj.owner_.nFrames_,...
                @(x) ismember(x,1:obj.owner_.nFrames_));
            ip.addParamValue('output',outputList{1},@(x) all(ismember(x,outputList)));
            ip.addParameter('useCache',true,@islogical);
            ip.addParameter('noStackRequired',false,@islogical) % 
            ip.parse(obj,varargin{:})
            iFrame = ip.Results.iFrame;
            
            % Data loading
            persistent dMapMap dMapMapRef lastFinishTime
            % Data loading
            output = ip.Results.output;
            noStackRequired = ip.Results.noStackRequired; % this variable is used to empty
            % tMapMap, anticipating it will be very large, creating
            % out-of-memory error. But if the full tMapMap is already
            % created, we'll use that.
            if ischar(output), output = {output}; end
            iOut = cellfun(@(x) strcmp(x,output),outputList);
%             s = load(obj.outFilePaths_{iOut},output{:});
%             s = cached.load(obj.outFilePaths_{iOut}, '-useCache', ip.Results.useCache, output{:});
            if ismember(output,outputList(1))
%                 iOut=1;
                s = cached.load(obj.outFilePaths_{iOut}, '-useCache', ip.Results.useCache, output{1});
                varargout{1}=s.(output{1})(iFrame);
            elseif ismember(output,outputList(2:4))
                iOut=2;
                if isempty(lastFinishTime)
                    lastFinishTime = clock; % assigning current time.. This will be definitely different from obj.finishTime_
                end
                if ~all(obj.finishTime_==lastFinishTime) % We have initialize maps if the process is updated or different
                    dMapMap = [];
                    dMapMapRef = [];
                end
                if (strcmp(output,'dMap') || strcmp(output,'dMapUnshifted')) && (isempty(dMapMap) ...
                        || size(dMapMap,3)<iFrame(end) ...
                        || (size(dMapMap,3)>=iFrame && ~any(any(dMapMap(:,:,iFrame(end)))))) ...
                        || strcmp(output,'dMapRef') && (isempty(dMapMapRef) ...
                        || size(dMapMapRef,3)<iFrame(end) ...
                        || (size(dMapMapRef,3)>=iFrame && ~any(any(dMapMapRef(:,:,iFrame(end))))))
                    try
                        s = load(obj.outFilePaths_{iOut});
                        dMapObj = s.(outputList{2});
                        fString = ['%0' num2str(floor(log10(obj.owner_.nFrames_))+1) '.f'];
                        numStr = @(frame) num2str(frame,fString);
                        outputDir = fullfile(obj.funParams_.OutputDirectory,'displMaps');
                        outFileDMap = @(frame) [outputDir filesep 'displMap' numStr(frame) '.mat'];
                        if ~isstruct(dMapObj)
                            % This means that the data is is stored in an old
                            % way. (cell array). 
                            disp('The displacement map is stored in a huge cell array format with no compression.')
                            disp('Reformating and re-saving the individual maps using compression..')
                            s = load(obj.outFilePaths_{iOut});
                            mkdir(outputDir);
                            for ii = obj.owner_.nFrames_:-1:1
                                cur_dMap = s.dMap{ii};
                                cur_dMapX = s.dMapX{ii};
                                cur_dMapY = s.dMapY{ii};
                                dMapMap(:,:,ii) = cur_dMap;
                                save(outFileDMap(ii),'cur_dMap','cur_dMapX','cur_dMapY'); % I removed v7.3 option to save the space,'-v7.3');
                                progressText((obj.owner_.nFrames_-ii)/obj.owner_.nFrames_,'One-time displacement map saving') % Update text
                            end
                            dMap.outFileTMap = @(frame) [outputDir filesep 'displMap' numStr(frame) '.mat'];
                            dMap.eachTMapName = 'cur_dMap';
                            dMap.outputDir = outputDir;
                            dMapX=dMap; dMapX.eachTMapName = 'cur_dMapX';
                            dMapY=dMap; dMapY.eachTMapName = 'cur_dMapY';
                            save(obj.outFilePaths_{iOut},'dMap','dMapX','dMapY'); % need to be updated for faster loading. SH 20141106
                            lastFinishTime = obj.finishTime_;
                        elseif isfield(dMapObj, 'eachDMapName')
%                             for ii=obj.owner_.nFrames_:-1:1
%                                 cur_dMapObj = load(outFileDMap(ii),dMapObj.eachDMapName);
%                                 dMapMap(:,:,ii) = cur_dMapObj.cur_dMap;
%                                 progressText((obj.owner_.nFrames_-ii)/obj.owner_.nFrames_,'One-time displacement map loading') % Update text
%                             end
                            for ii=iFrame
                                cur_dMapObj = load(outFileDMap(ii),dMapObj.eachDMapName);
                                curMap = cur_dMapObj.cur_dMap;
                                if ~noStackRequired
                                    dMapMap(:,:,ii) = curMap;
                                end
                            end
                            lastFinishTime = obj.finishTime_;
                        else % very new format
                            displFieldObj = cached.load(dMapObj.displFieldPath, '-useCache', ip.Results.useCache, 'displField');
                            displField = displFieldObj.displField;
                            [dMapIn, ~, ~, cropInfo] = generateHeatmapShifted(displField(iFrame),displField(iFrame),0);
                            pp=numel(iFrame)+1;
                            for ii=fliplr(iFrame)
                                pp=pp-1;
                                curMap = zeros(dMapObj.firstMaskSize);
                                curMap(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3)) = dMapIn{pp};
                                if ~noStackRequired
                                    dMapMap(:,:,ii) = curMap;
                                end
%                                     progressText((obj.owner_.nFrames_-ii)/obj.owner_.nFrames_,'One-time traction map loading') % Update text
                            end
                            lastFinishTime = obj.finishTime_;
                        end
                    catch
                        tfmPack = obj.owner_.packages_{obj.getPackageIndex};
                        dMapObj.displFieldPath = [tfmPack.outputDirectory_ filesep 'correctedDisplacementField' filesep 'displField.mat'];
                        displField = load(dMapObj.displFieldPath,'displField'); displField=displField.displField;
                        
                        SDCproc = tfmPack.processes_{1};
                        if isempty(SDCproc)
                            dMapObj.firstMaskSize = size(obj.owner_.channels_(1).loadImage(1));
                        else
                            dMapObj.firstMaskSize = size(SDCproc.loadChannelOutput(1,1));
                        end
                        
%                         if exist('dMapMap','var')
%                             clear dMapMap
%                         end
                        
                        if ismember(output,outputList([2 4])) 
                            [dMapIn, ~, ~, cropInfo] = generateHeatmapShifted(displField(iFrame),displField(iFrame),0);
                            pp=numel(iFrame)+1;
                            for ii=fliplr(iFrame)
                                pp=pp-1;
                                curMap = zeros(dMapObj.firstMaskSize);
                                curMap(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3)) = dMapIn{pp};
                                if ~noStackRequired
                                    dMapMap(:,:,ii) = curMap;
                                end
%                                 progressText((obj.owner_.nFrames_-ii)/obj.owner_.nFrames_,'One-time displacement map loading') % Update text
                            end
                        elseif strcmp(output,'dMapRef')
                            [dMapIn, ~, ~, cropInfo] = generateHeatmapShifted(displField(iFrame),0,0);
                            iTFMPack = obj.owner_.getPackageIndex('TFMPackage');
                            tfmPackageHere=obj.owner_.packages_{iTFMPack}; iSDCProc=1;
                            SDCProc=tfmPackageHere.processes_{iSDCProc};

                            pp=numel(iFrame)+1;
                            for ii=fliplr(iFrame)
                                pp=pp-1;
                                curMapRef = zeros(size(SDCProc.loadOutImage(1,1)));
                                curMapRef(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3)) = dMapIn{pp};
                                if ~noStackRequired
                                    dMapMapRef(:,:,ii) = curMapRef;
                                end
%                                 progressText((obj.owner_.nFrames_-ii)/obj.owner_.nFrames_,'One-time displacement map loading') % Update text
                            end
                        end
                        lastFinishTime = obj.finishTime_;
                    end
                end
                if strcmp(output,'dMap')
                    if noStackRequired
                        varargout{1} = curMap;
                    else
                        varargout{1}=dMapMap(:,:,iFrame);
                    end
                elseif strcmp(output,'dMapRef')
                    if noStackRequired
                        varargout{1} = curMapRef;
                    else
                        varargout{1}=dMapMapRef(:,:,iFrame);
                    end
                else %This is for unshifted (in the size of raw channels)
                    sampleRawChanImg = obj.owner_.channels_(1).loadImage(1);
                    ref_obj = imref2d(size(sampleRawChanImg));
                    iTFMPack = obj.owner_.getPackageIndex('TFMPackage');
                    tfmPackageHere=obj.owner_.packages_{iTFMPack}; iSDCProc=1;
                    SDCProc=tfmPackageHere.processes_{iSDCProc};

                    if ~isempty(SDCProc)
                        try
                            iBeadChan=SDCProc.funParams_.iBeadChannel;
                        catch
                            iBeadChan=1;
                        end
                        s = load(SDCProc.outFilePaths_{3,iBeadChan},'T');
                        T = s.T;
                        if length(iFrame)>1
                            if noStackRequired
                                curMap2 = zeros(size(curMap));
                            else
                                curMap2 = zeros(size(dMapMap(:,:,iFrame)));
                            end
                            for ii=fliplr(iFrame)
                                Tr = affine2d([1 0 0; 0 1 0; (T(ii, :)) 1]);
                                invTr = invert(Tr);
                                if noStackRequired
                                    curMap2 = imwarp(curMap,invTr,'OutputView',ref_obj);
                                else
                                    curMap2(:,:,ii) = imwarp(dMapMap(:,:,ii),invTr,'OutputView',ref_obj);
                                end
                            end
                            varargout{1} = curMap2;
                        else
                            Tr = affine2d([1 0 0; 0 1 0; fliplr(T(iFrame, :)) 1]);
                            invTr = invert(Tr);
                            if noStackRequired
                                if ~exist('curMap','var') % this is because tMapMap was there.
                                    curMap = dMapMap(:,:,iFrame);
                                end
                                varargout{1} = imwarp(curMap,invTr,'OutputView',ref_obj);
                            else
                                varargout{1} = imwarp(dMapMap(:,:,iFrame),invTr,'OutputView',ref_obj);
                            end
                        end
                    else
                        if noStackRequired
                            varargout{1} = curMap;
                        else
                            varargout{1}=dMapMap(:,:,iFrame);
                        end
                    end
                end
            end
            
%             if numel(iFrame)>1,
%                 for i=1:numel(output),
%                     varargout{i}=s.(output{i});
%                 end
%             else
%                 for i=1:numel(output),
%                     varargout{i}=s.(output{i})(iFrame);
%                 end
%             end
        end

        function setTractionMapLimits(obj,tMapLimits)
            obj.tMapLimits_ = tMapLimits;
        end
        
        function h=draw(obj,iFrame,varargin)
            % Function to draw process output (template method)
            
            if ~ismember('getDrawableOutput',methods(obj)), h=[]; return; end
            outputList = obj.getDrawableOutput();

            rendertMap = (contains(varargin(3),'dMap'));
            if rendertMap
                % Input parser
                ip = inputParser;
                ip.addRequired('obj',@(x) isa(x,'Process'));
%                 ip.addRequired('iChan',@isnumeric);
                ip.addRequired('iFrame',@isnumeric);
                ip.addParamValue('output',outputList(2).var,...
                    @(x) any(cellfun(@(y) isequal(x,y),{outputList.var})));
                ip.KeepUnmatched = true;
                ip.parse(obj,varargin{1:end})
                iFrame=ip.Results.iFrame;
                % recognize how big the movie is and determine if
                % noStackRequired is used or not. I will say if the movie
                % is larger than 512x512x300, we will call noStackRequired.
                nFrames = obj.owner_.nFrames_;
                movie3DSize = obj.owner_.imSize_(1)*obj.owner_.imSize_(2)*nFrames;
                thres3DSize = 512*512*299;
                if movie3DSize > thres3DSize
                    data=obj.loadChannelOutput('iFrame',iFrame,'output',ip.Results.output,'noStackRequired',true);
                else
                    data=obj.loadChannelOutput('iFrame',iFrame,'output',ip.Results.output);
                end
                if iscell(data), data = data{1}; end
            else                
                ip = inputParser;
                ip.addRequired('obj',@(x) isa(x,'Process'));
                ip.addRequired('iFrame',@isnumeric);
                ip.addParamValue('output',outputList(1).var,@(x) any(cellfun(@(y) isequal(x,y),{outputList.var})));
                ip.KeepUnmatched = true;
                ip.parse(obj,iFrame,varargin{:})
                data=obj.loadChannelOutput(iFrame,'output',ip.Results.output);
            end
            iOutput= find(cellfun(@(y) isequal(ip.Results.output,y),{outputList.var}));
            if ~isempty(outputList(iOutput).formatData),
                data=outputList(iOutput).formatData(data);
            end
            try
                assert(~isempty(obj.displayMethod_{iOutput,1}));
            catch ME
                obj.displayMethod_{iOutput,1}=...
                    outputList(iOutput).defaultDisplayMethod();
            end

            % Delegate to the corresponding method
            tag = ['process' num2str(obj.getIndex) '_output' num2str(iOutput)];
            drawArgs=reshape([fieldnames(ip.Unmatched) struct2cell(ip.Unmatched)]',...
                2*numel(fieldnames(ip.Unmatched)),1);
            h=obj.displayMethod_{iOutput}.draw(data,tag,drawArgs{:});
        end                
        function output = getDrawableOutput(obj)
            output(1).name='Displacement field';
            output(1).var='displField';
            output(1).formatData=@(x) [x.pos x.vec];
            output(1).type='movieOverlay';
            output(1).defaultDisplayMethod=@(x) VectorFieldDisplay('Color','b');
            
            output(2).name='Displacement map (in cell coord)';
            output(2).var='dMap';
            output(2).formatData=[];
            output(2).type='image';
            output(2).defaultDisplayMethod=@(x)ImageDisplay('Colormap','jet',...
                'Colorbar','on','Units',obj.getUnits,'CLim',obj.tMapLimits_);

            output(3).name='Displacement map (in ref coord)';
            output(3).var='dMapRef';
            output(3).formatData=[];
            output(3).type='image';
            output(3).defaultDisplayMethod=@(x) ImageDisplay('Colormap','jet','Colorbar','on','Units',obj.getUnits,'CLim',obj.tMapLimits_);

            output(4).name='Displacement map unshifted';
            output(4).var='dMapUnshifted';
            output(4).formatData=[];
            output(4).type='image';
            output(4).defaultDisplayMethod=@(x) ImageDisplay('Colormap','jet',...
                'Colorbar','on','Units',obj.getUnits,'CLim',obj.tMapLimits_);
        end
    end
    
    methods (Static)
        function name =getName()
            name = 'Displacement Field Correction';
        end
        function h = GUI()
            h= @displacementFieldCorrectionProcessGUI;
        end
                
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.OutputDirectory = [outputDir  filesep 'correctedDisplacementField'];
            funParams.doRotReg=0;
            funParams.outlierThreshold = 2;
            funParams.fillVectors=0;
            funParams.angDiffThreshold=1;
            funParams.magDiffThreshold=2;
        end
        function units = getUnits(varargin)
            units = 'Displacement (Pix)';
        end
    end
end