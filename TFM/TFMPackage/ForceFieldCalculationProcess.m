classdef ForceFieldCalculationProcess < DataProcessingProcess
    % Concrete process for calculating a force field
    %
    % Sebastien Besson, Aug 2011
    properties (SetAccess = protected)  
        tMapLimits_
        dELimits_
        distBeadMapLimits_
    end
    
    methods
        function obj = ForceFieldCalculationProcess(owner,varargin)
            
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
                super_args{2} = ForceFieldCalculationProcess.getName;
                super_args{3} = @calculateMovieForceField;
                if isempty(funParams)
                    funParams=ForceFieldCalculationProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            
            obj = obj@DataProcessingProcess(super_args{:});
            
        end
        
        function status = checkChannelOutput(obj,varargin)
            
            status = logical(exist(obj.outFilePaths_{1},'file'));
            
        end
        
        function varargout = loadChannelOutput(obj,varargin)
            
            outputList = {'forceField','forceFieldShifted','forceFieldShiftedColor','forceFieldUnshifted',...
                'tMap','tMapX','tMapY','tMapUnshifted'};
            ip =inputParser;
            ip.addRequired('obj',@(x) isa(x,'ForceFieldCalculationProcess'));
            ip.addOptional('iFrame',1:obj.owner_.nFrames_,@(x) all(obj.checkFrameNum(x)));
%             ip.addOptional('iOut',1,@isnumeric);
%             ip.addOptional('iFrame',1:obj.owner_.nFrames_,@(x) ismember(x,1:obj.owner_.nFrames_));
%             ip.addParamValue('output',outputList{1},@(x) all(ismember(x,outputList)));
            ip.addParameter('output',outputList{1},@(x) all(ismember(x,outputList)));
            ip.addParameter('useCache',true,@islogical);
            ip.addParameter('noStackRequired',false,@islogical) % 
            ip.parse(obj,varargin{:})
            iFrame = ip.Results.iFrame;
            noStackRequired = ip.Results.noStackRequired; % this variable is used to empty
            % tMapMap, anticipating it will be very large, creating
            % out-of-memory error. But if the full tMapMap is already
            % created, we'll use that.
%             iOut = ip.Results.iOut;
            
            persistent tMapMap tMapMapX tMapMapY lastFinishTime
            
            % Data loading
            output = ip.Results.output;
            if ischar(output), output = {output}; end
            if ismember(output,outputList(1:4))
                iOut=1;
                if ismember(output,outputList(1:3)) 
                    if ismember(output,outputList(1:2)) 
                        s = cached.load(obj.outFilePaths_{iOut}, '-useCache', ip.Results.useCache, output{1});
                        varargout{1}=s.(output{1})(iFrame);
                    else
                        s = cached.load(obj.outFilePaths_{iOut}, '-useCache', ip.Results.useCache, outputList{2});
                        varargout{1}=s.(outputList{2})(iFrame);
                    end
                else %This is for unshifted
                    s = cached.load(obj.outFilePaths_{iOut}, '-useCache', ip.Results.useCache, outputList{2});
                    curField = s.(outputList{2})(iFrame);
                    numPos = size(curField.pos,1);
                    iTFMPack = obj.owner_.getPackageIndex('TFMPackage');
                    tfmPackageHere=obj.owner_.packages_{iTFMPack}; iSDCProc=1;
                    SDCProc=tfmPackageHere.processes_{iSDCProc};
                    if ~isempty(SDCProc)
                        s = load(SDCProc.outFilePaths_{3,1},'T');
                        T = s.T;
                        Tr = affine2d([1 0 0; 0 1 0; fliplr(T(1, :)) 1]);
                        invTr = invert(Tr);
                        for ii=1:numPos
                            curField.pos(ii,1)=curField.pos(ii,1)+invTr.T(3,1);
                            curField.pos(ii,2)=curField.pos(ii,2)+invTr.T(3,2);
                        end
                        varargout{1}=curField;
                    else
                        varargout{1}=curField;
                    end
                end
%             elseif strcmp(output,outputList{5})
%                 [OutputDirectory,tMapFolder] = fileparts(obj.outFilePaths_{2});
%                 % Set up the output directories
%                 outputDir = fullfile(OutputDirectory,tMapFolder);
%                 outFileTMap=@(frame) [outputDir filesep 'tractionMap' numStr(frame) '.mat'];
            elseif ismember(output,outputList(5:8))
                iOut=2;
                if isempty(lastFinishTime)
                    lastFinishTime = clock; % assigning current time.. This will be definitely different from obj.finishTime_
                end
                if ~all(obj.finishTime_==lastFinishTime) % We have initialize maps if the process is updated or different
                    tMapMap = [];
                end
                if isempty(tMapMap) || size(tMapMap,3)<length(iFrame) ... %length(iFrame)==1 || 
                        || size(tMapMap,3)<iFrame(end) ...
                        || (size(tMapMap,3)>=iFrame(end) && ~any(any(tMapMap(:,:,iFrame(end))))) 
                    try
                        s = load(obj.outFilePaths_{iOut}); %This will need to be changed if one really wants to see tMapX or tMapY
                        fString = ['%0' num2str(floor(log10(obj.owner_.nFrames_))+1) '.f'];
                        numStr = @(frame) num2str(frame,fString);
                        outputDir = fullfile(obj.funParams_.OutputDirectory,'tractionMaps');
                        outFileTMap = @(frame) [outputDir filesep 'tractionMap' numStr(frame) '.mat'];
                        if isfield(s,outputList{5}) || isfield(s,'tMapObj')
                            tMapObj = s.(outputList{5});
                            if ~isstruct(tMapObj)
                                mkdir(outputDir);
                                % This means that the data is is stored in an old
                                % way. (cell array). 
                                disp('The traction map is stored in a huge cell array format with no compression.')
                                disp('Reformating and re-saving the individual maps using compression..')
                                s = load(obj.outFilePaths_{iOut});
                                for ii = obj.owner_.nFrames_:-1:1
                                    cur_tMap = s.tMap{ii};
                                    %backward compatibility
                                    if isfield(s,'tMapX')
                                        cur_tMapX = s.tMapX{ii};
                                        cur_tMapY = s.tMapY{ii};
                                    else
                                        cur_tMapX = [];
                                        cur_tMapY = []; %Later this can be changed to the code that actually generates tMapX and Y.
                                    end   
                                    if ~noStackRequired
                                        tMapMap(:,:,ii) = cur_tMap;
                                    else
                                        tMapMap = [];
                                    end
                                    if ii==1 && strcmpi(obj.funParams_.method,'FastBEM')
                                        try
                                            cur_fCfdMap = s.fCfdMap;
                                        catch
                                            cur_fCfdMap = s.fCfdMap{ii};
                                        end
                                        save(outFileTMap(ii),'cur_tMap','cur_tMapX','cur_tMapY','cur_fCfdMap'); % I removed v7.3 option to save the space,'-v7.3');
                                    else
                                        save(outFileTMap(ii),'cur_tMap','cur_tMapX','cur_tMapY'); % I removed v7.3 option to save the space,'-v7.3');
                                    end
                                    progressText((obj.owner_.nFrames_-ii)/obj.owner_.nFrames_,'One-time traction map saving') % Update text
                                end
                                tMap.outFileTMap = @(frame) [outputDir filesep 'tractionMap' numStr(frame) '.mat'];
                                tMap.eachTMapName = 'cur_tMap';
                                tMap.outputDir = outputDir;
                                tMapX=tMap; tMapX.eachTMapName = 'cur_tMapX';
                                tMapY=tMap; tMapY.eachTMapName = 'cur_tMapY';
                                if strcmpi(obj.funParams_.method,'FastBEM')
                                    fCfdMap=tMap; fCfdMap.eachTMapName = 'cur_fCfdMap';
                                    save(obj.outFilePaths_{iOut},'tMap','tMapX','tMapY','fCfdMap'); % need to be updated for faster loading. SH 20141106
                                else
                                    save(obj.outFilePaths_{iOut},'tMap','tMapX','tMapY'); % need to be updated for faster loading. SH 20141106
                                end
                            elseif isfield(tMapObj,'eachTMapName')
%                                 for ii=obj.owner_.nFrames_:-1:1
%                                     cur_tMapObj = load(outFileTMap(ii),tMapObj.eachTMapName);
%                                     tMapMap(:,:,ii) = cur_tMapObj.cur_tMap;
%                                     progressText((obj.owner_.nFrames_-ii)/obj.owner_.nFrames_,'One-time traction map loading') % Update text
%                                 end
                                for ii=iFrame
                                    cur_tMapObj = load(outFileTMap(ii),tMapObj.eachTMapName);
                                    curMap = cur_tMapObj.cur_tMap;
                                    if ~noStackRequired
                                        tMapMap(:,:,ii) = curMap;
                                    end
                                end
                            else % very new format
                                forceFieldObj = cached.load(tMapObj.forceFieldPath, '-useCache', ip.Results.useCache, outputList{1});
                                forceField=forceFieldObj.forceField;
                                displFieldObj = cached.load(tMapObj.displFieldPath, '-useCache', ip.Results.useCache, 'displField');
                                displField = displFieldObj.displField;
                                [tMapIn, ~, ~, cropInfo, tMapXIn, tMapYIn] = generateHeatmapShifted(forceField(iFrame),displField(iFrame),0); %,iFrame);
                                pp=numel(iFrame)+1;
                                for ii=fliplr(iFrame)
                                    pp=pp-1;
                                    curMap = zeros(tMapObj.firstMaskSize);
                                    curMap(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3)) = tMapIn{pp};
                                    if ~noStackRequired
                                        tMapMap(:,:,ii) = curMap;
                                    else
                                        tMapMap = [];
                                    end
                                    if ismember(output,outputList(6:7))
                                        tMapMapX(:,:,ii) = zeros(tMapObj.firstMaskSize);
                                        tMapMapX(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3),ii) = tMapXIn{pp};
                                        tMapMapY(:,:,ii) = zeros(tMapObj.firstMaskSize);
                                        tMapMapY(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3),ii) = tMapYIn{pp};
                                    end
                                    progressText((obj.owner_.nFrames_-ii)/obj.owner_.nFrames_,'One-time traction map loading') % Update text
                                end
                                
                            end
                            lastFinishTime = obj.finishTime_;
                        end
                    catch
                        tfmPack = obj.owner_.packages_{obj.getPackageIndex};
                        tMapObj.forceFieldPath = [tfmPack.outputDirectory_ filesep 'forceField' filesep 'forceField.mat'];
                        tMapObj.displFieldPath = [tfmPack.outputDirectory_ filesep 'correctedDisplacementField' filesep 'displField.mat'];
                        forceField = load(tMapObj.forceFieldPath,'forceField'); forceField=forceField.forceField;
                        displField = load(tMapObj.displFieldPath,'displField'); displField=displField.displField;
                        
                        SDCproc = tfmPack.processes_{1};
                        if isempty(SDCproc)
                            tMapObj.firstMaskSize = size(obj.owner_.channels_(1).loadImage(1));
                        else
                            tMapObj.firstMaskSize = size(SDCproc.loadChannelOutput(1,1));
                        end
                        
                        [tMapIn, ~, ~, cropInfo] = generateHeatmapShifted(forceField(iFrame),displField(iFrame),0);
                        pp=numel(iFrame)+1;
                        for ii=fliplr(iFrame)
                            pp=pp-1;
                            curMap = zeros(tMapObj.firstMaskSize);
                            curMap(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3)) = tMapIn{pp};
                            if ~noStackRequired
                                tMapMap(:,:,ii) = curMap;
                            end
                            progressText((obj.owner_.nFrames_-ii)/obj.owner_.nFrames_,'One-time traction map loading') % Update text
                        end
                        lastFinishTime = obj.finishTime_;
                    end
                end
                if ismember(output,outputList(5:7)) 
                    if strcmp(output,outputList(5))
                        if noStackRequired
                            varargout{1} = curMap;
                        else
                            varargout{1}=tMapMap(:,:,iFrame);
                        end
                    elseif strcmp(output,outputList(6))
                        varargout{1}=tMapMapX(:,:,iFrame);
                    elseif strcmp(output,outputList(7))
                        varargout{1}=tMapMapY(:,:,iFrame);
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
                                curMap2 = zeros(size(tMapMap(:,:,iFrame)));
                            end
                            for ii=fliplr(iFrame)
                                Tr = affine2d([1 0 0; 0 1 0; (T(ii, :)) 1]);
                                invTr = invert(Tr);
                                if noStackRequired
                                    curMap2 = imwarp(curMap,invTr,'OutputView',ref_obj);
                                else
                                    curMap2(:,:,ii) = imwarp(tMapMap(:,:,ii),invTr,'OutputView',ref_obj);
                                end
                            end
                            varargout{1} = curMap2;
                        else
                            Tr = affine2d([1 0 0; 0 1 0; fliplr(T(iFrame, :)) 1]);
                            invTr = invert(Tr);
                            if noStackRequired
                                if ~exist('curMap','var') % this is because tMapMap was there.
                                    curMap = tMapMap(:,:,iFrame);
                                end
                                varargout{1} = imwarp(curMap,invTr,'OutputView',ref_obj);
                            else
                                varargout{1} = imwarp(tMapMap(:,:,iFrame),invTr,'OutputView',ref_obj);
                            end
                        end
                    else
                        varargout{1}=tMapMap(:,:,iFrame);
                    end
                end
%                 s = load(obj.outFilePaths_{iOut},output{1});
            end
            
% %             if numel(iFrame)>1,
%             varargout{1}=s.(output{1})(iFrame);
% %             else
% %                 varargout{1}=s.(output{1});
% %             end
        end
                
        function h=draw(obj,varargin)
            % Function to draw process output
            
            outputList = obj.getDrawableOutput();
            drawLcurve = any(strcmpi('lcurve',varargin));
            rendertMap = any(strncmpi('tMap',varargin,4) | ...
                strncmpi('dErrMap',varargin,4) | strncmpi('distBeadMap',varargin,4) );
            if drawLcurve %Lcurve
                ip = inputParser;
                ip.addRequired('obj',@(x) isa(x,'Process'));
                ip.addParameter('output',outputList(1).var,...
                    @(x) any(cellfun(@(y) isequal(x,y),{outputList.var})));
                ip.KeepUnmatched = true;
                ip.parse(obj,varargin{:})
                data=obj.outFilePaths_{4,1};
            elseif rendertMap % forceMap
                % Input parser
                ip = inputParser;
                ip.addRequired('obj',@(x) isa(x,'Process'));
                ip.addRequired('iChan',@isnumeric);
                ip.addRequired('iFrame',@isnumeric);
                ip.addParameter('output',outputList(2).var,...
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
            else % forcefield
                % Input parser
                ip = inputParser;
                ip.addRequired('obj',@(x) isa(x,'Process'));
                ip.addRequired('iFrame',@isnumeric);
                ip.addParameter('output',outputList(1).var,...
                    @(x) any(cellfun(@(y) isequal(x,y),{outputList.var})));
                ip.KeepUnmatched = true;
                ip.parse(obj,varargin{1},varargin{2:end})
                iFrame=ip.Results.iFrame;
                
                data=obj.loadChannelOutput('iFrame',iFrame,'output',ip.Results.output);
            end
            iOutput= find(cellfun(@(y) isequal(ip.Results.output,y),{outputList.var}));
            if ~isempty(outputList(iOutput).formatData)
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
        
        function setTractionMapLimits(obj,tMapLimits)
            obj.tMapLimits_ = tMapLimits;
        end
        function setDisplErrMapLimits(obj,dELimits)
            obj.dELimits_ = dELimits;
        end
        function setDistBeadMapLimits(obj,distBeadMapLimits)
            obj.distBeadMapLimits_ = distBeadMapLimits;
        end
        
        function output = getDrawableOutput(obj)
            output(1).name='Force  field';
            output(1).var='forceField';
            output(1).formatData=@(x) [x.pos x.vec(:,1)/nanmean((x.vec(:,1).^2+x.vec(:,2).^2).^0.5) x.vec(:,2)/nanmean((x.vec(:,1).^2+x.vec(:,2).^2).^0.5)];
            output(1).type='movieOverlay';
%             output(1).defaultDisplayMethod=@(x) VectorFieldDisplay('Color','r');
            output(1).defaultDisplayMethod=@(x) VectorFieldDisplay('Color',[75/255 0/255 130/255]);
            
            output(2).name='Traction map';
            output(2).var='tMap';
            output(2).formatData=[];
            output(2).type='image';
            output(2).defaultDisplayMethod=@(x)ImageDisplay('Colormap','jet',...
                'Colorbar','on','Units',obj.getUnits,'CLim',obj.tMapLimits_);

            output(3).name='Force field shifted';
            output(3).var='forceFieldShifted';
            output(3).formatData=@(x) [x.pos x.vec(:,1)/mean((x.vec(:,1).^2+x.vec(:,2).^2).^0.5) x.vec(:,2)/mean((x.vec(:,1).^2+x.vec(:,2).^2).^0.5)];
            output(3).type='movieOverlay';
            output(3).defaultDisplayMethod=@(x) VectorFieldDisplay('Color',[175/255 30/255 230/255]);

            output(4).name='Force field shifted (c)';
            output(4).var='forceFieldShiftedColor';
            output(4).formatData=@(x) [x.pos x.vec(:,1) x.vec(:,2)];
            output(4).type='movieOverlay';
            output(4).defaultDisplayMethod=@(x) VectorFieldDisplay('Colormap',jet,'Linewidth',1);
            
%             if ~strcmp(obj.funParams_.method,'FTTC')

            output(5).name='Lcurve';
            output(5).var='lcurve';
            output(5).formatData=[];
            output(5).type='movieGraph';
            output(5).defaultDisplayMethod=@FigFileDisplay;

            output(6).name='Traction map unshifted';
            output(6).var='tMapUnshifted';
            output(6).formatData=[];
            output(6).type='image';
            output(6).defaultDisplayMethod=@(x)ImageDisplay('Colormap','jet',...
                'Colorbar','on','Units',obj.getUnits,'CLim',obj.tMapLimits_);

            output(7).name='Force field SDC unshifted';
            output(7).var='forceFieldUnshifted';
            output(7).formatData=@(x) [x.pos x.vec(:,1)/mean((x.vec(:,1).^2+x.vec(:,2).^2).^0.5) x.vec(:,2)/mean((x.vec(:,1).^2+x.vec(:,2).^2).^0.5)];
            output(7).type='movieOverlay';
            output(7).defaultDisplayMethod=@(x) VectorFieldDisplay('Color',[125/255 50/255 210/255]);
                %% TODO -- Ensure outputs are generated and available for display
                % output(6).name='Prediction Err map';
                % output(6).var='dErrMap';
                % output(6).formatData=[];
                % output(6).type='image';
                % output(6).defaultDisplayMethod=@(x)ImageDisplay('Colormap','jet',...
                %     'Colorbar','on','Units',obj.getUnits,'CLim',obj.dELimits_);

                % output(7).name='Map of distance to bead';
                % output(7).var='distBeadMap';
                % output(7).formatData=[];
                % output(7).type='image';
                % output(7).defaultDisplayMethod=@(x)ImageDisplay('Colormap','jet',...
                %     'Colorbar','on','Units',obj.getUnits,'CLim',obj.distBeadMapLimits_);


%             end                
        end
        
        
    end
    methods (Static)
        function name =getName()
            name = 'Force Field Calculation';
        end
        function h = GUI()
            h= @forceFieldCalculationProcessGUI;
        end
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.OutputDirectory = [outputDir  filesep 'forceField'];
            funParams.YoungModulus = 8000;
            funParams.PoissonRatio = .5;
            funParams.method = 'FastBEM';
            funParams.meshPtsFwdSol = 4096;
            funParams.regParam=1e-4;
            funParams.solMethodBEM='1NormReg';
            funParams.basisClassTblPath='';
            funParams.LcurveFactor=10;
            funParams.thickness=34000;
            funParams.useLcurve=true;
            funParams.useLcurveEveryFrame=false;
            funParams.lastToFirst=false;
            funParams.lcornerOptimal='optimal';
            funParams.tolx=0.2;
        end
        function units = getUnits(varargin)
            units = 'Traction (Pa)';
        end
    end
end