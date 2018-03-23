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
            ip.parse(obj,varargin{:})
            iFrame = ip.Results.iFrame;
%             iOut = ip.Results.iOut;
            
            persistent tMapMap lastFinishTime
            % Data loading
            output = ip.Results.output;
            if ischar(output), output = {output}; end
            if ismember(output,outputList(1:4))
                iOut=1;
                if ismember(output,outputList(1:3)) 
                    s = cached.load(obj.outFilePaths_{iOut}, '-useCache', ip.Results.useCache, output{1});
                    varargout{1}=s.(output{1})(iFrame);
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
                if isempty(tMapMap) || ~all(obj.finishTime_==lastFinishTime)
                    s = load(obj.outFilePaths_{iOut},outputList{5}); %This will need to be changed if one really wants to see tMapX or tMapY
                    tMapObj = s.(outputList{5});
                    fString = ['%0' num2str(floor(log10(obj.owner_.nFrames_))+1) '.f'];
                    numStr = @(frame) num2str(frame,fString);
                    outputDir = fullfile(obj.funParams_.OutputDirectory,'tractionMaps');
                    mkdir(outputDir);
                    outFileTMap = @(frame) [outputDir filesep 'tractionMap' numStr(frame) '.mat'];
                    if ~isstruct(tMapObj)
                        % This means that the data is is stored in an old
                        % way. (cell array). 
                        disp('The traction map is stored in a huge cell array format with no compression.')
                        disp('Reformating and re-saving the individual maps using compression..')
                        s = load(obj.outFilePaths_{iOut});
                        for ii = obj.owner_.nFrames_:-1:1
                            cur_tMap = s.tMap{ii};
                            cur_tMapX = s.tMapX{ii};
                            cur_tMapY = s.tMapY{ii};
                            tMapMap(:,:,ii) = cur_tMap;
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
                    else
                        for ii=obj.owner_.nFrames_:-1:1
                            cur_tMapObj = load(outFileTMap(ii),tMapObj.eachTMapName);
                            tMapMap(:,:,ii) = cur_tMapObj.cur_tMap;
                            progressText((obj.owner_.nFrames_-ii)/obj.owner_.nFrames_,'One-time traction map loading') % Update text
                        end
                    end
                    lastFinishTime = obj.finishTime_;
                end
                if ismember(output,outputList(5:7)) 
                    varargout{1}=tMapMap(:,:,iFrame);
                else %This is for unshifted
                    curMap=tMapMap(:,:,iFrame);
                    ref_obj = imref2d(size(curMap));
                    iTFMPack = obj.owner_.getPackageIndex('TFMPackage');
                    tfmPackageHere=obj.owner_.packages_{iTFMPack}; iSDCProc=1;
                    SDCProc=tfmPackageHere.processes_{iSDCProc};
                    if ~isempty(SDCProc)
                        s = load(SDCProc.outFilePaths_{3,1},'T');
                        T = s.T;
                        Tr = affine2d([1 0 0; 0 1 0; fliplr(T(1, :)) 1]);
                        invTr = invert(Tr);
                        varargout{1} = imwarp(curMap,invTr,'OutputView',ref_obj);
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
                data=obj.loadChannelOutput('iFrame',iFrame,'output',ip.Results.output);
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
            funParams.lastToFirst=false;
            funParams.lcornerOptimal='optimal';
            funParams.tolx=0.2;
        end
        function units = getUnits(varargin)
            units = 'Traction (Pa)';
        end
    end
end