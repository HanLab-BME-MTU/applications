classdef DisplacementFieldCalculationProcess < ImageAnalysisProcess
    % Concrete class for a displacement field calculation process
    %
    % Sebastien Besson, Aug 2011
    properties (SetAccess = protected)  
        tMapLimits_
    end
    
    methods
        function obj = DisplacementFieldCalculationProcess(owner,varargin)
            
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
                super_args{2} = DisplacementFieldCalculationProcess.getName;
                super_args{3} = @calculateMovieDisplacementField;
                if isempty(funParams)
                    funParams=DisplacementFieldCalculationProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4}=funParams;
            end
            
            obj = obj@ImageAnalysisProcess(super_args{:});
        end
        
        function status = checkChannelOutput(obj,varargin)
            status = logical(exist(obj.outFilePaths_{1},'file'));
        end
        
        function sanityCheck(obj)
            sanityCheck@ImageAnalysisProcess(obj);
            channelIndex = obj.funParams_.ChannelIndex;
            assert(isscalar(channelIndex), 'lccb:Process:sanityCheck',...
                'A single bead channel must be selected to run this process');
            psfSigma = obj.owner_.channels_(channelIndex).psfSigma_;
            assert(~isempty(psfSigma), 'lccb:Process:sanityCheck',...
                ['The beads channel does not have a valid '...
                'standard deviation of the Gaussian point-spread function.']);
        end
        
        function varargout = loadChannelOutput(obj,varargin)
            
            outputList = {'displField','dMap'};
            ip =inputParser;
            ip.addRequired('obj',@(x) isa(x,'DisplacementFieldCalculationProcess'));
            ip.addOptional('iFrame',1:obj.owner_.nFrames_,...
                @(x) ismember(x,1:obj.owner_.nFrames_));
            ip.addParamValue('output',outputList{1},@(x) all(ismember(x,outputList)));
            ip.addParameter('useCache',true,@islogical);
            ip.parse(obj,varargin{:})
            iFrame = ip.Results.iFrame;
            
            persistent dMapMap lastFinishTime
            % Data loading
            output = ip.Results.output;
            if ischar(output), output = {output}; end
            iOut = cellfun(@(x) strcmp(x,output),outputList);
%             s = load(obj.outFilePaths_{iOut},output{:});
%             s = cached.load(obj.outFilePaths_{iOut}, '-useCache', ip.Results.useCache, output{:});
            if ismember(output,outputList(1))
%                 iOut=1;
                s = cached.load(obj.outFilePaths_{iOut}, '-useCache', ip.Results.useCache, output{1});
                varargout{1}=s.(output{1})(iFrame);
            elseif ismember(output,outputList(2))
%                 iOut=2;
                if isempty(lastFinishTime)
                    lastFinishTime = clock; % assigning current time.. This will be definitely different from obj.finishTime_
                end
                if isempty(dMapMap) || ~all(obj.finishTime_==lastFinishTime)
                    s = load(obj.outFilePaths_{iOut},output{1});
                    dMapObj = s.(output{1});
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
                        for ii=obj.owner_.nFrames_:-1:1
                            cur_dMapObj = load(outFileDMap(ii),dMapObj.eachDMapName);
                            dMapMap(:,:,ii) = cur_dMapObj.cur_dMap;
                            progressText((obj.owner_.nFrames_-ii)/obj.owner_.nFrames_,'One-time displacement map loading') % Update text
                        end
                        lastFinishTime = obj.finishTime_;
                    else % very new format
                        displField = load(dMapObj.displFieldPath,'displField'); displField=displField.displField;
                        [dMapIn, ~, ~, cropInfo] = generateHeatmapShifted(displField,displField,0);
                        for ii=obj.owner_.nFrames_:-1:1
                            dMapMap(:,:,ii) = zeros(dMapObj.firstMaskSize);
                            dMapMap(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3),ii) = dMapIn{ii};
                            progressText((obj.owner_.nFrames_-ii)/obj.owner_.nFrames_,'One-time traction map loading') % Update text
                        end
                    end
                end
                varargout{1}=dMapMap(:,:,iFrame);
            end
            
%             varargout = cell(numel(output),1);
%             if numel(iFrame)>1
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

            rendertMap = any(strcmpi('dMap',varargin));
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
                data=obj.loadChannelOutput('iFrame',iFrame,'output',ip.Results.output);
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
                data=outputList(iOutput,1).formatData(data);
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
            output(1).defaultDisplayMethod=@(x) VectorFieldDisplay('Color','r');

            output(2).name='Displacement map';
            output(2).var='dMap';
            output(2).formatData=[];
            output(2).type='image';
            output(2).defaultDisplayMethod=@(x) ImageDisplay('Colormap','jet','Colorbar','on','Units',obj.getUnits,'CLim',obj.tMapLimits_);
        end
        
    end
    methods (Static)
        function name =getName()
            name = 'Displacement Field Calculation';
        end
        function h = GUI()
            h= @displacementFieldCalculationProcessGUI;
        end
        
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex = 1;
            funParams.OutputDirectory = [outputDir  filesep 'displacementField'];
            funParams.referenceFramePath='';
            funParams.alpha=.05;
            funParams.minCorLength = 21;
            funParams.maxFlowSpeed =20;
            funParams.highRes = true;
            funParams.mode = 'fast';
            funParams.useGrid = false;
            funParams.noFlowOutwardOnBorder = true;
            funParams.lastToFirst=false;
            funParams.addNonLocMaxBeads = false;
            funParams.trackSuccessively = false;
            funParams.sigCrit = 0.5;
        end
        function units = getUnits(varargin)
            units = 'Displacement (Pix)';
        end
    end
end