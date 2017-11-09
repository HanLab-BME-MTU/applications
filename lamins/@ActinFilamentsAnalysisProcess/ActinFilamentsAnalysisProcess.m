classdef ActinFilamentsAnalysisProcess < DetectionProcess & NonSingularProcess
    %ActinFilamentsAnalysisProcess Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = ActinFilamentsAnalysisProcess(owner,varargin)
%             if(nargin < 1)
%                 % Allow empty creation
%                 return;
%             end
            ip = inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('funParams', ...
                ActinFilamentsAnalysisProcess.getDefaultParams(owner), ...
                @isstruct);
            ip.parse(owner,varargin{:});
            
            obj = obj@DetectionProcess(owner, ... 
                'ActinFilamentsAnalysisProcess', ... % name
                @actinFun, ... % funName
                ip.Results.funParams ...
                );
            
            obj.outFilePaths_ = ip.Results.funParams.outFilePaths;
        end
        function varargout = loadChannelOutput(obj,iChan,varargin)
            % Input check
            outputList = {'skeleton'};
            ip =inputParser;
            ip.StructExpand = true;
            ip.addRequired('iChan',@(x) isscalar(x) && obj.checkChanNum(x));
            ip.addOptional('iFrame',1,@(x) all(obj.checkFrameNum(x)));
            ip.addParamValue('iZ',1);
            ip.addParamValue('useCache',true,@islogical);
            ip.addParamValue('output',[], ...
                @(x) all(ismember(x,outputList)) ...
                );
            ip.parse(iChan,varargin{:})
            iFrame = ip.Results.iFrame;
            iZ = ip.Results.iZ;
            output = ip.Results.output;
            if(isempty(output))
                params = obj.getParameters();
                if(isfield(params,'defaultOutput'))
                    output = params.defaultOutput;
                else
                    output = 'nms_skel_bridged_skel_sc1';
                end
            end
            if ischar(output),output={output}; end
            
            varargout = cell(1,length(output));
            
            selectCell = cell(1,length(output));
            for o=1:length(output)
                [splits,tokens] = regexp(output{o},'{([0-9]+)}$|_sc([0-9]+)$','split','tokens');
                if(~isempty(tokens))
                    output{o} = splits{1};
                    selectCell{o} = str2double(tokens{1});
                end               

            end
            
            params = obj.getParameters();
            template = [obj.outFilePaths_{iChan} filesep 'actin_skeleton_c%02d_t%03d_z%03d.mat'];
            file = sprintf(template,iChan,ip.Results.iFrame,ip.Results.iZ);
            try
                s = cached.load(file, '-useCache', ip.Results.useCache, output{:});
            catch err
                varargout{1} = {zeros(obj.owner_.imSize_)};
                return;
            end
            
            switch(output{1})
                otherwise
                    varargout{1}=s.(output{1});
            end

            for o = 1:length(output)
                if(~isempty(selectCell{o}))
                    varargout{o} = varargout{o}{selectCell{o}};
                end
            end
        end
        function output = getDrawableOutput(obj)
            output(1).name='Meshwork';
            output(1).var='skeleton';
            output(1).formatData=@(x) x{2}.getEdgeXY;
            output(1).type='overlay';
            output(1).defaultDisplayMethod=@(x) LineDisplay('Marker','none',...
                'Color',[1 0 0]);
        end
        function status = checkChannelOutput(obj,iChan)
            
           %Checks if the selected channels have valid output images          
           nChanTot = numel(obj.owner_.channels_);
           if nargin < 2 || isempty(iChan), iChan = 1:nChanTot; end
           assert(all(obj.checkChanNum(iChan)));
           status =  arrayfun(@(x) exist(obj.outFilePaths_{1,x},'dir') && ...
               ~isempty(dir([obj.outFilePaths_{1,x} filesep '*.mat'])),iChan);
        end
    end
    methods (Static)
        function funParams =  getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',[owner.outputDirectory_,filesep,'ActinFilamentsAnalysisProcess'],@ischar);
            ip.parse(owner, varargin{:})
%             outputDir=ip.Results.outputDir;
            
            % Set default parameters
            chanNumStr = cellfun(@num2str,num2cell(1:length(owner.channels_)),'UniformOutput',false);
            funParams.ImageProcessingProcessIndex = owner.getProcessIndex(MinimalBridgingProcess.empty,1,false);
            funParams.OutputDirectory = ip.Results.outputDir;
            funParams.outFilePaths = strcat(ip.Results.outputDir, ...
                filesep,'Channel_',chanNumStr);
            funParams.t = 1:owner.nFrames_;
            funParams.c = 1:length(owner.channels_);
            funParams.z = 1:owner.zSize_;
            funParams.defaultOutput = 'skeleton';

        end
        function name = getName()
            name = 'Actin Filaments';
        end
    end
    
end
function actinFun(process)
    [MD,process] = process.getOwnerAndProcess('ActinFilamentsAnalysisProcess',true);
    for c = 1:length(process.outFilePaths_)
        if(~exist(process.outFilePaths_{c},'dir'))
            mkdir(process.outFilePaths_{c});
        end
    end
    params = process.getParameters();
    
    out.params = params;
       
    numImages = length(params.c)*length(params.t)*length(params.z);
    counter = 0;
    
    for c = params.c
        template = [process.outFilePaths_{c} filesep 'actin_skeleton_c%02d_t%03d_z%03d.mat'];
        for t = params.t
            for z = params.z
                progressText(counter/numImages,sprintf('Analyzing actin filaments c=%02d, t=%03d, z=%03d',c,t,z));

                bridged =  process.owner_.processes_{params.ImageProcessingProcessIndex}.loadChannelOutput(c,t,'iZ',z,'output','nms_skel_bridged_skel');
                
                if(~iscell(bridged))
                    bridged = {bridged};
                end
                
                for k = 1:length(bridged)
                    S = lamins.classes.Skeleton(bridged{k});
                    % Need to do this or need to edit
                    % Skeleton.connectedVertices to deal with short edges
                    S.convertShortEdgesToVertices(2);
                    S.reduceVerticesToPoints;
                    out.skeleton{k} = S;
                end
                save(sprintf(template,c,t,z),'-struct','out');
                counter = counter+1;
                progressText(counter/numImages,sprintf('Analyzing actin filaments c=%02d, t=%03d, z=%03d',c,t,z));
            end
        end
    end
end
