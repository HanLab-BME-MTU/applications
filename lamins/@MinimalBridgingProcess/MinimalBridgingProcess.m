classdef MinimalBridgingProcess < ImageProcessingProcess & NonSingularProcess
    %MinimumBridgingProcess Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = MinimalBridgingProcess(owner,varargin)
%             if(nargin < 1)
%                 % Allow empty creation
%                 return;
%             end
            ip = inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('funParams', ...
                MinimalBridgingProcess.getDefaultParams(owner), ...
                @isstruct);
            ip.parse(owner,varargin{:});
            
            obj = obj@ImageProcessingProcess(owner, ... 
                'MinimalBridgingProcess', ... % name
                @doMinimalBriding, ... % funName
                ip.Results.funParams, ... % funParams
                owner.getChannelPaths, ... % inFilePaths_
                ip.Results.funParams.outFilePaths ... % outFilePaths_
                );
        end
        function varargout = loadChannelOutput(obj,iChan,varargin)
            % Input check
            outputList = {'nms_skel_bridged_skel','bridges'};
            ip =inputParser;
            ip.StructExpand = true;
            ip.addRequired('iChan',@(x) isscalar(x) && obj.checkChanNum(x));
            ip.addOptional('iFrame',1,@(x) all(obj.checkFrameNum(x)));
            ip.addParamValue('iZ',1);
            ip.addParamValue('useCache',true,@islogical);
            ip.addParamValue('output',[], ...
                @(x) all(ismember(x,outputList)) || ...
                ~isempty(regexp(x,'^nms_skel_bridged_skel', 'once')) || ...
                ~isempty(regexp(x,'^bridges')) ...
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
            
%             switch(output{1})
%                 case 'nlms_overlay'
%                     output{1} = 'nlms';
%             end
            selectCell = cell(1,length(output));
            doMip = false(1,length(output));
            for o=1:length(output)
                [splits,tokens] = regexp(output{o},'{([0-9]+)}$|_sc([0-9]+)$','split','tokens');
                if(~isempty(tokens))
                    output{o} = splits{1};
                    selectCell{o} = str2double(tokens{1});
                end
                
%                 if(~isempty(regexp(output{o},'_mip$','once')))
%                     doMip(o) = true;
%                 end
            end

            if(~isempty(regexp(output{1},'^nms_skel_bridged_skel', 'once')))
                output{1} = 'nms_skel_bridged_skel';
%             elseif(~isempty(regexp(output{1},'^maxima_f', 'once')))
%                 output{1} = 'maxima';
%                 output{end+1} = 'nlms';
%             elseif(~isempty(regexp(output{1},'^maxima_value', 'once')))
%                 output{1} = 'maxima_value';
%             elseif(~isempty(regexp(output{1},'^maxima', 'once')))
%                 output{1} = 'maxima';
%             elseif(~isempty(regexp(output{1},'^maxRes', 'once')))
%                 output{1} = 'maxRes';
            end
            
            
            % Data loading
            % load outFilePaths_{1,iChan}
            %
%             params = obj.getParameters();
            template = [obj.outFilePaths_{iChan} filesep 'briding_c%02d_t%03d_z%03d.mat'];
            file = sprintf(template,iChan,ip.Results.iFrame,ip.Results.iZ);
            try
                s = cached.load(file, '-useCache', ip.Results.useCache, output{:});
            catch err
%                 switch(err.identifier)
% %                     case 'MATLAB:nonExistentField'
% %                         s = load(file);
% %                         if(~isfield(s,'S4'))
% %                             % Upgrade from S3 to S4
% %                             f = ~cellfun('isempty',s.S3);
% %                             s.S4(f) = cellfun(@copy,s.S3(f),'UniformOutput',false);
% %                             cellfun(@mergeEdgesAtObsoleteVertices,s.S4(f));
% %                             save(file,'-struct','s');
% %                             s = cached.load(file,'-useCache',false,output{:});
% %                         end
%                     otherwise
%                         rethrow(err)
%                 end
                varargout{1} = {zeros(obj.owner_.imSize_)};
                return;
            end
            
%             if(isempty(ip.Results.output))
%                 if(isfield(params,'output'))
%                     output = params.output;
%                 else
%                     output = outputList{1};
%                 end
%             end
            switch(output{1})
%                 case 'vertexMovieInfo'
%                     [X,Y] = cellfun(@getVertexXY,data,'UniformOutput',false);
%                     occMap = cellfun(@getEdgeEndpointOccupancyMap,data,'UniformOutput',false);
%                     A = cellfun(@(occMap,S) occMap(vertcat(S.vertices.PixelIdxList{:})),occMap,data,'UniformOutput',false);
%                     varargout{1} = makeMovieInfo(X,Y,A);
%                 case 'edgeMovieInfo'
%                     [X,Y] = cellfun(@getEdgeMidPointXY,data,'UniformOutput',false);
%                     rp = cellfun(@(S) regionprops(S.edges,'Area'),data,'UniformOutput',false);
%                     A = cellfun(@(rp) vertcat(rp.Area),rp,'UniformOutput',false);
%                     varargout{1} = makeMovieInfo(X,Y,A);
%                 case 'faceMovieInfo'
%                     [X,Y] = cellfun(@getFaceCentroidXY,data,'UniformOutput',false);
%                     rp = cellfun(@(S) regionprops(S.faces,'Area'),data,'UniformOutput',false);
%                     A = cellfun(@(rp) vertcat(rp.Area),rp,'UniformOutput',false);
%                     varargout{1} = makeMovieInfo(X,Y,A);
%                 case 'nlms_overlay'
%                     varargout{1}=s.nlms;
                otherwise
                    varargout{1}=s.(output{1});
            end
%             if(length(output) > 1 && strcmp(output{1},'maxima') && strcmp(output{2},'nlms'))
%                 varargout{1} = {s.maxima,s.nlms};
%             elseif(doMip)
%                 varargout{1} = cellfun(@(out) nanmax(out,[],3),varargout{1},'UniformOutput',false);
%             end
            for o = 1:length(output)
%                 if(doMip(o))
%                     varargout{1} = cellfun(@(out) nanmax(out,[],3),varargout{1},'UniformOutput',false);
%                 end
                if(~isempty(selectCell{o}))
                    varargout{1} = varargout{1}{selectCell{o}};
                end
            end
        end

        function output = getDrawableOutput(obj)
%             nOutput = 4;
            
%             output(1).name='NLMS Overlay';
%             output(1).var='nlms_overlay';
%             output(1).formatData=@formatNLMSForOverlay;
%             output(1).type='image';
%             output(1).defaultDisplayMethod=@ImageOverlayDisplay;
            params = obj.getParameters();
            MD = obj.getOwner();
            orientationProcess = MD.processes_{params.orientationSpaceProcessIndex};
            orientationParams = orientationProcess.getParameters();
            
%             m = 1;
%             
% %                 output(m).name=sprintf('Width Image K=%d, m=%d',params.responseAngularOrder(m),params.maximaAngularOrder(m));
%                 output(m).name=sprintf('Width Image');
%                 output(m).var=sprintf('width');
%                 output(m).formatData=@(x) mat2gray(nanmax(x,[],3));
%                 output(m).type='image';
%                 output(m).defaultDisplayMethod=@ImageDisplay;
                
                orientationParams.maximaAngularOrder;
                
            output(length(orientationParams.maximaAngularOrder)*2) = struct();
                
            for k = 1:length(orientationParams.maximaAngularOrder)
                output(k).name=sprintf('Bridged Image K=%d, m=%d', ...
                    orientationParams.responseAngularOrder(k), ...
                    orientationParams.maximaAngularOrder(k));
                output(k).var=sprintf('nms_skel_bridged_skel_o%d',k);
                output(k).formatData=@(x) double(x{k}(:,:,1));
                output(k).type='image';
                output(k).defaultDisplayMethod= ...
                    @ImageDisplay;
            end
            
            kk = k;
            
            for k = 1:length(orientationParams.maximaAngularOrder)
                output(k+kk).name=sprintf('Bridged Image K=%d, m=%d', ...
                    orientationParams.responseAngularOrder(k), ...
                    orientationParams.maximaAngularOrder(k));
                output(k+kk).var=sprintf('nms_skel_bridged_skel_k%d',k);
                output(k+kk).formatData=@(x) formatMaskForOverlay(double(x{k}(:,:,1)),[1 0 0]);
                output(k+kk).type='overlay';
                output(k+kk).defaultDisplayMethod= ...
                    @ImageOverlayDisplay;
            end
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
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',[owner.outputDirectory_,filesep,'MinimalBridgingProcess'],@ischar);
            ip.parse(owner, varargin{:})
%             outputDir=ip.Results.outputDir;
            
            % Set default parameters
            
%             funParams.filter = OrientationSpaceFilter(1/2/pi/2,1/2/pi/2,8);
            chanNumStr = cellfun(@num2str,num2cell(1:length(owner.channels_)),'UniformOutput',false);
            funParams.OutputDirectory = ip.Results.outputDir;
            funParams.outFilePaths = strcat(ip.Results.outputDir, ...
                filesep,'Channel_',chanNumStr);
%             funParams.responseAngularOrder = [8 3 5];
%             funParams.maximaAngularOrder = [8 8 5];
            funParams.t = 1:owner.nFrames_;
            funParams.c = 1:length(owner.channels_);
            funParams.z = 1:owner.zSize_;
            funParams.defaultOutput = 'nms_skel_bridged_skel{1}';
            warning off lccb:process
            funParams.orientationSpaceProcessIndex = owner.getProcessIndex(OrientationSpaceProcess.empty,1,false);
            funParams.thresholdProcessIndex = owner.getProcessIndex(ThresholdProcess.empty,1,false);
            warning on lccb:process
%             funParams.ChannelIndex = 1:numel(owner.channels_);
%             funParams.Levelsofsteerablefilters = 2;
%             funParams.BaseSteerableFilterSigma = 1;
%             funParams.ImageFlattenFlag = 2;
%             % sub-sample number, since often VIF images are taken at a
%             % lower sample rate than the other channel, so use this number
%             % to save some time.
%             funParams.Sub_Sample_Num = 1;
        end
        function name = getName()
            name = 'MinimalBridgingProcess';
        end
        function func = GUI(varargin)
            func = @cliGUI;
        end
    end
    
end
function overlayData = formatMaskForOverlay(mask,color)
    overlayData = zeros(size(mask,1),size(mask,2),4);
%     color = [1 0 0];
    overlayData(:,:,1:3) = bsxfun(@times,mask,shiftdim(color,-1));
    overlayData(:,:,4) = mask;
%     overlayData(:,:,4) = nlms_mip > 125;
end
function doMinimalBriding(process)
    gcp;

%     I = double(I);
%     F = OrientationSpaceRidgeFilter(1./2/pi./2,[],8,'none');
%     R = I*F;
%     R3 = R.getResponseAtOrderFT(3);
%     [original.maxima,~,original.maximaV] = R.getRidgeOrientationLocalMaxima;
%     nlms = R3.nonLocalMaximaSuppressionPrecise(original.maxima);

    MD = process.getOwner();
    params = process.getParameters();
    out.params = params;
    orientationSpaceProcess = MD.processes_{params.orientationSpaceProcessIndex};
    thresholdProcess = MD.processes_{params.thresholdProcessIndex};
    
    responseAngularOrder = orientationSpaceProcess.funParams_.responseAngularOrder;

    
    numImages = length(params.c)*length(params.t)*length(params.z)*length(responseAngularOrder);
%     numImages = length(params.c)*length(params.t)*length(params.z)*(length(out.uMaximaOrder)+length(params.responseAngularOrder)*3);
    counter = 0;
    
    for c = params.c
        template = [process.outFilePaths_{c} filesep 'briding_c%02d_t%03d_z%03d.mat'];
        thresholdFile = sprintf([thresholdProcess.funParams_.OutputDirectory filesep 'threshold_values_for_channel_%d.mat'],c);
        matData = load(thresholdFile);
        thresholdValues = matData.thresholdValues;
        for t = params.t
            for z = params.z
                T = thresholdValues(t,z);
                progressText(counter/numImages,sprintf('Analyzing bridges c=%02d, t=%03d, z=%03d',c,t,z));

                nlms_by_order = orientationSpaceProcess.loadChannelOutput(c,t,'iZ',z,'output','nlms');
                
                for k = 1:length(responseAngularOrder)
    
                    nlms = nlms_by_order{k};

                    nms = nlms(:,:,1);
                    nlms_mip3 = nanmax(nlms,[],3);

                %     T = thresholdOtsu(nlms_mip3);
                    % T = 0;
                    nlms_binary = nlms_mip3 > T;
                    nms_binary = nms > T;
                    % non_nms_binary = nlms_binary & ~nms_binary;


                    nms_skel = bwmorph(nms_binary,'skel',Inf);
                    nms_skel_thresh = nms_skel.* nlms_binary;

                    nms_skel_thresh_wo_bp = lamins.functions.bwRemoveBranchPoints(nms_skel_thresh);

                    nlms_fragments = nlms_binary & ~nms_skel_thresh_wo_bp;
                    nlms_fragments_cc = bwconncomp(nlms_fragments);

                    nms_skel_wo_bp_cc = bwconncomp(nms_skel_thresh_wo_bp);

                    I = MD.channels_(c).loadImage(t,z);
                    out.bridges{k} = lamins.functions.minimalBridge(nlms_fragments_cc,nms_skel_wo_bp_cc,I);

                    nms_skel_bridged = nms_skel | out.bridges{k};
                    out.nms_skel_bridged_skel{k} = bwmorph(nms_skel_bridged,'skel',Inf);
                    
                    counter = counter+1;
                    progressText(counter/numImages,sprintf('Analyzing bridges c=%02d, t=%03d, z=%03d',c,t,z));
                end
                
                save(sprintf(template,c,t,z),'-struct','out');
            end
        end
    end

end

