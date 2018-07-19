classdef OrientationSpacePointDetectionProcess < DetectionProcess
    %OrientationSpacePointDetectionProcess Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        imageOutput_ = [];
    end
    
    methods
        
        function obj = OrientationSpacePointDetectionProcess(owner,varargin)
%             if(nargin < 1)
%                 % Allow empty creation
%                 return;
%             end
            ip = inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('funParams', ...
                OrientationSpacePointDetectionProcess.getDefaultParams(owner), ...
                @isstruct);
            ip.parse(owner,varargin{:});
            
            obj = obj@DetectionProcess(owner, ... 
                'OrientationSpacePointDetectionProcess', ... % name
                @saveOrientationSpacePointDetection, ... % funName
                ip.Results.funParams ... % funParams
                );
            
            obj.outFilePaths_ = ip.Results.funParams.outFilePaths;
%                 owner.getChannelPaths, ... % inFilePaths_
%                 ip.Results.funParams.outFilePaths ... % outFilePaths_
%                 );
        end
        function varargout = loadChannelOutput(obj,iChan,varargin)
            % Input check
            outputList = {'points','localMaxima'};
            ip =inputParser;
            ip.StructExpand = true;
            ip.addRequired('iChan',@(x) isscalar(x) && obj.checkChanNum(x));
            ip.addOptional('iFrame',1,@(x) all(obj.checkFrameNum(x)));
            ip.addParamValue('iZ',1);
            ip.addParamValue('useCache',true,@islogical);
            ip.addParamValue('output',[], ...
                @(x) all(ismember(x,outputList)) || ...
                ~isempty(regexp(x,'^maxMeanDiff', 'once')) || ...
                ~isempty(regexp(x,'^maxMeanDiffRatio', 'once')) || ...
                ~isempty(regexp(x,'^maxRes', 'once')) ...
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
                    output = 'nlms_mip_sc1';
                end
            end
            if ischar(output),output={output}; end
            
%             switch(output{1})
%                 case 'nlms_overlay'
%                     output{1} = 'nlms';
%             end
            selectCell = cell(1,length(output));
            doNZ = false(1,length(output));
            for o=1:length(output)
                if(~isempty(regexp(output{o},'\.nz$','once')))
                    doNZ(o) = true;
                    output{o} = output{o}(1:end-3);
                end
                
                [splits,tokens] = regexp(output{o},'{([0-9]+)}$|_sc([0-9]+)$','split','tokens');
                if(~isempty(tokens))
                    output{o} = splits{1};
                    selectCell{o} = str2double(tokens{1});
                end
                
                
            end
            

            if(~isempty(regexp(output{1},'^maxMeanDiffRatioAll','once')))
                output{1} = 'maxMeanDiffRatioAll';
            elseif(~isempty(regexp(output{1},'^maxMeanDiffRatio','once')))
                output{1} = 'maxMeanDiffRatio';
            elseif(~isempty(regexp(output{1},'^maxMeanDiffAll','once')))
                output{1} = 'maxMeanDiffAll';
            elseif(~isempty(regexp(output{1},'^maxMeanDiff','once')))
                output{1} = 'maxMeanDiff';
            end

%             if(~isempty(regexp(output{1},'^nlms', 'once')))
%                 output{1} = 'nlms';
%             elseif(~isempty(regexp(output{1},'^maxima_f', 'once')))
%                 output{1} = 'maxima';
%                 output{end+1} = 'nlms';
%                 selectCell{end+1} = selectCell{1};
%                 doMip(end+1) = doMip(1);
%             elseif(~isempty(regexp(output{1},'^maxima_value', 'once')))
%                 output{1} = 'maxima_value';
%             elseif(~isempty(regexp(output{1},'^maxima', 'once')))
%                 output{1} = 'maxima';
%             elseif(~isempty(regexp(output{1},'^maxRes', 'once')))
%                 output{1} = 'maxRes';
%             end
            
            varargout = cell(1,length(output));
            
            
            % Data loading
            % load outFilePaths_{1,iChan}
            %
            params = obj.getParameters();
            template = [obj.outFilePaths_{iChan} filesep 'points_c%02d_t%03d_z%03d.mat'];
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
                case 'nlms_overlay'
                    varargout{1}=s.nlms;
                otherwise
                    varargout{1}=s.(output{1});
            end
            if(length(output) > 1 && strcmp(output{1},'maxima') && strcmp(output{2},'nlms'))
                varargout{1} = {s.maxima,s.nlms};
            end
            for o = 1:length(output)
                if(doNZ(o))
                    varargout{o} = cellfun(@(out) out(out ~= 0),varargout{1},'UniformOutput',false);
                end
                if(~isempty(selectCell{o}))
                    varargout{o} = varargout{o}{selectCell{o}};
                end
            end
        end
        function outIm = loadOutImage(obj,iChan,iFrame,varargin)
            if(isempty(obj.imageOutput_))
                output = obj.funParams_.defaultOutput;
            else
                output = obj.imageOutput_;
            end
            outIm=obj.loadChannelOutput(iChan,iFrame,'output',output,varargin{:});
        end
        function output = getDrawableOutput(obj)
            nOutput = 4;
            
%             output(1).name='NLMS Overlay';
%             output(1).var='nlms_overlay';
%             output(1).formatData=@formatNLMSForOverlay;
%             output(1).type='image';
%             output(1).defaultDisplayMethod=@ImageOverlayDisplay;
            params = obj.getParameters();
            
            output(1).name = '2D LM, Mean Orientation Response';
            output(1).var = 'localMaxima';
            output(1).formatData = @(x) mat2gray(x);
            output(1).type = 'image';
            output(1).defaultDisplayMethod = @ImageDisplay;
            

            output(2).name = '2D Local Maxima, Points';
            output(2).var = 'points';
            output(2).formatData = @(p) [p.x p.y];
            output(2).type = 'overlay';
            output(2).defaultDisplayMethod = @(iChan) LineDisplay('Marker','o','LineStyle','none');
            
            
            orientationProc = obj.owner_.getProcess(params.orientationSpaceProcessIndex);
            orientationProcParams = orientationProc.getParameters();
            
            mm = 2;
            
            for m = 1:length(orientationProcParams.responseAngularOrder)
                output(m+mm).name = sprintf('2DLM, Max-Mean K%dm%d',orientationProcParams.responseAngularOrder(m),orientationProcParams.maximaAngularOrder(m));
                output(m+mm).var = sprintf('maxMeanDiff_m%d',m);
                output(m+mm).formatData = @(x) mat2gray(x{m});
                output(m+mm).type = 'image';
                output(m+mm).defaultDisplayMethod = @ImageDisplay;
            end
            
            mm = mm + length(orientationProcParams.responseAngularOrder);
            
            for m = 1:length(orientationProcParams.responseAngularOrder)
                output(m+mm).name = sprintf('2DLM, Max-Mean DR K%dm%d',orientationProcParams.responseAngularOrder(m),orientationProcParams.maximaAngularOrder(m));
                output(m+mm).var = sprintf('maxMeanDiffRatio_m%d',m);
                output(m+mm).formatData = @(x) mat2gray(x{m});
                output(m+mm).type = 'image';
                output(m+mm).defaultDisplayMethod = @ImageDisplay;
            end
            
            mm = mm + length(orientationProcParams.responseAngularOrder);
            
            for m = 1:length(orientationProcParams.responseAngularOrder)
                output(m+mm).name = sprintf('2DLM, Max-Mean DR All K%dm%d',orientationProcParams.responseAngularOrder(m),orientationProcParams.maximaAngularOrder(m));
                output(m+mm).var = sprintf('maxMeanDiffRatioAll_a%d',m);
                output(m+mm).formatData = @(x) mat2gray(x{m});
                output(m+mm).type = 'image';
                output(m+mm).defaultDisplayMethod = @ImageDisplay;
            end
            
            mm = mm + length(orientationProcParams.responseAngularOrder);
            
            for m = 1:length(orientationProcParams.responseAngularOrder)
                output(m+mm).name = sprintf('2DLM, Max-Mean DR K%dm%d',orientationProcParams.responseAngularOrder(m),orientationProcParams.maximaAngularOrder(m));
                output(m+mm).var = sprintf('maxMeanDiffRatio_h%d',m);
                output(m+mm).formatData = @(x) x{m}(x{m} ~= 0);
                output(m+mm).type = 'graph';
                output(m+mm).defaultDisplayMethod = @HistogramDisplay;
            end

            
%                         [uMaximaAngularOrder,~,ic] = unique(params.maximaAngularOrder);
% 
%             
%             output(length(params.responseAngularOrder)*3+length(uMaximaAngularOrder)) = struct();
%             
%             for m = 1:length(params.responseAngularOrder)
%                 output(m).name=sprintf('NLMS Image K=%d, m=%d',params.responseAngularOrder(m),params.maximaAngularOrder(m));
%                 output(m).var=sprintf('nlms_m%d',m);
% %                 output(m).formatData=@(x) mat2gray(nanmax(x{m},[],3));
%                 output(m).formatData=@(x) formatNLMSForImage(x{m});
%                 output(m).type='image';
%                 output(m).defaultDisplayMethod=@(varargin) ImageDisplay('Colormap','Parula');
%             end
%             
%             mm = length(params.responseAngularOrder);
%             cm = hsv(mm);
%             
%             for m = 1:length(params.responseAngularOrder)
%                 output(m+mm).name=sprintf('NLMS Overlay K=%d, m=%d',params.responseAngularOrder(m),params.maximaAngularOrder(m));
%                 output(m+mm).var=sprintf('nlms_o%d',m);
%                 output(m+mm).formatData=@(x) formatNLMSForOverlay(x{m},cm(m,:));
%                 output(m+mm).type='overlay';
%                 output(m+mm).defaultDisplayMethod=@ImageOverlayDisplay;
%             end
%             
%             mm = length(params.responseAngularOrder)*2;
%             
%             for m = 1:length(params.responseAngularOrder)
%                 output(m+mm).name=sprintf('Max Response K=%d, m=%d',params.responseAngularOrder(m),params.maximaAngularOrder(m));
%                 output(m+mm).var=sprintf('maxRes_m%d',m);
%                 output(m+mm).formatData=@(x) mat2gray(real(x{m}));
%                 output(m+mm).type='image';
%                 output(m+mm).defaultDisplayMethod=@ImageDisplay;
%             end
%             
%             mm = length(params.responseAngularOrder)*3;
%             
%             
%             for m = 1:length(uMaximaAngularOrder)
%                 output(m+mm).name=sprintf('Max Orientation m=%d',uMaximaAngularOrder(m));
%                 output(m+mm).var=sprintf('maxima_m%d',m);
%                 output(m+mm).formatData=@(x) x{m}(:,:,1)/pi;
%                 output(m+mm).type='image';
%                 output(m+mm).defaultDisplayMethod=@(varargin) ImageDisplay('Colormap','hsv');
%             end
%             
%             mm = mm + length(uMaximaAngularOrder);
%             cm = hsv(length(params.responseAngularOrder));
%             
%             for m = 1:length(params.responseAngularOrder)
%                 output(m+mm).name=sprintf('Orientation Field K=%d, m=%d',params.responseAngularOrder(m),params.maximaAngularOrder(m));
%                 output(m+mm).var=sprintf('maxima_f%d',m);
%                 output(m+mm).formatData=@(x) maximaToOrientationField(x{1}{ic(m)},x{2}{m});
%                 output(m+mm).type='overlay';
%                 output(m+mm).defaultDisplayMethod=@(varargin) VectorFieldDisplay('Color',cm(m,:));
%             end
%             
%             mm = mm+m;
%             m = 1;
%                 output(m+mm).name=sprintf('Mean Response');
%                 output(m+mm).var=sprintf('meanResponse');
%                 output(m+mm).formatData=@(x) mat2gray(real(x));
%                 output(m+mm).type='image';
%                 output(m+mm).defaultDisplayMethod=@ImageDisplay;
%                 
%             output = flip(output);

            
%             colors = parula(numel(obj.owner_.channels_)*nOutput);
%             output(1).name='Meshwork';
%             output(1).var='S4';
%             output(1).formatData=@(x) x.getEdgeXY;
%             output(1).type='overlay';
%             output(1).defaultDisplayMethod=@(x) LineDisplay('Marker','none',...
%                 'Color',colors((x-1)*nOutput+1,:));
%             
%             output(2).name = 'Junctions';
%             output(2).var='S4v';
%             output(2).formatData=@(x) x.getVertexXY;
%             output(2).type='overlay';
%             output(2).defaultDisplayMethod=@(x) LineDisplay('Marker','x',...
%                 'LineStyle','none','Color',colors((x-1)*nOutput+2,:));
%             
%             output(3).name = 'Edge Midpoints';
%             output(3).var='S4e';
%             output(3).formatData=@(x) x.getEdgeMidPointXY;
%             output(3).type='overlay';
%             output(3).defaultDisplayMethod=@(x) LineDisplay('Marker','o',...
%                 'LineStyle','none','Color',colors((x-1)*nOutput+3,:));
%             
%             output(4).name = 'Faces';
%             output(4).var='S4f';
%             output(4).formatData=@(x) x.getFaceCentroidXY;
%             output(4).type='overlay';
%             output(4).defaultDisplayMethod=@(x) LineDisplay('Marker','^',...
%                 'LineStyle','none','Color',colors((x-1)*nOutput+4,:));
%             
%             output(1).name='Steerable Filter Response';
%             output(1).var='MAX_st_res';
%             output(1).formatData=@mat2gray;
%             output(1).type='image';
%             output(1).defaultDisplayMethod=@ImageDisplay;
%             
%             output(2).name='Non-Maximum Suppressed';
%             output(2).var='nms';
%             output(2).formatData=@mat2gray;
%             output(2).type='image';
%             output(2).defaultDisplayMethod=@ImageDisplay;
%             
%             output(3).name='Orientation Map';
%             output(3).var='orienation_map';
%             output(3).formatData=@mat2gray;
%             output(3).type='image';
%             output(3).defaultDisplayMethod=@ImageDisplay;
%             
%             output(4).name='Scale Map';
%             output(4).var='scaleMap';
%             output(4).formatData=@mat2gray;
%             output(4).type='image';
%             output(4).defaultDisplayMethod=@ImageDisplay;
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
        function out = GUI()
            out = @cliGUI;
        end
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',[owner.outputDirectory_,filesep,'OrientationSpacePointDetectionProcess'],@ischar);
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
            funParams.defaultOutput = 'points';
            funParams.orientationSpaceProcessIndex = owner.getProcessIndex(OrientationSpaceProcess.empty);
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
            name = 'OrientationSpacePointDetectionProcess';
        end
    end
    
end
function saveOrientationSpacePointDetection(process)
    [MD,process] = process.getOwnerAndProcess('OrientationSpacePointDetectionProcess',true);
    for c = 1:length(process.outFilePaths_)
        if(~exist(process.outFilePaths_{c},'dir'))
            mkdir(process.outFilePaths_{c});
        end
    end
    params = process.getParameters();
%     filter = params.filter;
    out.params = params;
    
    orientationProc = MD.getProcess(params.orientationSpaceProcessIndex);

    
%     [out.uMaximaOrder,~,out.uMaximaOrderMap] = unique(params.maximaAngularOrder);
    
%     progressText(0,sprintf('Analyzing c=%02d, t=%03d, z=%03d',params.c(1),params.t(1),params.z(1)));
    numImages = length(params.c)*length(params.t)*length(params.z)*2;
    counter = 0;
    
    for c = params.c
        template = [process.outFilePaths_{c} filesep 'points_c%02d_t%03d_z%03d.mat'];
        for t = params.t
            for z = params.z
                progressText(counter/numImages,sprintf('Analyzing Points in Orientation Space c=%02d, t=%03d, z=%03d',c,t,z));
                out.ctz = [c t z];
                
                  orientationMean = real(orientationProc.loadChannelOutput(c,t,'iZ',z,'output','meanResponse'));
                  orientationAbsoluteMaximum = orientationProc.loadChannelOutput(c,t,'iZ',z,'output','maxima_value_mip');
                  out.localMaxima = locmax2d(orientationMean,3);
                  [out.points.y,out.points.x] = find(out.localMaxima);
                  
                  nz = out.localMaxima ~= 0;
                  
                  for m = 1:length(orientationAbsoluteMaximum)
                      out.maxMeanDiff{m} = zeros(size(out.localMaxima));
                      out.maxMeanDiff{m}(nz) = orientationAbsoluteMaximum{m}(nz) - orientationMean(nz);
                      out.maxMeanDiffAll{m} = orientationAbsoluteMaximum{m} - orientationMean;
                      out.maxMeanDiffRatio{m} = zeros(size(out.localMaxima));
                      out.maxMeanDiffRatio{m}(nz) = out.maxMeanDiff{m}(nz)./out.localMaxima(nz);
                      out.maxMeanDiffRatioAll{m} = out.maxMeanDiffAll{m}./orientationMean;
                  end

%                 out.maxima = cell(1,length(out.uMaximaOrder));
%                 out.nlms = cell(1,length(params.responseAngularOrder));
%                 
%                 I = double(MD.channels_(c).loadImage(t,z));
%                 response = filter * I;
%                 out.meanResponse = mean(response.a,3);
%                 for u = 1:length(out.uMaximaOrder)
%                     tempResponse = response.getResponseAtOrderFT(out.uMaximaOrder(u));
%                     out.maxima{u} = tempResponse.getRidgeOrientationLocalMaxima;
%                     counter = counter + 1;
%                     progressText(counter/numImages,sprintf('Analyzing Orientation c=%02d, t=%03d, z=%03d',c,t,z));
%                 end
%                 lastOrder = out.uMaximaOrder(u);
%                 for m = 1:length(params.responseAngularOrder)
%                     if(params.responseAngularOrder(m) ~= lastOrder)
%                         tempResponse = response.getResponseAtOrderFT(params.responseAngularOrder(m));
%                         lastOrder = params.responseAngularOrder(m);
%                     end
%                     
%                     out.nlms{m} = tempResponse.nonLocalMaximaSuppressionPrecise(out.maxima{out.uMaximaOrderMap(m)});
%                     counter = counter + 1;
%                     progressText(counter/numImages,sprintf('Analyzing Orientation c=%02d, t=%03d, z=%03d',c,t,z));
%                     
%                     out.maxRes{m} = tempResponse.getMaxResponseFT(out.maxima{out.uMaximaOrderMap(m)});
%                     counter = counter + 1;
%                     progressText(counter/numImages,sprintf('Analyzing Orientation c=%02d, t=%03d, z=%03d',c,t,z));
%                     
%                     out.maxima_value{m} = tempResponse.interpft1(out.maxima{out.uMaximaOrderMap(m)});
                    counter = counter + 1;
                    progressText(counter/numImages,sprintf('Analyzing Points in Orientation Space c=%02d, t=%03d, z=%03d',c,t,z));
%                     
%                 end
                               
                save(sprintf(template,c,t,z),'-struct','out');
            end;
        end;
    end;
end
function imageData = formatNLMSForImage(nlms)
    nlms_mip = nanmax(real(nlms),[],3);
    imageData = mat2gray(nlms_mip);
    imageData(nlms_mip == 0) = NaN;
end
function overlayData = formatNLMSForOverlay(nlms,color)
    overlayData = zeros(size(nlms,1),size(nlms,2),4);
%     color = [1 0 0];
    nlms_mip = nanmax(nlms,[],3);
    overlayData(:,:,1:3) = bsxfun(@times,mat2gray(nlms_mip),shiftdim(color,-1));
    overlayData(:,:,4) = mat2gray(nlms_mip);
%     overlayData(:,:,4) = nlms_mip > 125;
end
function data = maximaToOrientationField(maxima,nlms)
    [X,Y] = meshgrid(1:size(maxima,1),1:size(maxima,2),1:size(maxima,3));
    nlms(nlms < 0) = 0;
    data(:,1) = X(:);
    data(:,2) = Y(:);
    data(:,3) = cos(maxima(:));
    data(:,4) = sin(maxima(:));
    data(:,3:4) = bsxfun(@times,data(:,3:4),mat2gray(nlms(:)));
    data = data(nlms(:) > 0,:);
end