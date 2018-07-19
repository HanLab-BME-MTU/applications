classdef LaminsNLMSAnalysisProcess < DetectionProcess
    %LaminsNLMSAnalysisProcess Analyzes lamin curvilinear structures
    
    properties
    end
    
    methods
        function obj = LaminsNLMSAnalysisProcess(owner)
%             obj = obj@ImageAnalysisProcess(owner,'LaminsAnalysisProcess',@lamins.functions.analyzeLaminsForProcess);
              obj = obj@DetectionProcess(owner, ...
                  'LaminsNLMSAnalysisProcess', ...
                  @lamins.functions.analyzeLaminsForProcessWithNLMS, ...
                  struct());
            
        end
        function run(obj,varargin)
            basePath = [obj.owner_.outputDirectory_ filesep 'LaminsNLMSAnalysisProcess'];
            
            changed = false;
            for c=1:length(obj.owner_.channels_)
                outPath = [basePath filesep 'channel_' num2str(c)];
                if(isempty(obj.outFilePaths_{c}))
                    obj.setOutFilePaths(outPath,c);
                    changed = true;
                end
            end
            if(changed)
                obj.owner_.save();
            end
            
            params = obj.getParameters();
            if(~isfield(params,'outFilePaths') || isempty(params.outFilePaths))
                params.outFilePaths = obj.outFilePaths_;
            end
            params.process = obj;
            
            % Call the superclass runner
            obj.run@ImageAnalysisProcess(params,varargin{:});
        end
        function varargout = loadChannelOutput(obj,iChan,varargin)
            
            % Input check
            outputList = {'S','S2','S3','S4','S4v','S4e','S4f'};
            ip =inputParser;
            ip.StructExpand = true;
            ip.addRequired('iChan',@(x) isscalar(x) && obj.checkChanNum(x));
            ip.addOptional('iFrame',1,@(x) all(obj.checkFrameNum(x)));
            ip.addParamValue('iZ',1);
            ip.addParamValue('useCache',true,@islogical);
            ip.addParamValue('output',[],@(x) all(ismember(x,outputList)));
            ip.parse(iChan,varargin{:})
            iFrame = ip.Results.iFrame;
            iZ = ip.Results.iZ;
            output = ip.Results.output;
            output = 'S4';
            if ischar(output),output={output}; end
            
            % Data loading
            % load outFilePaths_{1,iChan}
            %
            params = obj.getParameters();
            file = [obj.outFilePaths_{1,iChan} filesep 'skeletons_' params.analysisDate '.mat'];
            try
                s = cached.load(file, '-useCache', ip.Results.useCache, output{:});
            catch err
                switch(err.identifier)
                    case 'MATLAB:nonExistentField'
                        s = load(file);
                        if(~isfield(s,'S4'))
                            % Upgrade from S3 to S4
                            f = ~cellfun('isempty',s.S3);
                            s.S4(f) = cellfun(@copy,s.S3(f),'UniformOutput',false);
                            cellfun(@mergeEdgesAtObsoleteVertices,s.S4(f));
                            save(file,'-struct','s');
                            s = cached.load(file,'-useCache',false,output{:});
                        end
                    otherwise
                        rethrow(err)
                end
            end
            
            tz = sub2ind([obj.owner_.nFrames_ obj.owner_.zSize_],iFrame,iZ);
            if(isempty(ip.Results.output))
                if(isfield(params,'output'))
                    output = params.output;
                else
                    varargout{1}=s.(output{1}){tz};
                    return;
                end
            end
            data = s.S4;
            switch(output)
                case 'vertexMovieInfo'
                    [X,Y] = cellfun(@getVertexXY,data,'UniformOutput',false);
                    occMap = cellfun(@getEdgeEndpointOccupancyMap,data,'UniformOutput',false);
                    A = cellfun(@(occMap,S) occMap(vertcat(S.vertices.PixelIdxList{:})),occMap,data,'UniformOutput',false);
                    varargout{1} = makeMovieInfo(X,Y,A);
                case 'edgeMovieInfo'
                    [X,Y] = cellfun(@getEdgeMidPointXY,data,'UniformOutput',false);
                    rp = cellfun(@(S) regionprops(S.edges,'Area'),data,'UniformOutput',false);
                    A = cellfun(@(rp) vertcat(rp.Area),rp,'UniformOutput',false);
                    varargout{1} = makeMovieInfo(X,Y,A);
                case 'faceMovieInfo'
                    [X,Y] = cellfun(@getFaceCentroidXY,data,'UniformOutput',false);
                    rp = cellfun(@(S) regionprops(S.faces,'Area'),data,'UniformOutput',false);
                    A = cellfun(@(rp) vertcat(rp.Area),rp,'UniformOutput',false);
                    varargout{1} = makeMovieInfo(X,Y,A);
                otherwise
                    varargout{1}=s.(output{1}){tz};
            end
        end
        function output = getDrawableOutput(obj)
            nOutput = 4;
            colors = parula(numel(obj.owner_.channels_)*nOutput);
            output(1).name='Meshwork';
            output(1).var='S4';
            output(1).formatData=@(x) x.getEdgeXY;
            output(1).type='overlay';
            output(1).defaultDisplayMethod=@(x) LineDisplay('Marker','none',...
                'Color',colors((x-1)*nOutput+1,:));
            
            output(2).name = 'Junctions';
            output(2).var='S4v';
            output(2).formatData=@(x) x.getVertexXY;
            output(2).type='overlay';
            output(2).defaultDisplayMethod=@(x) LineDisplay('Marker','x',...
                'LineStyle','none','Color',colors((x-1)*nOutput+2,:));
            
            output(3).name = 'Edge Midpoints';
            output(3).var='S4e';
            output(3).formatData=@(x) x.getEdgeMidPointXY;
            output(3).type='overlay';
            output(3).defaultDisplayMethod=@(x) LineDisplay('Marker','o',...
                'LineStyle','none','Color',colors((x-1)*nOutput+3,:));
            
            output(4).name = 'Faces';
            output(4).var='S4f';
            output(4).formatData=@(x) x.getFaceCentroidXY;
            output(4).type='overlay';
            output(4).defaultDisplayMethod=@(x) LineDisplay('Marker','^',...
                'LineStyle','none','Color',colors((x-1)*nOutput+4,:));
        end 
%         function status = checkChannelOutput(obj,varargin)
%             status = true;
%         end
    end
    methods (Static)
        function name = getName(varargin)
            name = 'Lamins Analysis Process';
        end
        function params = getDefaultParams(varargin)
            if(nargin > 0)
                if(isa(varargin{1},'MovieData'))
                    varargin = varargin(2:end);
                end
            end
            ip = lamins.functions.analyzeLaminsForProcessParameters();
            ip.parse(varargin{:});
            params = ip.Results;
        end
    end
    
end

function movieInfo = makeMovieInfo(X,Y,A)
    movieInfo(numel(X)) = struct('xCoord',[],'yCoord',[],'zCoord',[],'amp',[]);
%                     movieInfo(numel(data)) = struct('xCoord',[],'yCoord',[],'amp',[]);
    for i = 1:numel(X)
        X{i}(end,2) = 0;
        Y{i}(end,2) = 0;
        Z{i} = zeros(size(X{i}));
        A{i}(end,2) = 0;
    end
    [movieInfo.xCoord] = X{:};
    [movieInfo.yCoord] = Y{:};
    [movieInfo.zCoord] = Z{:};
    [movieInfo.amp] = A{:};
end