classdef LaminsData
    %LaminsData Represents a movieData containing information about Lamins
    %A, B1, and B2
    %   Detailed explanation goes here
    
    properties
        movieData;
        cellReader;
        params;
        outputDirectory_;
        thumbs;
    end
    
    methods
        function obj= LaminsData(MD)
            if(nargin > 0)
                obj.movieData = MD;
                obj.cellReader = CellReader(CachedReader(MD));
                obj.params = laminParams(MD);
                obj.outputDirectory_ =  [MD.outputDirectory_ filesep 'lamins'];
            end
        end
        function corrmatrix = calcCorrMatrix(obj,cache)
            if(nargin < 2)
                cache = true;
            end
            filename = [ obj.outputDirectory_ filesep 'corrmatrix.mat'];
            if(cache && exist(filename,'file'))
                C = load(filename);
                corrmatrix = C.corrmatrix;
            else
                % not cached
                corrmatrix = calcCorrMatrix(obj.movieData);
                save(filename,'corrmatrix');
            end
        end
        function plotCorrMatrix(obj)
            order = obj.params.channels.order;
            corrmatrix = obj.calcCorrMatrix();
            corrmatrix = corrmatrix(order,:,order,:);
            plotcorrmatrix(corrmatrix,obj.movieData.getFilename());
        end
        function thumbs = getThumbs(obj,scale,subs)
            if(nargin < 2 && ~isempty(obj.thumbs))
                thumbs = obj.thumbs;
                return;
            end
            if(nargin < 2 || isempty(scale))
                scale = 0.25;
            end
            I = obj.cellReader;
            if(nargin < 3)
                % show all planes
                I = I(:,:,obj.params.goodZ);
            else
                I = I(subs{:});
            end
            order = obj.params.channels.order;
            thumbs = vertcat( ...
                imadjust(horzcat(I{order(1),:})), ...
                imadjust(horzcat(I{order(2),:})), ...
                imadjust(horzcat(I{order(3),:})), ...
                imadjust(horzcat(I{order(4),:})) ...
            );
            thumbs = imresize(thumbs,scale);
        end
        function showThumbs(obj)
            thumbs = obj.getThumbs();
            imtool(thumbs,'colormap',jet(256));
        end
        function S = getSkeleton(obj,varargin)
            I = obj.cellReader{varargin{:}};
            S = getLaminSkeleton(I);
        end
        function [edges_cc,vertices_cc, pairs, skel] = getSkeletonGraph(obj,varargin)
            skel = obj.getSkeleton(varargin{:});
            branch_pts = bwmorph(skel,'branchpoints');
            end_pts = bwmorph(skel,'endpoints');
            vertices = branch_pts | end_pts;
            [edges_cc, vertices_cc, pairs] = bwtrace(skel,vertices);
        end
        function showSkeleton(obj,varargin)
            [edges_cc,vertices_cc, pairs, skel] = obj.getSkeletonGraph(varargin{:});
            figure;
            imshow(labelmatrix(edges_cc),[]);
            hold on;
            cm = jet(edges_cc.NumObjects);
            cm = cm(randperm(edges_cc.NumObjects),:);
            cm = [ [0 0 0]; cm];
            colormap(cm);
            for i=1:edges_cc.NumObjects
                [y,x] = ind2sub(edges_cc.ImageSize,[vertices_cc.PixelIdxList{edges_cc.vertices{i}}]);
                line(x,y,'Color',cm(i+1,:));
            end
            [r,c] = ind2sub(vertices_cc.ImageSize,[vertices_cc.PixelIdxList{:}]);
            scatter(c,r);
            figure;
            histogram(edges_cc.lengths);
        end
        function [mask,thresh] = getNuclearMask(obj,varargin)
            I = horzcat(obj.cellReader{varargin{:}});
            [mask, thresh] = generateMask(I);
        end
        % plot the cummulative distribution functions of the whole image,
        % masked, and non-masked areas for comparison
        % if thresh is given, show the cutoffs
        function plotMaskedCdf(obj,thresh,varargin)
            [mask,maskThresh] = obj.getNuclearMask(varargin{:});
            %I = cat(3,obj.cellReader{varargin{:}});
            I = obj.cellReader{varargin{:}};
            I;
        end
    end
    
end

