classdef Skeleton < hgsetget
    properties ( Transient )
        edges
        vertices
    end
    properties
        bw
    end
    methods
        function obj = Skeleton(bw)
            if(nargin > 0)
                if(isa(bw,'lamins.classes.LaminsImage'))
                    obj.bw = bw.extendedSkel;
                else
                    obj.bw = bw;
                end
            end
        end
        function e = get.edges(obj)
            e = obj.edges;
        end
        function set.edges(obj,e)
            obj.edges = e;
            % create label matrix
            obj.edges.lm = labelmatrix(e);
            edgesbw = obj.edges.lm > 0;
            % sort PixelIdxList along the edges
            endpts = find(bwmorph(edgesbw,'endpoints'));
            edges_startpt = zeros(1,obj.edges.NumObjects);
            edges_startpt(obj.edges.lm(endpts)) = endpts;
            geo = bwdistgeodesic(edgesbw,edges_startpt(edges_startpt ~= 0));
            I = obj.edges.PixelIdxList;
            I = cellfun(@(x) sortrows([geo(x) x]),I,'UniformOutput',false);
            I = cellfun(@(x) x(:,2),I,'UniformOutput',false);
            obj.edges.PixelIdxList = I;
        end
        function v = get.vertices(obj)
            v = obj.vertices;
        end
        function set.vertices(obj,v)
            obj.vertices = v;
            obj.vertices.lm = labelmatrix(v);
        end
        function set.bw(obj,bw)
            obj.bw = bw;
            % check skeletonization
            obj.bw = bwmorph(obj.bw,'skel',Inf);
            n8 = obj.countNeighbors;
            % edges only have 2 8-connected neighbors
            obj.edges = bwconncomp(n8 == 2);
            % vertices have more than 2 8-connected neighbors
            obj.vertices = bwconncomp(n8 > 2);
            disp('set.bw');
        end
        function v = connectedVertices(obj,e)
            if(nargin < 2)
                e = 1:obj.edges.NumObjects;
            end
            v = zeros(length(e),2);
            dilated = imdilate(obj.vertices.lm,strel('square',3));
            for i=1:length(e)
                v(i,:) = dilated(obj.edges.PixelIdxList{e(i)}([1 end]));
            end
        end
        function e = connectedEdges(obj,v)
            if(nargin < 2)
                v = 1:obj.vertices.NumObjects;
            end
            %e = cell(1,length(v));
            dilated = imdilate(obj.vertices.lm,strel('square',3));
            e = regionprops(dilated,obj.edges.lm,'PixelValues');
            e = { e(v).PixelValues };
            e = cellfun(@(x) unique(x(x ~= 0)),e,'UniformOutput',false);
        end
        function S = deleteEdges(obj,e)
            bw = obj.bw;
            for i = 1:length(e)
                bw(obj.edges.PixelIdxList{e(i)}) = 0;
            end
            S = Skeleton(bw);
        end
        function E = edgeIntensities(obj,I)
            rp = regionprops(obj.edges,I,'PixelValues');
            E = { rp.PixelValues };
        end
        function E = edgeIntensitiesMatrix(obj,I,width)
            Ecell = edgeIntensities(obj,I);
            if(nargin < 3)
                width = max(cellfun(@length,Ecell));
            end
            E = zeros(length(Ecell),width);
            for i=1:length(Ecell)
                if(nargin < 3)
                    E(i,1:length(Ecell{i})) = Ecell{i};
                elseif(length(Ecell{i}) > 1)
                    E(i,:) = interp1(double(Ecell{i}),linspace(1,length(Ecell{i}),width));
                else
                    E(i,:) = Ecell{i};
                end
            end
        end
        function neighbors = countNeighbors(obj,conn)
            if(nargin < 2)
                conn = 8;
            end
            if(conn == 4)
                filter = [ 0 1 0 ; 1 0 1; 0 1 0];
            elseif(conn == 8)
                filter = [ 1 1 1 ; 1 0 1; 1 1 1];
            end
            neighbors = imfilter(double(obj.bw),filter).*obj.bw;
        end
        function convertShortEdgesToVertices(obj,maxLength)
            import lamins.functions.*;
            lengths = cellfun(@length,obj.edges.PixelIdxList);
            shortEdges = filtercc(obj.edges,lengths <= maxLength);
            shortEdges.lm = labelmatrix(shortEdges);
            obj.vertices = bwconncomp(shortEdges.lm > 0 | obj.vertices.lm > 0);
            obj.edges = filtercc(obj.edges,lengths > maxLength);
        end
        function rgb = showGraph(obj,I)
            rgb = zeros([size(obj.bw) 3]);
            rgb(:,:,1) = obj.edges.lm > 0;
            rgb(:,:,2) = obj.vertices.lm > 0;
            if(nargin > 1)
                I = double(I);
                rgbI = rgb;
                rgbI(:,:,3) = I .* ~any(rgb,3);
                rgbI(:,:,1) = I .* ~rgb(:,:,2);
                rgbI(:,:,2) = I .* ~rgb(:,:,1);
                rgb = rgbI;
            end
            imshow(rgb);
        end
    end
end
