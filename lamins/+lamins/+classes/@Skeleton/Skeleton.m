classdef Skeleton < hgsetget &  matlab.mixin.Copyable
    % Class to contain vertices, edges, and faces of a skeleton meshwork
    properties ( Transient )
        edges
        vertices
        assumeOrdered = false;
        faces
    end
    properties
        bw
    end
    methods
        function obj = Skeleton(bw)
            if(nargin > 0)
                if(isa(bw,'lamins.classes.LaminsImage'))
                    obj.bw = bw.auditedSkel;
%                     obj.bw = bw.extendedSkel;
                else
                    obj.bw = bw;
                end
            end
        end
        function e = get.edges(obj)
            if(isempty(obj.edges))
                n8 = obj.countNeighbors;
                % edges only have 2 8-connected neighbors
                obj.edges = bwconncomp(n8 == 2);
            end
            e = obj.edges;
        end
        function set.edges(obj,e)
            obj.edges = e;
            % create label matrix
            obj.edges.lm = labelmatrix(e);
            if(~obj.assumeOrdered)
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
                obj.edges.SortedPixels = true;
                obj.assumeOrdered = true;
            end
        end
        function deleteEdges(obj,e)
            import connectedComponents.*;
            bw = obj.bw;
            for i = 1:length(e)
                bw(obj.edges.PixelIdxList{e(i)}) = 0;
            end
            obj.bw = bw;
            
            filter = true(1,length(obj.edges.PixelIdxList));
            filter(e) = false;
            obj.edges = ccFilter(obj.edges,filter);
        end
        function deleteEdgeLoops(obj)
            % an edge loop is an edge that starts and ends at the same
            % place
            edgeLoops = cellfun(@(x) x(1) == x(end),S.edges.PixelIdxList);
            obj.deleteEdges(edgeLoops);
        end
        function deleteFaces(obj,f)
            import connectedComponents.*;
            bw = obj.bw;
            for i = 1:length(f)
                bw(obj.faces.PixelIdxList{f(i)}) = 1;
            end
            obj.bw = bw;
            
            filter = true(1,length(obj.faces.PixelIdxList));
            filter(f) = false;
            obj.faces = ccFilter(obj.faces,filter);
        end
        function addEdgeBetweenPoints(obj,p1,p2)
            import connectedComponents.*;
            if(isscalar(p1))
                [p1(1),p1(2)] = ind2sub(obj.edges.ImageSize,p1);
            end
            if(isscalar(p2))
                [p2(1),p2(2)] = ind2sub(obj.edges.ImageSize,p2);
            end
            P = bresenham(p1,p2,8);
            obj.edges = ccAppend(obj.edges,single(P));
        end
        function simplifyLoops(obj)
            import connectedComponents.*;
            import lamins.functions.*;
            faceEdges = obj.faceEdges;
            faceLoops = cellfun(@length,faceEdges) == 2;
            faceEdges = faceEdges(faceLoops);
            for fl = 1 : length(faceEdges)
                endpts = getEdgeEndpoints(ccFilter(obj.edges,faceEdges{fl}));
                endpts = unique(endpts);
                obj.addEdgeBetweenPoints(endpts(1),endpts(2));
            end
            obj.deleteEdges(unique([faceEdges{:}]));
            obj.bw = [];
            obj.faces = [];
        end
        function v = get.vertices(obj)
            if(isempty(obj.vertices))
                n8 = obj.countNeighbors;
                % vertices have more than 2 8-connected neighbors
                obj.vertices = bwconncomp(n8 > 2);
            end
            v = obj.vertices;
        end
        function set.vertices(obj,v)
            obj.vertices = v;
            obj.vertices.lm = labelmatrix(v);
        end
        function set.bw(obj,bw)
            % check skeletonization
            if(isempty(bw))
                obj.bw = bw;
            else
                obj.bw = bwmorph(bw,'skel',Inf);
            end
        end
        function bw = get.bw(obj)
            if(isempty(obj.bw))
                % if empty, generate bw from edge
                obj.bw = zeros(1024);
                obj.bw(vertcat(obj.edges.PixelIdxList{:})) = 1;
            end
            bw = obj.bw;
        end
        function faces =get.faces(obj)
            if(isempty(obj.faces))
                faces = bwconncomp(~bwmorph(obj.edges.lm > 0 | obj.vertices.lm > 0,'diag'));
                faces.lm = labelmatrix(faces);
                filter = true(1,faces.NumObjects);
                filter( nonzeros(faces.lm(~imfill(obj.bw,'holes'))) ) = false;
                faces = connectedComponents.ccFilter(faces,filter);
                faces.lm = labelmatrix(faces);
                obj.faces = faces;
            end
            faces = obj.faces;
        end
        function [edgeIndices,faceIndices] = faceEdges(obj,f)
            import connectedComponents.*;
            import lamins.functions.*;
            
            dilatedFaces = obj.faces;
            if(nargin > 1)
                dilatedFaces = ccFilter(dilatedFaces,f);
            end
            dilatedFaces = ccDilate(dilatedFaces,strel('square',5));
            
            
            lm = obj.edges.lm;
            lm([obj.vertices.PixelIdxList{:}]) = 0;
            % identify edges indices which overlap the dilated face
            edgeIndices = cellfun(@(x) unique(nonzeros(lm(x))),dilatedFaces.PixelIdxList,'UniformOutput',false);
            map = zeros(1,prod(obj.edges.ImageSize));
            edges = obj.edges;
            
            % obtain the linear index of the endpoints of the overlapping edges
            for faceIdx = 1:dilatedFaces.NumObjects
                edgeCandidates = ccFilter(edges,edgeIndices{faceIdx});
                
                endpts = getEdgeEndpoints(edgeCandidates);
                
                % each endpoint should be repeated exactly twice
                [mult , uniq] = getMultiplicityInt(endpts);
%                 map = sparse(1,double(uniq),mult);
%                 map = zeros(1,max(uniq));
                map(uniq) = mult;
%                 map = accumarray(endpts(:),1)';
                % remove hanging endpoints that do not intersect the other
                % edges
                goodEdges = ~any(map(endpts) == 1,2);
                
                % only consider remaining good edges
                edgeIndices{faceIdx} = edgeIndices{faceIdx}(goodEdges,:);
                edgeCandidates = ccFilter(edgeCandidates,goodEdges);
                
                midpoints = getEdgeMidpoints(edgeCandidates);
                goodEdges = ismember( midpoints,  dilatedFaces.PixelIdxList{faceIdx} );
                
                edgeIndices{faceIdx} = edgeIndices{faceIdx}(goodEdges,:);
%                 endpts = endpts(goodEdges,:);
            end
            if(nargout > 1)
                % calculate a map from edges to faces
                edgeIdx = vertcat( edgeIndices{:} );
                steps =  cumsum(cellfun(@length,edgeIndices))+1;
                % get multiplicity in case a face has no edges
                [r,u] = getMultiplicityInt([1 steps(1:end-1)]);
                runs = zeros(size(edgeIdx));
                runs(u) = r;
                faceIdx = cumsum(runs);
                % in case the last face has no edges
                faceIdx = faceIdx(1:length(edgeIdx));
                faceIndices = accumarray(edgeIdx,faceIdx,[],@(x) {x});
            end
        end
        function A = getEdgeAdjacency(obj)
            % get edge adjacency matrix where adjacency occurs when edges
            % share a face
            FE = obj.faceEdges;
            A = false(obj.edges.NumObjects);
            for i=1:length(FE)
                A(FE{i},FE{i}) = 1;
            end
            A = A & ~eye(size(A));
        end
        function A = getFaceAdjacency(obj)
            % get face adjacency matrix where adjacency occurs when faces
            % share an edge
            [FE,EF] = obj.faceEdges;
            A = false(obj.faces.NumObjects);
            for i=1:length(EF)
                A(EF{i},EF{i}) = 1;
            end
            A = A & ~eye(size(A));
        end
        function rp = getEdgeProperties(obj,I,e)
            import connectedComponents.*;
            if(nargin < 3)
                e = 1:obj.edges.NumObjects;
                edges = obj.edges;
            else
                edges = ccFilter(obj.edges, e);
            end
            rp = regionprops(edges,I,'MinIntensity','MaxIntensity','MeanIntensity','Orientation','Area','WeightedCentroid');
            mid = cellfun(@(x) I( x(ceil(end/2)) ), edges.PixelIdxList, 'UniformOutput', false);
            [rp.MiddleIntensity] = mid{:};
            firstQuartile = cellfun(@(x) prctile( I(x) ,25), edges.PixelIdxList, 'UniformOutput' , false);
            medianIntensity = cellfun(@(x) prctile( I(x) ,50), edges.PixelIdxList, 'UniformOutput' , false);
            thirdQuartile = cellfun(@(x) prctile( I(x) ,75), edges.PixelIdxList, 'UniformOutput' , false);
            [rp.FirstQuartileIntensity]  = firstQuartile{:};
            [rp.MedianIntensity] = medianIntensity{:};
            [rp.ThirdQuartileIntensity] = thirdQuartile{:};
        end
        function rp = getFaceProperties(obj,I,f)
            import connectedComponents.*;
            if(nargin < 3)
                f = 1:obj.faces.NumObjects;
                faces = obj.faces;
            else
                faces = ccFilter(obj.faces, f);
            end
            I = double(I);
            rp = regionprops(faces,I,'MinIntensity','MaxIntensity','MeanIntensity','Area','Centroid','WeightedCentroid');
            D = bwdist(obj.bw);
            distances = cellfun(@(x) D(x) , faces.PixelIdxList, 'UniformOutput', false);
            [rp.Distances] = distances{:};
            
            [uniqueDistances, ~ , UtoIn] = cellfun(@(x) unique(x), {rp.Distances}, 'UniformOutput', false);
            [rp.UniqueDistances] = uniqueDistances{:};
            
            nDistance = cellfun(@(uidx) accumarray(uidx,1), UtoIn, 'UniformOutput', false);
            [rp.N_vs_Distance] = nDistance{:};
            
            meanIntensityDistance = cellfun(@(idx,uidx) accumarray(uidx,I(idx),[],@mean),faces.PixelIdxList,UtoIn , 'UniformOutput', false);
            [rp.MeanIntensity_vs_Distance] = meanIntensityDistance{:};
            
            stdIntensityDistance = cellfun(@(idx,uidx) accumarray(uidx,I(idx),[],@std),faces.PixelIdxList,UtoIn , 'UniformOutput', false);
            [rp.StdIntensity_vs_Distance] = stdIntensityDistance{:};
            
            distanceWeightedIntensity = cellfun(@(x) sum(D(x).*I(x))./sum(D(x)), faces.PixelIdxList , 'UniformOutput', false);
            [rp.DistanceWeightedIntensity] = distanceWeightedIntensity{:};
        end
        function plotFaceProperties(obj,I,fidx,varargin)
            % f should be scalar
            rp = getFaceProperties(obj,I,fidx);
            for f=1:length(fidx)
                plot(rp(f).UniqueDistances,rp(f).MeanIntensity_vs_Distance,varargin{:});
                errorbar(rp(f).UniqueDistances,rp(f).MeanIntensity_vs_Distance,rp(f).StdIntensity_vs_Distance./rp(f).N_vs_Distance);
                hold on;
            end
            hold off;
        end
        function v = connectedVertices(obj,e)
            if(nargin < 2)
                e = 1:obj.edges.NumObjects;
            end
            v = zeros(length(e),2);
            se = strel('square',3);
%             import connectedComponents.*;
%             cc{1} = obj.edges(e);
%             cc{2} = obj.edges(e);
%             cc{1}.PixelIdxList = cellfun(@(x) x(1)  , obj.edges.PixelIdxList ,'UniformOutput',false);
%             cc{2}.PixelIdxList = cellfun(@(x) x(end), obj.edges.PixelIdxList ,'UniformOutput',false);
%             cc{1} = ccDilate(cc{1},se);
%             cc{2} = ccDilate(cc{2},se);
            dilated = imdilate(obj.vertices.lm,se);
            E = obj.edges;
            for i=1:length(e)
                v(i,:) = dilated(E.PixelIdxList{e(i)}([1 end]));
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
        function varargout = connectedEndPoints(obj,v)
            if(nargin < 2)
                v = 1:obj.vertices.NumObjects;
            end
            %e = cell(1,length(v));
            import connectedComponents.*;
            se = strel('square',3);
            dilated = ccDilate(obj.vertices,se);
            e = regionprops(dilated,obj.edges.lm,'PixelValues','PixelIdxList');
%             dilated = imdilate(obj.vertices.lm,strel('square',3));
%             e = regionprops(dilated,obj.edges.lm,'PixelValues');
            
            [varargout{1:nargout}] = arrayfun(@(s) ind2sub(obj.vertices.ImageSize,s.PixelIdxList(s.PixelValues > 0)),e(v),'UniformOutput',false);
        end
%         function S = deleteEdges(obj,e)
%             bw = obj.bw;
%             for i = 1:length(e)
%                 bw(obj.edges.PixelIdxList{e(i)}) = 0;
%             end
%             S = Skeleton(bw);
%         end
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
        function [E,V] = reduceVerticesToPoints(obj)
            rp = regionprops(obj.vertices,'Centroid');
            centroids = vertcat(rp.Centroid);
            centroids = round(centroids);
            centroidIdx = sub2ind(obj.vertices.ImageSize,centroids(:,2),centroids(:,1));
            v = obj.connectedVertices;
            E = obj.edges;
            V = obj.vertices;
            for i=1:E.NumObjects
                if(all(v(i,:) ~= 0))
                    [r,c] = ind2sub([1024 1024],obj.edges.PixelIdxList{i}([1 end]));
                    pre  = bresenham( centroids(v(i,1), [2 1]) , [r(1) c(1)]);
                    post = bresenham([r(2) c(2)], centroids(v(i,2),[2 1]) );
                    pre = sub2ind(E.ImageSize,pre(:,1),pre(:,2));
                    post = sub2ind(E.ImageSize,post(:,1),post(:,2));
                    E.PixelIdxList{i} = [ pre(1:end-1); obj.edges.PixelIdxList{i}; post(2:end) ];
                end
            end
            V.PixelIdxList = mat2cell(centroidIdx(:),ones(obj.vertices.NumObjects,1));
            obj.assumeOrdered = true;
            obj.edges = E;
            obj.vertices = V;
        end
        function drawLines(obj,e,edgeColor)
            if(nargin < 2 || isempty(e))
                e = 1:obj.edges.NumObjects;
            end
            if(islogical(e))
                e = find(e);
            end
            if(nargin < 3 || isempty(edgeColor))
                edgeColor = 'r';
            end
            rp = regionprops(obj.vertices,'Centroid');
            v = obj.connectedVertices;
            for i=e(:)'
                [r,c] = ind2sub([1024 1024],obj.edges.PixelIdxList{i});
                if(all(v(i,:) ~= 0))
                    r = [rp(v(i,1)).Centroid(2) ; r ;  rp(v(i,2)).Centroid(2)];
                    c = [rp(v(i,1)).Centroid(1) ; c ;  rp(v(i,2)).Centroid(1)];
                end
                line(c,r,'Color',edgeColor);
            end
%             [R,C] = obj.connectedEndPoints;
%             for i=1:obj.vertices.NumObjects
%                 for j = 1:length(R{i})
%                     line([rp(i).Centroid(1) C{i}(j)], [rp(i).Centroid(2) R{i}(j)],'Color',vertexColor);
%                 end
%             end
        end
        function drawEdgesAsLines(obj,e,edgeColor)
            % Draw just the edges (not the vertices as above)
            if(nargin < 2 || isempty(e))
                e = 1:obj.edges.NumObjects;
            end
            if(nargin < 3 || isempty(edgeColor))
                edgeColor = 'r';
            end
            for i=e(:)'
                [r,c] = ind2sub([1024 1024],obj.edges.PixelIdxList{i});
                line(c,r,'Color',edgeColor);
            end
        end
        function imshow(obj)
            showGraph(obj);
        end
        function [e,f] = cleanup(obj)
            obj.convertShortEdgesToVertices(2);
            obj.reduceVerticesToPoints;
            % Cleans up the edges and faces of the skeleton
            % 1. Removes edges that have no faces
            % 2. Removes faces that have no edges
            % Future considerations: Remove small faces
            [~,EF] = obj.faceEdges;
            e = find(cellfun(@length,EF) == 0);
            e = [ e ; (length(EF)+1:obj.edges.NumObjects)'];
            obj.deleteEdges(e);
            FE = obj.faceEdges;
            f = find(cellfun(@length,FE) == 0);
            obj.deleteFaces(f);
        end
        function score = getEdgeScore(obj,e)
        end
        function score = getFaceScore(obj,f)
        end
        function score = getVertexScore(obj,v)
        end
    end
end
