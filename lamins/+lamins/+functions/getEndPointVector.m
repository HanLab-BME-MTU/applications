function vector = getEndPointVector(skel,method,varargin)
    import lamins.functions.*;

        if(nargin < 2)
            method = 'local';
        end        

        %% Identify end points and branch points
        endpts = bwmorph(skel,'endpoints');
        branchpts = bwmorph(skel,'branchpoints');

        switch(method)
            case 'local'
                %% Determine angle of segment near end points
                % The dilation method would fail with nearby points
                % labeled_endpts = threshed.lm;
                % labeled_endpts(~endpts) = 0;
                % dilated_labeled_endpts = imdilate(labeled_endpts,strel('disk',5));
                % endpt_neighborhood_mask = dilated_labeled_endpts == threshed.lm;
                endpt_neighborhood_mask = bwdistgeodesic(skel,endpts) < 6;
                endpt_neighborhood = skel;
                endpt_neighborhood(~endpt_neighborhood_mask) = 0;
                cc = bwconncomp(endpt_neighborhood,8);


                rpcentroid = regionprops(cc,'Centroid');
                centroids = vertcat(rpcentroid.Centroid);
                if(~isempty(centroids))
                    centroid_idx = sub2ind(size(skel),round(centroids(:,1)),round(centroids(:,2)));
                    centroid_matrix = propmatrix(cc,centroid_idx);
                    [centroid.c,centroid.r] = ind2sub(size(skel),centroid_matrix(find(endpts)));
                    % find and centroid coordinates are reversed
                    [endpt.r,endpt.c] = find(endpts);
                    local.u = endpt.c - centroid.c;
                    local.v = endpt.r - centroid.r;
                    local.n = sqrt(local.u.^2 +local.v.^2);
                    local.x = endpt.c;
                    local.y = endpt.r;

                    vector = [local.u local.v];
                else
                    vector = [];
                end

            case 'filter'
                % need theta from steerable filter
                theta = varargin{1};

                idx = find(endpts);
                x = c;
                y = r;
                u = sin(theta(idx));
                v = -cos(theta(idx));
    
                vector = [u v];


            case 'segment'
                %% Breakup skeleton into segments and determine vectors between vertices
                branchpt.idx = find(branchpts);
                endpt.idx = find(endpts);

                vertices = branchpts | endpts;
                [edges_cc, vertices_cc, pairs] = bwtrace(skel,vertices);

                vlm = labelmatrix(vertices_cc);
                vfilt = cellfun(@length,vertices_cc.edges(vlm(endpt.idx))) == 1;
                v_endpts = vlm(endpt.idx(vfilt));

                edge_vertices = [edges_cc.vertices{[vertices_cc.edges{v_endpts}]'}];
                % reorder the vertices to have end points last
                edge_vertices = [edge_vertices(edge_vertices ~= repmat(v_endpts',2,1)) v_endpts];
                edge_vertices = reshape([vertices_cc.PixelIdxList{edge_vertices(:)}],[],2);

                [seg.row, seg.col ] = ind2sub(vertices_cc.ImageSize,edge_vertices);

                seg.u = seg.col(:,2)-seg.col(:,1);
                seg.v = seg.row(:,2)-seg.row(:,1);
                seg.n = sqrt(seg.u.^2 +seg.v.^2);

                vector =[ u v ];
        end
end % getEndPointAngle

