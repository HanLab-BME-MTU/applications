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
    optimalHistogram(edges_cc.lengths);
end

