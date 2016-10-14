function [ bridges, allway] = minimalBridge( fragment_cc, segment_cc, I )
%minimalBridge Summary of this function goes here
%   Detailed explanation goes here

I(I < 0) = 0;
I = imcomplement(mat2gray(double(I)));


import lamins.functions.*;

fragment_dilated_cc = connectedComponents.ccDilate(fragment_cc,ones(3));
segment_label = labelmatrix(segment_cc);
fragment_label = labelmatrix(fragment_cc);
bridges = zeros(fragment_cc.ImageSize);

for f = 1:fragment_cc.NumObjects
    L = fragment_label == f;
    segments = unique(segment_label(fragment_dilated_cc.PixelIdxList{f}));
    segments = segments(segments ~= 0);
    segment_dist = cell(1,length(segments));
    for s = 1:length(segments)
        Lseg = segment_label == segments(s);
%         segment_dist{s} = bwdistgeodesic(L | Lseg, Lseg,'quasi-euclidean');
        A = I;
        A(~L & ~Lseg) = Inf;
        segment_dist{s} = graydist(A, Lseg,'chessboard');
    end
    if(length(segments) > 1)
        np = nchoosek(length(segments),2);
        pair_bridge = cell(1,np);
        mindist = zeros(1,np);
        p = 0;
%         center = imregionalmin(nan2inf(sum(cat(3,segment_dist{:}),3)));
%         center_dist = bwdistgeodesic(L | center, center,'quasi-euclidean');
%         
%         for s = 1:length(segments)
%             D = round((segment_dist{s}+center_dist)*8)/8;
%             bridges = bridges | imregionalmin(nan2inf(D));
%         end
        for sx = 1:length(segments)-1
            for sy = sx+1:length(segments)
                p = p + 1;
                D = segment_dist{sx}+segment_dist{sy};
                D = round(D*100)/100;
%                 D = round((segment_dist{sx}+segment_dist{sy})*8)/8;
                mindist(p) = nanmin(D(:));
                pair_bridge{p} = D == mindist(p);
%                 pair_bridge{p} = imregionalmin(nan2inf(D));
%                 mindist(p) = unique(D(pair_bridge{p}));
            end
        end
        % Find the minimum spanning tree
        % That is the tree with the minimum distance weight that connects
        % all the segmments
        min_span_tree = graphminspantree(sparse(squareform(mindist)));
        % Map square indices to pair number
        ind2p = squareform(1:np);
        pair_bridge = pair_bridge(ind2p(find(min_span_tree)));
        bridges = bridges | any(cat(3,pair_bridge{:},zeros(201)),3);
    end
end


end

