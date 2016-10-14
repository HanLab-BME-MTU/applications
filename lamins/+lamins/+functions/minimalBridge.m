function [ bridges, allway] = minimalBridge( fragment_cc, segment_cc )
%minimalBridge Summary of this function goes here
%   Detailed explanation goes here

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
        segment_dist{s} = bwdistgeodesic(L | Lseg, Lseg,'quasi-euclidean');
    end
%     pair_bridge = cell(1,nchoosek(length(segments),2));
%     p = 0;
    if(length(segments) > 1)
%         center = imregionalmin(nan2inf(sum(cat(3,segment_dist{:}),3)));
%         center_dist = bwdistgeodesic(L | center, center,'quasi-euclidean');
%         
%         for s = 1:length(segments)
%             D = round((segment_dist{s}+center_dist)*8)/8;
%             bridges = bridges | imregionalmin(nan2inf(D));
%         end
        for sx = 1:length(segments)-1
            for sy = sx+1:length(segments)
    %             p = p + 1;
                D = round((segment_dist{sx}+segment_dist{sy})*8)/8;
                bridges = bridges | imregionalmin(nan2inf(D));
            end
        end
    end
end


end

