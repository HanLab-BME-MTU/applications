function [ bridges, allway] = minimalBridge( fragment_cc, segment_cc, I )
%minimalBridge Find parsimonious connections over fragments that connect
%segments as guided by image intensity
%
% INPUT
% fragment_cc - connected components structure representing fragments
% segment_cc  - connected components structure representing segments
% I - image or response map
%
% OUTPUT
% bridges - a minimum subset of fragment_cc pixels that connect segment_cc such
% that segment_cc is well connected

%% Image map processing
% Suppress negative pixels
I(I < 0) = 0;
% Take the complement such that high intensity pixels have low cost
I = imcomplement(mat2gray(double(I)));


import lamins.functions.*;

% Dilate fragments individually to detect which segments are adjacent to
% each fragment
fragment_dilated_cc = connectedComponents.ccDilate(fragment_cc,ones(3));
segment_label = labelmatrix(segment_cc);
fragment_label = labelmatrix(fragment_cc);
% Bridges will be a binary map consisting of a subset of fragment_cc pixels
bridges = false(fragment_cc.ImageSize);

% Loop over each fragment
parfor f = 1:fragment_cc.NumObjects
    % L is true when the pixels are in fragment f
    L = fragment_label == f;
    % Obtain the list of segments adjacent to fragment f
    segments = unique(segment_label(fragment_dilated_cc.PixelIdxList{f}));
    % Don't include the background
    segments = segments(segments ~= 0);
    % Determine the geodesic distance along the fragment from each segment
    segment_dist = cell(1,length(segments));
    for s = 1:length(segments)
        % Lseg is true when the segment 
        Lseg = segment_label == segments(s);
%         segment_dist{s} = bwdistgeodesic(L | Lseg, Lseg,'quasi-euclidean');
        A = I;
        % The cost of traveling outside the fragment or segment is infinite
        % See bwdistgeodesic source code
        A(~L & ~Lseg) = Inf;
        % chessboard: The distance for diagonals is the same as up or down
        segment_dist{s} = graydist(A, Lseg,'chessboard');
    end
    % If fragment is adjacent to more than one segment, then ...
    if(length(segments) > 1)
        % Measure the distance along the fragment between each pair of
        % segments
        % Determine number of segment pairs
        np = nchoosek(length(segments),2);
        % Store the connection between each segment
        pair_bridge = cell(1,np);
        % For each pair, determine the distance of the smallest path
        mindist = zeros(1,np);
        % Counter to track pair number
        p = 0;
        % Loop over each pair of segments
        for sx = 1:length(segments)-1
            for sy = sx+1:length(segments)
                p = p + 1;
                % Minimum path distance will be revealed by the sum of the
                % distances from each segment in the pair
                D = segment_dist{sx}+segment_dist{sy};
                % Round to the nearest 100th
                D = round(D*100)/100;
                % Distance between each segment is the minimum of the
                % summed distance
                mindist(p) = nanmin(D(:));
                % The bridge is the path through the fragment where the
                % summed distance is equal the minimum
                pair_bridge{p} = D == mindist(p);
            end
        end
        % Find the minimum spanning tree
        % That is the tree with the minimum distance weight that connects
        % all the segmments
        min_span_tree = graphminspantree(sparse(squareform(mindist)));
        % Map square indices to pair number
        ind2p = squareform(1:np);
        % Select the bridges that provide the minimal connectivity
        pair_bridge = pair_bridge(ind2p(min_span_tree > 0));
%         pair_bridge = pair_bridge(ind2p(find(min_span_tree)));
        
        bridges = bridges | any(cat(3,pair_bridge{:},zeros(size(I))),3);
    end
end


end

