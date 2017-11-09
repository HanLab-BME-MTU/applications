function BL = cutInteriorPolygon(CC, distMax, tol)
% This function cuts from corner to corner, looking for lumps
% characterized by a large difference in euclidean versus geodesic distance 
% between corners.
%
%   CC = connected component structure with corners as the odd CC's
%   tRad = radius of lumps 
%   tol = tolerance for difference in distance
%   BL = a label matrix with 1 is non-corner, 

    dim = CC.ImageSize;
    num = CC.NumObjects;
    nC = ceil(num/2);% number of corners
    geod = zeros(nC,1); % geodisic distance between each pair of
    stC = zeros(nC,2); % arrays of start and end points for each corner
    epC = zeros(nC,2);
   
    for j = 1:nC
        k = 2*(j-1) + 1;
      
        [y1, z1] = ind2sub(dim, CC.PixelIdxList{k}(1)); % start point
        stC(j,:) = [y1, z1];
        [y1, z1] = ind2sub(dim, CC.PixelIdxList{k}(end)); % end point
        epC(j,:) = [y1, z1];
        
        % geod(j) = length(CC.PixelIdxList{k+1}); this is chessboard
        % distance, need the eulerian
        tot = 0;
        temp = [CC.PixelIdxList{k}(end);
                CC.PixelIdxList{k+1};
                CC.PixelIdxList{mod(k+1, num)+1}(1)]; % list of indices along a non corner
        nPix = length(temp);
        for k = 2:nPix
            [y1, z1] = ind2sub(dim, temp(k-1));
            [y2, z2] = ind2sub(dim, temp(k));
            tot = tot + sqrt((y2-y1)^2+(z2-z1)^2);
        end
        geod(j) = tot;
    end
    
    % This may need to be changed when the output is chosen
    BL = zeros(dim);
    for j = 1:num
        BL(CC.PixelIdxList{j}) = mod(j,2)+1;
    end
    
    NN = 12; % number of neighbors *in front* to consider
    temp = zeros(NN,2);
    for j = 1:nC
        for k = 1:NN
            NNIdx = mod((j:(j+k-1))-1,nC) + 1;
            temp(k,1) = norm(epC(j,:)-stC(mod(NNIdx(end),nC)+1,:));
            temp(k,2) = sum(geod(NNIdx));
        end
        temp = temp(temp(:,2)<distMax,:); % truncate past a certain euclidian distance
        
        if ~isempty(temp)
            % threshold the ratio of geodesic to euclidian distance
            temp = temp(:,2)./temp(:,1);
            [maxRatio, cornerIdx] = max(temp); % naiive selection
            if maxRatio > tol*.5*pi
                p0 = epC(j,:);
                p1 = stC(mod(j+cornerIdx-1,nC)+1,:);
                p = bresenham(p0, p1); % draws a line between p0 and p1
                BL(sub2ind(dim,p(:,1),p(:,2))) = 3;
            end
        end
    end
return