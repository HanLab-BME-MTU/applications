function [CC, flag] = cncvCorner(BW, tol)
% This function identifies  concave corners of a single boundary image with
% angles of 90+tol or smaller
%
%   B0 = counter clockwise binary image of the boundary
%   tol = angle tolerance
%   CC = connected component structure where the 1st CC is the 1st corner,
%   the 2nd CC is the pixels separating the 1st and 2nd corner, the 3rd CC
%   is the 2nd corner, and so on

    % Should check that B0 is counter clockwise
    B0 = bwboundaries(BW);
    B0 = B0{1};
    dim = size(BW);
    nB = size(B0,1);
    CC = struct('Connectivity', 8, 'ImageSize', dim, 'NumObjects', 0, ...
         'PixelIdxList', []);
    
    m = 5; % number of points away that we test
    B = [B0((end-m+1):end,:); B0]; % add dummy points at the beginning

    tthresh = pi/2 + tol;
    k = m+1; % find a starting point
    past = 1;
    present = 1;
    temp = zeros(2); % 2 vectors describing 3 points for each test
    % This loop finds start pixel of the first corner
    while past == present || present == 0
        
        if k == nB
            flag = 1;
            return;
        end
        
        k = k + 1;
        past = present;
        temp(1,:) = B(k-m,:) - B(k,:);
        temp(2,:) = B(k+m,:) - B(k,:);
        n1 = norm(temp(1,:));
        n2 = norm(temp(2,:));
        skew = det(temp)/(n1*n2);
        present = skew < 0 && acos(skew) > tthresh;
    end

    start = k - m;
    B = [B0(start:end,:); B0(1:(start-1),:)];
    B = [B((end-m+1):end,:); B; B(1:m,:)]; % add dummy points on each end
    Blin = sub2ind(dim, B(:,1), B(:,2));
    
    past = 1;
    startIdx = 1+m; % starting index of the current connected component
    num = 0; % number of CC's
    for j = startIdx:(m+nB)
        temp(1,:) = B(j-m,:) - B(j,:);
        temp(2,:) = B(j+m,:) - B(j,:);
        n1 = norm(temp(1,:));
        n2 = norm(temp(2,:));
        skew = det(temp)/(n1*n2);
        present = skew < 0 && acos(skew) > tthresh; % concave or not
        if past ~= present
            num = num + 1;
            CC.PixelIdxList{num} = Blin(startIdx:(j-1));
            startIdx = j;
        end
        past = present;
    end
    % Add the last component 
    num = num + 1;
    CC.PixelIdxList{num} = Blin(startIdx:j);
    CC.NumObjects = num;
    
    flag = 0;
end