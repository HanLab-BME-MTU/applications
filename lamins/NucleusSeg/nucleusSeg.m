function BW = nucleusSeg(I)
% This function executes the current process for nucleus segmentation given
% image I

    % Segment out the shadows
    BW1 = shadowSeg(I);
    % Threshold on distance to nearest neighbor
    [BW2, tRad] = clusterSeg(BW1);
    BW3 = getForeground(BW2, 8);
    BW3 = imfill(BW3, 'holes');
    %figure; imshow(BW1, [], 'InitialMagnification', 'fit')
    %figure; imshow(BW3, [], 'InitialMagnification', 'fit')
    
    tolC = pi/4;
    [CCcr, flag] = cncvCorner(BW3, tolC); % find concave corners
    
    if flag
        BW = BW1;
    else
        tolP = .9;
        BL = cutInteriorPolygon(CCcr, tRad, tolP);
        %figure; imshow(BL, [], 'InitialMagnification', 'fit')
    
        % Invert the new border and take the largest CC
        BW4 = getForeground(BL, 4);
        %figure; imshow(BW4, [], 'InitialMagnification', 'fit')
    
        BW = BW1;
        BW(~BW4) = 0;
        %figure; imshow(BW, [], 'InitialMagnification', 'fit')
    end
    
    BW = bwconvhull(BW);
    %figure; imshow(BW, [], 'InitialMagnification', 'fit')
end