function BW = nucleusSeg(I)
% This function executes the current process for nucleus segmentation given
% image I

    % scale pixel intensities to avoid loss of precision
    Isc = imadjust(I, stretchlim(I,0), []);
    BW1 = shadowSeg(Isc); % Segment out the shadows
    figure; imshow(BW1, [], 'InitialMagnification', 'fit')

    BW1f = getForeground(~imfill(BW1, 'holes'), 8);
    tolI = 1.05;
    figure; imshow(BW1f, [], 'InitialMagnification', 'fit')
    
    % check if BW1f larger than the original mask and has tolI increase in
    % mean pixel intensity - should idenitfy if 
    if 0==1 %sum(BW1f(:)) > sum(BW1(:)) && mean(Isc(BW1f(:))) > tolI*mean(Isc(BW1(:)))
        BW = BW1; 
        BW(~BW1f) = 1; % delete background
        BW = getForeground(BW, 8);
        BW = imfill(BW, 'holes'); 
        %figure; imshow(BW, [], 'InitialMagnification', 'fit')
        BW = bwconvhull(BW)*2;
        return
    end
     
    % Threshold on distance to nearest neighbor
    BW2 = clusterSeg(BW1);
    BW3 = getForeground(BW2, 8);
    BW3 = imfill(BW3, 'holes');
    figure; imshow(BW2, [], 'InitialMagnification', 'fit')
    figure; imshow(BW3, [], 'InitialMagnification', 'fit')
    
    tolC = pi/4;
    [CCcr, flag] = cncvCorner(BW3, tolC); % find concave corners
    
    if flag
        BW = BW1;
    else
        tolP = .9;
        distMax = 200;
        BL = cutInteriorPolygon(CCcr, distMax, tolP);
        figure; imshow(BL, [], 'InitialMagnification', 'fit')
    
        % Invert the new border and take the largest CC
        BW4 = getForeground(BL, 4);
        figure; imshow(BW4, [], 'InitialMagnification', 'fit')
    
        BW = BW1;
        BW(~BW4) = 0;
        figure; imshow(BW, [], 'InitialMagnification', 'fit')
    end
    
    BW = bwconvhull(BW)*1.0;
    figure; imshow(BW, [], 'InitialMagnification', 'fit')
end