function BW = nucleusSeg(I)
% This function executes the current process for nucleus segmentation given
% image I

    % scale pixel intensities to avoid loss of precision
    Isc = imadjust(I, stretchlim(I,0), []);
    BW1 = shadowSeg(Isc); % Segment out the shadows
    %figure; imshow(BW1, [], 'InitialMagnification', 'fit')

    se = strel('disk', 10);
    BW1c = imclose(BW1, se);
    BW1f = getForeground(~imfill(BW1c, 'holes'), 8);
    tolI = 1.10;
    %figure; imshow(BW1f, [], 'InitialMagnification', 'fit')
    
    
    backI = double(Isc(~BW1f(:))); % background intensities of filled image
    [N, X] = optimalHistogram(backI);
    [dummy, ti] = cutFirstHistMode(N, X, 0); % find the right hand tail
    lol = (mean(backI(backI>ti)) - ti)/65535; % look for bright pixels
    
    
    % check if BW1f larger than the original mask and has tolI increase in
    % mean pixel intensity - should idenitfy if 
    if sum(BW1f(:)) > sum(BW1(:)) && mean(Isc(BW1f(:))) > tolI*mean(Isc(BW1(:)))
        BW = BW1;
        BW(~BW1f) = 1; % delete background
        BW = getForeground(BW, 8);
        BW = imfill(BW, 'holes'); 
        BW = bwconvhull(BW)*lol;
        %figure; imshow(BW, [], 'InitialMagnification', 'fit')
        return
    end
     
    % Threshold on distance to nearest neighbor
    BW2 = distSeg(BW1);
    BW3 = getForeground(BW2, 8);
    BW3 = imfill(BW3, 'holes');
    %figure; imshow(BW2, [], 'InitialMagnification', 'fit')
    %figure; imshow(BW3, [], 'InitialMagnification', 'fit')
    
    tolC = pi/4;
    [CCcr, flag] = cncvCorner(BW3, tolC); % find concave cornershel
    
    if flag
        BW = BW1;
    else
        tolP = .9;
        distMax = 150;
        BL = cutInteriorPolygon(CCcr, distMax, tolP);
        %figure; imshow(BL, [], 'InitialMagnification', 'fit')
    
        % Invert the new border and take the largest CC
        BW4 = getForeground(BL, 4);
        %figure; imshow(BW4, [], 'InitialMagnification', 'fit')
    
        BW = BW1;
        BW(~BW4) = 0;
        %figure; imshow(BW, [], 'InitialMagnification', 'fit')
    end
    
    BW = bwconvhull(BW)*lol;
    %figure; imshow(BW, [], 'InitialMagnification', 'fit')
end