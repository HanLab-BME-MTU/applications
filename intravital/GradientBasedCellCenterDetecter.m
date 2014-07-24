function [ imCenterProbability ] = GradientBasedCellCenterDetecter( im, radrange, mask )
% Tries to detect/enhance the core regions of cells using gradient based
% hough voting
%
    imsize = size(im);
    imSmooth = filterGauss2D(im,3);
    [ gx, gy ] = gradient( imSmooth );
    gmag = sqrt( gx.^2 + gy.^2 );        
    maxgrad = max(gmag(:));
    gradThresh = 0.05 * maxgrad;
    
    imCenterProbability = zeros( size(im) );

    if ~exist( 'mask', 'var' )       
        mask = ones( size(im) );
    end
    
    pixind = find( mask > 0 );
    [yind,xind] = ind2sub( size(mask), pixind );
    
    for i = 1:numel(pixind)        
        
        cur_gmag = gmag( pixind(i) );
                
        if cur_gmag <= 0
            continue;
        end
        
        gvec = [ gx(pixind(i)), gy(pixind(i)) ] / cur_gmag;
        
        ptCur = [xind(i), yind(i)];
        ptCenter = ptCur + 0.5 * gvec * range(radrange);
        
        sigma = range(radrange) / 2.5;
        v = @(pt) ( (sigma^2 * 2 * pi)^-(numel(pt)/2) * exp( -(norm(pt - ptCenter) / sigma)^2 )  );
        
        ptRayList = repmat(ptCur, [range(radrange) + 1, 1]) + (radrange(1):radrange(2))' * gvec;
        
        for r = radrange(1):radrange(2)

            ptNew = round( ptCur + r * gvec );

            if any(ptNew < 1) || any(ptNew > imsize(2:-1:1))
                continue;
            else
                imCenterProbability(ptNew(2),ptNew(1)) = imCenterProbability(ptNew(2),ptNew(1)) + cur_gmag * v(ptNew) / maxgrad;
            end
            
        end        
    end
    
end