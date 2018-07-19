function [sumVec,sumMag,method,inpos,invec]=integrateForceField(fpos,fvec,bwMask,pixSize_mu,gridSpacing)
% forceField (fpos,fvec) and bwMask have to be in the same coordinate
% system!
if nargin > 4 && ~isempty(gridSpacing)
    method='noIntp';
else
    method='pixIntp';
    gridSpacing=1;
end

size1= max(fliplr(fpos));
size2=size(bwMask);
ImgSize=max(vertcat(size1,size2));

% extend the mask if smaller than the image:
if numel(bwMask)<prod(ImgSize)
    bwMask(ImgSize(1),ImgSize(2))=0;
end
%check_mat1= logical(cellMask);

if strcmp(method,'noIntp')
    if sum(sum(abs(fpos-round(fpos))))==0
        [inpos,invec]=findVectorFieldInMask(fpos,fvec,bwMask);
    else
        display('Use inpolygon');
        checkVector = inpolygon(fpos(:,1),fpos(:,2),bwMask(:,1),bwMask(:,2));
        
        inpos=fpos(checkVector,:);
        invec=fvec(checkVector,:);
    end
else
    % interpolate the force field on all pixel positions:
    fx_intp = TriScatteredInterp(fpos(:,1),fpos(:,2),fvec(:,1));
    fy_intp = TriScatteredInterp(fpos(:,1),fpos(:,2),fvec(:,2));
    
    % convert the bw-mask into x-y-positions:
    [ypos,xpos] = ind2sub(ImgSize,find(bwMask));
    
    % interpolate the force field on these positions:
    inpos = horzcat(xpos,ypos);
    invec(:,1)= fx_intp(xpos,ypos);
    invec(:,2)= fy_intp(xpos,ypos);
end

factor_Pa2nN=gridSpacing^2*pixSize_mu^2*10^(-3);

% This is the sum of all forces. In
% constrForceField.cell.stats.resForce.vec we store however the reaction
% force which '-' this value:
sumVec = sum(invec,1)*factor_Pa2nN;

sumMag = sum(sqrt(sum(invec.^2,2)))*factor_Pa2nN;
