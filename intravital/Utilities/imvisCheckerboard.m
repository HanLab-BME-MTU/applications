function [r] = imvisCheckerboard(im1,im2,blockSize,weights,color)
% color is not working yet

if(~exist('weights','var') || isempty(weights))
    weights = [1 1];
end
if(~exist('color','var'))
    color = 0;
end

[m n] = size(im1);

[X Y] = meshgrid(1:n,1:m);

cx = floor(X./blockSize);
cx = rem(cx,2);

cy = floor(Y./blockSize);
cy = rem(cy,2);

c = xor(cx,cy);

ind0 = find(c==0);
ind1 = find(c==1);

if ~color
    r = double(c);
    r(ind0) = weights(1)*im1(ind0);
    r(ind1) = weights(2)*im2(ind1);
else
    r1 = 0*double(c);
    r1(ind0) = weights(1)*im1(ind0);   
    r2 = 0*double(c);
    r2(ind1) = weights(2)*im2(ind1);
    r(:,:,1) = r1;
    r(:,:,2) = r2;
    r(:,:,3) = c*0;    
end
