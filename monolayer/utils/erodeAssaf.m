function [Aer] = erodeAssaf(A,maskSize)

if nargin < 2
    maskSize = 10;
end

debug = 0;

mask = ones(maskSize);
se1 = strel(mask);

Aer = imerode(A,se1);

if debug
    figure; colormap(gray);
    subplot(2,2,1); imagesc(A); title('org');
    subplot(2,2,2); imagesc(Aer); title(sprintf('erode %d pixels',maskSize));    
end