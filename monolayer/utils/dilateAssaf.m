function [Adil] = dilateAssaf(A,maskSize)

if nargin < 2
    maskSize = 10;
end

debug = 0;

mask = ones(maskSize);
se1 = strel(mask);

Adil = imdilate(A,se1);

if debug
    figure; colormap(gray);
    subplot(2,2,1); imagesc(A); title('org');
    subplot(2,2,2); imagesc(Adil); title(sprintf('dilated %d pixels',maskSize));    
end