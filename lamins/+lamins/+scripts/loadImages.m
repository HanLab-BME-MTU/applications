% load('project/ag_080712wt_Reconstructed 3/ag_080712wt_Reconstructed 3.mat');
% get the bioformats reader
reader = MD2.getReader();
% retrieve number of channels and z slices
nc = reader.getSizeC();
nz = reader.getSizeZ();

% load all the images into a cell array nc x nz
I = cell(nc,nz);
for c = 1:nc
    for z = 1:nz
        I{c,z} = reader.loadImage(c,1,z);
    end
end

maxima = cellfun(@(I) max(I(:)),I);
minima = cellfun(@(I) min(I(:)),I);

% thumbs = vertcat(imadjust(horzcat(I{1,:})),imadjust(horzcat(I{2,:})),imadjust(horzcat(I{3,:})),imadjust(horzcat(I{4,:})));
% smallthumbs = imresize(thumbs,0.25);
% 
% imwrite(im2uint8(thumbs),jet(256),'thumbs2.png');
% imwrite(im2uint8(imresize(thumbs,0.25)),jet(256),'smallthumbs2.png');

for i=1:20; for j=1:20; m(i,j) = corr2(I{2,i},I{4,j}); end; end;
figure; imagesc(m);
