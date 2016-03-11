function W = wrap_affine(im, T)

tform = maketform('affine',T');

W = imtransform(im, tform, 'YData', [1,size(im,1)],'XData', [1,size(im,2)]);

end