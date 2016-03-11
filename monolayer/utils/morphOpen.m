function [openIm] = morphOpen(im, maskSize)
imDil = dilateAssaf(im,maskSize);
openIm = erodeAssaf(imDil,maskSize);
end