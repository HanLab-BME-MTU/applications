function [closeIm] = morphClose(im, maskSize)
imEr = erodeAssaf(im,maskSize);
closeIm = dilateAssaf(imEr,maskSize);
end