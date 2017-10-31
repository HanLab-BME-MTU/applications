function [mask, vI, fvI, bfvI imgOut] = segCellLCH(img, varargin)
% simple script to pre-process the cells for Deep learning.


I = mat2gray(img);
vI = stdfilt(I);
fvI = imfilter(vI, fspecial('gaussian', 7,3));
bfvI = imbinarize(fvI, .02);
bfvI = imdilate(bfvI, strel('disk', 2));
mask = imclose(bfvI, strel('disk', 7));


% [mask vI fvI bfvI]=segCellLCH(imread('./14-May-2017_atcc_s06_t120_x998_y1586_t130_f6.png'));
% figure;subplot(2,2,1);imshow(bfvI,[]);subplot(2,2,2);imshow(vI,[]);subplot(2,2,3);imshow(fvI,[]);subplot(2,2,4);imshow(imgOut,[])

imgOut=mask.*I;
figure;subplot(2,2,1);imshow(bfvI,[]);subplot(2,2,2);imshow(vI,[]);subplot(2,2,3);imshow(fvI,[]);subplot(2,2,4);imshow(imgOut,[])


% Works roughly for 128x128 downsampling
% vI = stdfilt(I);
% fvI = imfilter(vI, fspecial('gaussian', 7,3));
% bfvI = imbinarize(fvI, .02);
% bfvI = imdilate(bfvI, strel('disk', 3));
% mask = imclose(bfvI, strel('disk', 7));