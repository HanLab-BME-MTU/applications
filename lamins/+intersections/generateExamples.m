function [ I, Ig2, Ig3 ] = generateExamples( sz, sigma, padding, noise )
%generateExamples Summary of this function goes here
%   Detailed explanation goes here

if(nargin < 1)
    sz = [101 101];
end
if(nargin < 2)
    sigma = 2;
end
if(nargin < 3)
    padding = [77 77];
end
if(nargin < 4)
    noise = 1e-2;
end

mkdir('orig');
mkdir('g2');
mkdir('g3');



for ii = 1:180
%     I = intersections.drawRadialLines([0 ii]/180*pi);
    I = intersections.drawTwoLines([0 ii]/180*pi);
    I = padarray(I,padding);
    imwrite(I,sprintf('orig/%02d.png',ii));
    Ig2 = imnoise(mat2gray(imgaussfilt(I,2)),'gaussian',0,noise);
    imwrite(Ig2,sprintf('g2/%02d_g2.png',ii));
    Ig3 = imnoise(mat2gray(imgaussfilt(I,3)),'gaussian',0,noise);
    imwrite(Ig3,sprintf('g3/%02d_g3.png',ii));
end

end