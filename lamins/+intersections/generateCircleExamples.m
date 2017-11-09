function [ I, Ig2, Ig3 ] = generateCircleExamples( sz, sigma, padding, noise, mode )
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
if(nargin < 5)
    mode = 1;
else
    switch(mode)
        case 'or'
            mode = 1;
        case 'add'
            mode = 2;
        otherwise
            if(~isnumeric(mode))
                warning('intersections.drawRadialLines:Did not understand mode parmeter');
                mode = 1;
            end
    end
end

singleFile = true;

if(~singleFile)
    mkdir('orig');
    mkdir('g2');
    mkdir('g3');
end

coords = orientationSpace.getFrequencySpaceCoordinates(201);
r = fftshift(coords.f);
C = abs(r-0.25) < 0.003;

for ii = 1:100
    switch(mode)
        case 1
            I = double(C | circshift(C,ii));
        case 2
            C = C + circshift(C,ii);
    end
%     I = intersections.drawRadialLines([0 ii]/180*pi);
%     I = intersections.drawTwoLines([0 ii]/180*pi);
    I = padarray(I,padding);
        
    Ig2 = imnoise(mat2gray(imgaussfilt(I,2)),'gaussian',0,noise);
    
    Ig3 = imnoise(mat2gray(imgaussfilt(I,3)),'gaussian',0,noise);
    if(~singleFile)
        imwrite(I,sprintf('orig/%02d.png',ii));
        imwrite(Ig2,sprintf('g2/%02d_g2.png',ii));
        imwrite(Ig3,sprintf('g3/%02d_g3.png',ii));
    else
        imwrite(I,'orig.tiff','WriteMode','append');
        imwrite(I,'g2.tiff','WriteMode','append');
        imwrite(I,'g3.tiff','WriteMode','append');
    end
end

end
