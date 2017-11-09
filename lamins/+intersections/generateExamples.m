function [ I, Ig2, Ig3 ] = generateExamples( varargin )
%generateExamples Summary of this function goes here
%   Detailed explanation goes here

ip = inputParser;
ip.addParameter('size',[101 101]);
ip.addParameter('sigma',2);
ip.addParameter('padding',[77 77]);
ip.addParameter('noise',1e-1);
ip.addParameter('singleFile',true);
ip.addParameter('mode','add');
ip.addParameter('type','twoRadial');
ip.parse(varargin{:});

% if(nargin < 1)
%     sz = [101 101];
% end
% if(nargin < 2)
%     sigma = 2;
% end
% if(nargin < 3)
%     padding = [77 77];
% end
% if(nargin < 4)
%     noise = 1e-1;
% %     noise = 0;
% end

singleFile = true;
mode = 'add';

if(~singleFile)
    mkdir('orig');
    mkdir('g2');
    mkdir('g3');
else
    if(exist('orig.tiff','file'))
        movefile('orig.tiff',['orig_' datestr(datetime('now'),30) '.tif']);
        movefile('g2.tiff',['g2' datestr(datetime('now'),30) '.tif']);
        movefile('g3.tiff',['g3_' datestr(datetime('now'),30) '.tif']);
    end
end

switch(ip.Results.type)
    case 'twoCircles'
        coords = orientationSpace.getFrequencySpaceCoordinates(ip.Results.size);
        r = fftshift(coords.f);
        C = abs(r-0.25) < 0.003;
end



for ii = 1:180
    switch(ip.Results.type)
    case 'twoRadial'
        % Two Radial Lines
        I = intersections.drawRadialLines([0 ii]/180*pi,ip.Results.size,ip.Results.mode);
    case 'twoSymmetric'
        % Two symmetric lines
        I = intersections.drawTwoLines([0 ii]/180*pi,ip.Results.size,ip.Results.mode);
    case 'threeRadial'
        % Three Radial Lines
        I = intersections.drawRadialLines([0 ii -ii]/180*pi,ip.Results.size,ip.Results.mode);
    case 'twoCircles'
        switch(mode)
            case 'or'
                I = double(C | circshift(C,ii));
            case 'add'
                I = C + circshift(C,ii);
            otherwise
                error('Unknown mode');
        end
    otherwise
        error('Unknown Type');
    end
    I = padarray(I,ip.Results.padding);
    I = mat2gray(I);
    
    Ig2 = imnoise(mat2gray(imgaussfilt(I,2)),'gaussian',0,ip.Results.noise);
    
    Ig3 = imnoise(mat2gray(imgaussfilt(I,3)),'gaussian',0,ip.Results.noise);
    if(~singleFile)
        imwrite(I,sprintf('orig/%02d.png',ii));
        imwrite(Ig2,sprintf('g2/%02d_g2.png',ii));
        imwrite(Ig3,sprintf('g3/%02d_g3.png',ii));
    else
        imwrite(I,'orig.tiff','WriteMode','append');
        imwrite(Ig2,'g2.tiff','WriteMode','append');
        imwrite(Ig3,'g3.tiff','WriteMode','append');
    end
end

end