function [ I3, L ] = sesame( I, PSF, winding )
%sesame Do simulated structured illumination microscopy
%
% I, true image source with no PSF applied
% PSF, point spread function on object space
% winding, number of oscillations of structured illumination the field
%
% I3 

    if(nargin < 2)
        winding = 16;
    end

    % Do structured illumination
    [X,Y] = meshgrid(1:size(I,1),1:size(I,2));
    
    % angles to rotate illumination pattern
    % 2pi/3 4pi/3 0
    theta = [1 2 3]*2*pi/3;
    
    % X = sin(theta).*X - cos(theta).*Y;

    I3 = cell(3,3);
    L = cell(3,3);
    Xr = cell(1,3);
    Yr = cell(1,3);
    for angle_j=1:3
        Xr{angle_j} = cos(theta(angle_j)).*X - sin(theta(angle_j)).*Y;
        Yr{angle_j} = sin(theta(angle_j)).*X + cos(theta(angle_j)).*Y;
        for phase_i=1:3
            L{phase_i,angle_j} = (1+cos(Xr{angle_j}*2*pi/winding+phase_i*2*pi/3))/2;
            I3{phase_i,angle_j} = imfilter(I.*L{phase_i,angle_j},PSF);
        end;
    end

    figure;
    imshowpair(abs(fftshift(fft2(L{1}))),MTF);


end

