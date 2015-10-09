function [ filterKernel, radialFilter, angularFilter ] = steerableVanGinkelKernel( f_c, b_f, K, angle, N)
%steerableVanGinkelKernel
%
% Based on the thesis by Michael van Ginkel. Chapter 3
% "Image Analysis using Orientation Sapce based on Steerable Filters".
% Delft University of Technology. October 7th, 2002.
%
% f_c: maximum frequency for the radial filter
% b_f: frequency bandwidth for the radial filter
% K: number of rotation angles through 360 degrees

% Mark Kittisopikul, August 22nd, 2015
% Jaqaman Lab
% UT Southwestern

    import vanGinkel.*;

if(nargin < 3)
    K = 36;
end
if(nargin < 4)
    angle = 0;
end
if(nargin < 5)
    N = 512;
end

[X, Y] = meshgrid(-N:(N-1));
[theta, rho] = cart2pol(X,Y);

theta = fftshift(theta);
rho = fftshift(rho);


% theta_shifted = acos(cos(theta - angle + pi));
% theta = acos(cos(theta - angle));
% 
% theta = theta - angle;
theta = bsxfun(@minus,theta,shiftdim(angle(:),-2));
% theta_shifted = theta + pi;

theta = mod(theta+pi,2*pi)-pi;
% theta = atan2(sin(theta),cos(theta));
% theta_shifted = atan2(sin(theta_shifted),cos(theta_shifted));




f = rho / N / 2;

%% Radial part
% compute radial order, f_c = sqrt(K_f * b_f^2)
K_f = (f_c / b_f)^2;

% scale frequency
f_s = f / f_c;

% Equation 3.11
% Note -(f^2 - f_c^2)/(2*b_f^2) = (1 - (f/f_c)^2)/(2* b_f^2/f_c^2)
%                               = (1 - (f/f_c)^2)/(2 / K_f)
radialFilter = f_s.^K_f .* exp((1 - f_s.^2)*K_f/2);
% radialFilter2 = f_s^K_f .* exp(-(f.^2-f_c.^2)/2/b_f.^2);
% assertEqual(radialFilter,radialFilter2);

%% Angular part
% scale the angle

s_a = pi/(2*K+1);

theta_s = theta / s_a;
% theta_shifed_s = theta_shifted / s_a;

angularFilter = 2*exp(-theta_s.^2/2);
% angularFilter_shifted = 2*exp(-theta_shifed_s.^2/2);
angularFilter_shifted = angularFilter([1 end:-1:2],[1 end:-1:2],:);

filterKernel = bsxfun(@times,radialFilter .* 0.5,angularFilter + angularFilter_shifted);

% could we simplify this? neg/pos dividing line does not have to rotate
posMask = abs(theta) < pi/2;
negMask = ~posMask;
filterKernel(posMask) = filterKernel(posMask) + filterKernel(posMask)*1j;
filterKernel(negMask) = filterKernel(negMask) - filterKernel(negMask)*1j;


% shift = exp(1i*2*pi*(X/2+Y/2));
% filterKernel = radialFilter .* angularFilter;
% filterKernel = radialFilter .* angularFilter .* 2 .* ( abs(theta) < pi/2 );
% filterKernel = radialFilter .* angularFilter .*exp(1i*2*pi*(X/2+Y/2));
% filterKernel = radialFilter .* angularFilter + 1i * radialFilter .* angularFilter .* (-1).^(abs(theta) >= pi/2);
% filterKernel = filterKernel .* (f > 0);



end

