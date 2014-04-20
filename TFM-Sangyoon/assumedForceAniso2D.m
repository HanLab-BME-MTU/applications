function force = assumedForceAniso2D(j,x,y,xshift,yshift,wx,wy,d1,d2,forceType)
% assumedForceAniso2D takes grid x and y and make anisotropic gaussian
% distributed force field of which sources are xshift and yshift.
% input     :   x           grid of x coordinates
%               y           grid of y coordinates
%               xshift      x value of point source of force
%               yshift      y value of point source of force
%               wx          x component of force orientation
%               wy          y component of force orientation
%               forceType   type of force
%               ('pointForce','groupForce', or 'smoothForce')
%               FAsize      'largeFA' or 'smallFA'
%               dx,dy       diameter of FA in pixel in directions of normal
%               and tangential to (wx, wy) assuming that the length of FA
%               is in parallel with the direction of force
%               
% output    :   force       x or y grid of force distribution (if j=1 or 2,
% respectively)
%              
% Sangyoon Han Jan 2013
% if strcmp(FAsize,'largeFA')
%     std = 4.5;
%     adh_r = 4.1; % adhesion diameter in pixel
% else
%     std = 1.5;
%     adh_r = 1.1;
% end

% orientation of force
theta = atan2(wy,wx);
% if theta > pi/2

stdx = d1/2*1.1;
stdy = d2/2*1.1;
adh_rx = d1/2;
adh_ry = d2/2;
if j==1
    switch(forceType)
        case 'groupForce'
            force = anisoGaussian2DMatrix(xshift, yshift, stdy, stdx, theta, x, y);
            % cutting out outside of area
            force = wx*aniso2DMask(xshift, yshift, adh_ry, adh_rx, theta, x, y).*force;
        case 'pointForce'
            force=wx*(x==xshift || y==yshift);
        case 'smoothForce'
            force = anisoGaussian2DMatrix(xshift, yshift, stdy, stdx, theta, x, y);
            % cutting out outside of area
            force = wx*force;
    end
elseif j==2
    switch(forceType)
        case 'groupForce'
            force = anisoGaussian2DMatrix(xshift, yshift, stdy, stdx, theta, x, y);
            % cutting out outside of area
            force = wy*aniso2DMask(xshift, yshift, adh_ry, adh_rx, theta, x, y).*force;
        case 'pointForce'
            force=wy*(x==xshift || y==yshift);
        case 'smoothForce'
            force = anisoGaussian2DMatrix(xshift, yshift, stdy, stdx, theta, x, y);
            % cutting out outside of area
            force = wy*force;
    end
else
    error('please input 1 or 2 for j');
end

return

% test in 1D
x = 0:1:20;
xshift = 10;
std1d = 3.45;
w = 200;
fd = w*heaviside(std1d*1.1-sqrt((x-xshift).^2)).*exp(-((x-xshift).^2)/(2*std1d^2));
% fd = exp(-((x-xshift).^2)/(std^2));
% fd = x==xshift;
x_um = (x-xshift)*0.072;
figure,plot(x_um,fd)
     
% test in 2D
[x,y] = meshgrid(1:30,1:30);
xshift = 7;
yshift = 9;
std1d = 2;
wx = 10;
wy = 100;
amp = (wx^2+wy^2)^.5;
d1 = 2;
d2 = 10;
forceType = 'groupForce';
force = assumedForceAniso2D(1,x,y,xshift,yshift,wx,wy,d1,d2,forceType);


