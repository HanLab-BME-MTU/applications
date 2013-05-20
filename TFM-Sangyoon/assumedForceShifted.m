function force = assumedForceShifted(j,x,y,xshift,yshift,wx,wy,forceType,FAsize)
% This function assumedForceShifted takes grid x and y and make gaussian
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
% output    :   force       x or y grid of force distribution (if j=1 or 2,
% respectively)
%              
% Sangyoon Han Jan 2013
if strcmp(FAsize,'largeFA')
    std = 4.5;
    adh_r = 4.1; % adhesion diameter in pixel
else
    std = 1.5;
    adh_r = 1.1;
end
if j==1
    switch(forceType)
        case 'groupForce'
            force=wx*(heaviside(adh_r-sqrt((x-xshift).^2+(y-yshift).^2)).*...
                (exp(-((x-xshift).^2+(y-yshift).^2)/(2*std^2))));
        case 'pointForce'
            force=wx*(x==xshift || y==yshift);
        case 'smoothForce'
            force=wx*(exp(-((x-xshift).^2+(y-yshift).^2)/(2*std^2)));
    end
elseif j==2
    switch(forceType)
        case 'groupForce'
            force=wy*(heaviside(adh_r-sqrt((x-xshift).^2+(y-yshift).^2)).*...
                (exp(-((x-xshift).^2+(y-yshift).^2)/(2*std^2))));
        case 'pointForce'
            force=wy*(x==xshift || y==yshift);
        case 'smoothForce'
            force=wy*(exp(-((x-xshift).^2+(y-yshift).^2)/(2*std^2)));
    end
else
    error('please input 1 or 2 for j');
end

return

% test
x = 0:1:20;
xshift = 7;
std = 2;
w = 1300;
fd = w*heaviside(2.1-sqrt((x-xshift).^2)).*...
    exp(-((x-xshift).^2)/(2*std^2));
% fd = exp(-((x-xshift).^2)/(std^2));
% fd = x==xshift;
figure,plot(x,fd)
     