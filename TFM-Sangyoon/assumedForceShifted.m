function force = assumedForceShifted(j,x,y,xshift,yshift,wx,wy)
% This function assumedForceShifted takes grid x and y and make gaussian
% distributed force field of which sources are xshift and yshift.
% input     :   x           grid of x coordinates
%               y           grid of y coordinates
%               xshift      x value of point source of force
%               yshift      y value of point source of force
%               wx          x component of force orientation
%               wy          y component of force orientation
% output    :   force       x or y grid of force distribution (if j=1 or 2,
% respectively)
%              
% Sangyoon Han Jan 2013
std = sqrt(wx^2+wy^2);
if j==1
    force=wx*(heaviside(2.1-sqrt((x-xshift).^2+(y-yshift).^2)).*...
        (exp(-((x-xshift).^2+(y-yshift).^2)/(2*std^2))));
elseif j==2
    force=wy*(heaviside(2.1-sqrt((x-xshift).^2+(y-yshift).^2)).*...
        (exp(-((x-xshift).^2+(y-yshift).^2)/(2*std^2))));
else
    error('please input 1 or 2 for j');
end

return

%test
x=0:1:20;
xshift = 7;
std = 1.4;
fd = heaviside(3.1-sqrt((x-xshift).^2)).*...
    exp(-((x-xshift).^2)/(2*std^2));
figure,plot(x,fd)
     