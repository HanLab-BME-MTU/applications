%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  compute the coefficient function, mu(x,y) of the wave equation
%
%          u_xx + u_yy - (1/mu)u_tt = 0
%
%     x : column vector.
%     y : row vector.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function f = f_mu(t,x,y,A,xc,yc,wx,wy,rotXY)
function f = f_mu(t,x,y,A)

f = A + zeros(size(x));
%f = ones(size(x));
return;

x_rot = cos(rotXY)*(x-xc) - sin(rotXY)*(y-yc);
y_rot = sin(rotXY)*(x-xc) + cos(rotXY)*(y-yc);

x2 = (x_rot/wx).^2; 
y2 = (y_rot/wy).^2;

f = A*exp(-(x2+y2)/2);

f = 1+f;

