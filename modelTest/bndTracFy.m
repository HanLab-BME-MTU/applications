function fy = bndTracFy(x,ny,A)

wx = 300;
xc = 0;
x2 = ((x-xc)/wx).^2;

fy = -A*ny;
%fy = -A*exp(-x2/2).*ny;

