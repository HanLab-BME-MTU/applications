function fy = myoDFy(x,y,A)

wx = 80;
wy = 2.5;
xc = -100;
yc = 70;
x2 = ((x-xc)/wx).^2;
y2 = ((y-yc)/wy).^2;

fy = A*exp(-(x2+y2)/2);
