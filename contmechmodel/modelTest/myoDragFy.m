function fy = myoDragFy(x,y,A)

fy = zeros(size(y));

indL = find(x<0);
indR = find(x>0);

wxL = 80;
xcL = -350;
x2L = ((x(indL)-xcL)/wxL).^2;

fy(indL) = A*exp(-x2L/2);

wxR = 80;
xcR = 300;
x2R = ((x(indR)-xcR)/wxR).^2;

fy(indR) = A*exp(-x2R/2);
