function fx = myoDragFx(x,y,A)

fx = zeros(size(x));

indL = find(x<0);
indR = find(x>0);

wxL = 120;
xcL = -350;
x2L = ((x(indL)-xcL)/wxL).^2;

fx(indL) = A*exp(-x2L/2);

wxR = 120;
xcR = 300;
x2R = ((x(indR)-xcR)/wxR).^2;

fx(indR) = -A*exp(-x2R/2);
