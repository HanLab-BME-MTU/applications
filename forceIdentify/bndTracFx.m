function fx = bndTracFx(x,nx,A)

wx = 300;
xc = 0;
x2 = ((x-xc)/wx).^2;

fx = -A*nx;
%fx = -A*exp(-x2/2).*nx;

