function fx = myoDragFx(t,x,y)

fx = zeros(size(x));
I = find(x>-.5 & x<.5);
fx(I) = 1./(x(I).^2-.25).^2;
fx(I) = -8e8*x(I).*fx(I).*exp(-fx(I));

