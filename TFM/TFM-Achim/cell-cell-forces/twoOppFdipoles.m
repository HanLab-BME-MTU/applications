function [fx fy]=twoOppFdipoles(x,y,origin,s,r,d,sep)
% s: stress on the adhesions.
% r: radius of the adhesions.
% d: length of the dipole (distance between the two focal adhesions of one
%    dipole).
% sep: separation of the twod dipoles.

%           origin
% a1------a2  ||  a3------a4

% the centers are
c1=origin;

c1(1)=origin(1)-r-sep/2-d;
c2(1)=origin(1)-r-sep/2;
c3(1)=origin(1)+r+sep/2;
c4(1)=origin(1)+r+sep/2+d;

c1(2)=origin(2);
c2(2)=origin(2);
c3(2)=origin(2);
c4(2)=origin(2);

% The stress is zero everywhere:
fx=zeros(size(x));
fy=zeros(size(x));

% Set it to nonzero values only where the adhesions are:
check1=sqrt((x-c1(1)).^2+(y-c1(2)).^2)<r;
fx(check1)=s;

check2=sqrt((x-c2(1)).^2+(y-c2(2)).^2)<r;
fx(check2)=-s;

check3=sqrt((x-c3(1)).^2+(y-c3(2)).^2)<r;
fx(check3)=s;

check4=sqrt((x-c4(1)).^2+(y-c4(2)).^2)<r;
fx(check4)=-s;

return;





