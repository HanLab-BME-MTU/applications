function gauss=gauss_3D(amp,sigma,cent,XX,YY,ZZ)
% GAUSS3D	compute values of gaussian in 3D 
%
%    gauss=gauss3D(sigma,center,X,Y,Z)
%
%    INPUT: amp    (scalar) amplitude
%           sigma  of gauss mask [sigmaX sigmaY sigmaZ]
%           cent   3D vector with center position [0 0 0] 
%           Y,Y,Z  vectors of equal size containg the positions of for the gaussian
%
%    OUTPUT: gauss   vector values of gaussian @ X,Y,Z

% c: 12/08/02 dT


x=(XX-cent(1))./sigma(1);
y=(YY-cent(2))./sigma(2);
z=(ZZ-cent(3))./sigma(3);
ex=exp(-1/2*(x.^2));
ey=exp(-1/2*(y.^2));
ez=exp(-1/2*(z.^2));

gauss=amp*ex.*ey.*ez;    
