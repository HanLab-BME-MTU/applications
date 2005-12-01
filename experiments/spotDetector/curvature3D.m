function k=curvature3D(img,pos)
% CURVATURE3D compute the curvature k at pos
%
% SYNOPSIS  k=curvature3D(img,pnt)
%
% INPUT img   : 3D data
%       pos   : curvature at pos [x y z]
% OUTPUT  k   : curvature k

% c: 17/08/01	dT
delta=1;

[FX,FY,FZ] = gradient(img,delta);
[FXX,FXY,FXZ] = gradient(FX,delta);
[FYX,FYY,FYZ] = gradient(FY,delta);
[FZX,FZY,FZZ] = gradient(FZ,delta);
HE=[FXX(pos(1),pos(2),pos(3)) FXY(pos(1),pos(2),pos(3)) FXZ(pos(1),pos(2),pos(3));...
        FYX(pos(1),pos(2),pos(3)) FYY(pos(1),pos(2),pos(3)) FYZ(pos(1),pos(2),pos(3));...
        FZX(pos(1),pos(2),pos(3)) FZY(pos(1),pos(2),pos(3)) FZZ(pos(1),pos(2),pos(3))];
k=det(HE);
k=det(HE);