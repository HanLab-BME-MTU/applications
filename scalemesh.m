function mesh=scalemesh(geom,dim,varargin)
%SCALEMESH Scales geometry before meshing and then scales the mesh back.
%   MESH=SCALEMESH(GEOM,[SX,SY,SZ],...) scales the geometry by SX and SY
%   for 2D geometries and SX, SY and SZ for 3D geometries in the X-, Y-, and
%   Z-directions, respectively, before meshing it. After SZ
%   all property value pairs that are accepted by meshinit may be given.
%   SCALEMESH can be used if many elements are generated due to a
%   thin geometry or when the mesh generation fails due to large aspect
%   ratios in the geometry.
%
%   A geometry can be exported from the GUI and then meshed using SCALEMESH.
%   The resulting mesh can then be imported back into the GUI again.
%
%   Example:
%   geom=block3(1,1,0.1);
%   mesh=scalemesh(geom,[1,1,5],'hmax',0.5);
%   meshplot(mesh); axis equal;

%   T. Normark 11-29-01.
%   Copyright (c) 1994-2001 by COMSOL AB
%   $Revision: 1.14 $  $Date: 2001/03/13 15:53:52 $

d=size(dim);
if d(2)>2
    sx=dim(1);  % scaling x-direction
    sy=dim(2);  % scaling y-direction
    sz=dim(3);  % scaling z-direction
    rng=flgeomes(geom,1:flgeomnes(geom));
    geom1=scale(geom,sx,sy,sz);
    mesh=meshinit(geom1,varargin{:});
    mesh.p(1,:)=mesh.p(1,:)/sx;
    mesh.p(2,:)=mesh.p(2,:)/sy;
    mesh.p(3,:)=mesh.p(3,:)/sz;
    arcl=sqrt(sum((mesh.p(:,mesh.eg(2,:))-mesh.p(:,mesh.eg(1,:))).^2));
    neg=size(mesh.eg,2);
    ARCL=sparse(1:neg,mesh.eg(5,:),arcl);
    ARCL=cumsum(ARCL).*(ARCL>0);
    dist=sum(ARCL');
    mesh.eg(3,:)=dist-arcl;
    mesh.eg(4,:)=dist;
    endpoint=[find(diff(mesh.eg(5,:))) neg];
    mesh.eg(4,endpoint)=rng(2,:);
elseif d(2)>1
    sx=dim(1);  % scaling x-direction
    sy=dim(2);  % scaling y-direction
    rng=flgeomes(geom,1:flgeomnes(geom));
    geom1=scale(geom,sx,sy);
    mesh=meshinit(geom1,varargin{:});
    mesh.p(1,:)=mesh.p(1,:)/sx;
    mesh.p(2,:)=mesh.p(2,:)/sy;
else
    disp('You must specify a 2D or a 3D geometry')
end
