displField(1).plane=strain_z(1).level;

%displField(frame).plane(layer).pos;
%displField(frame).plane(layer).vec;
%displField(frame).plane(layer).z; z should be in pixel, can be real value.

% This is to make the displacement field sparse (to many points for my
% laptop):
% displField(1).pos(1:2:end,:)=[];
% displField(1).vec(1:2:end,:)=[];
% 
% displField(1).pos(1:2:end,:)=[];
% displField(1).vec(1:2:end,:)=[];
% 
% displField(1).pos(1:2:end,:)=[];
% displField(1).vec(1:2:end,:)=[];

% Take the first frame, this should be a loop in the future.
frame=1;
layer=1;

% First we get a regular grid with the right density of meshpoints. The
% created mesh will be centered in the field of view. Here, only the first
% frame will be used:
[xvec,yvec]=createHexGridFromDisplField(displField(frame).plane(layer));

% plot the grid points:
figure(1)
plot(xvec,yvec,'o');

forceMesh=createMeshAndBasisFastBEM(xvec,yvec,0);

E=10^4;
v=0.5;
L=10^(-5);
[fx fy x_out y_out M pos_u u] = BEM_force_reconstruction3D(displField(frame).plane,forceMesh,E,v,L,[],[],'fast');

% the forces still have to be scaled with the pixel size to give results in
% Pa.

figure(101)
quiver(x_out,y_out,fx,fy,'b')
hold on
quiver(displField(1).plane(1).pos(:,1),displField(1).plane(1).pos(:,2),displField(1).plane(1).vec(:,1),displField(1).plane(1).vec(:,2),'r')
hold off

% interpolate the force field a long two lines across a high force node:
[maxFmag imax]=max(sqrt(fx(:).^2+fy(:).^2));
xmax=x_out(imax);
ymax=y_out(imax);

xmax=250;
ymax=536;

direct_x=1;
direct_y=0;
width=300;
numPoints=1000;
plotSection(x_out, y_out, fx, fy, xmax, ymax, direct_x, direct_y, width, numPoints);

return;

L=10^(-6);
[fx fy x_out y_out sol_coef]=calcSolFromFwdMapFastBEM(M,u,forceMesh,L,[],[]);

figure(101)
quiver(x_out,y_out,fx,fy,'b')
hold on
quiver(displField(1).plane(1).pos(:,1),displField(1).plane(1).pos(:,2),displField(1).plane(1).vec(:,1),displField(1).plane(1).vec(:,2),'r')
hold off