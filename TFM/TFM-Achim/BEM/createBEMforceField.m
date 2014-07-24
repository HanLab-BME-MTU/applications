displField=strain;

% This is to make the displacement field sparse (too many points for my
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

% First we get a regular grid with the right density of meshpoints. The
% created mesh will be centered in the field of view. Here, only the first
% frame will be used:
[xvec,yvec]=createHexGridFromDisplField(displField);

% plot the grid points:
figure(1)
plot(xvec,yvec,'o');

forceMesh=createMeshAndBasisFastBEM(xvec,yvec,0);

E=10^4;
L=10^(-5);
meshPtsFwdSol=2^11;
[fx fy x_out y_out M pos_u u] = BEM_force_reconstruction(displField(frame).pos(:,1),displField(frame).pos(:,2),displField(frame).vec(:,1),displField(frame).vec(:,2),forceMesh,E,L,[],[],'fast',meshPtsFwdSol);

% the forces still have to be scaled with the pixel size to give results in
% Pa.

figure(100)
quiver(x_out,y_out,fx,fy,'b')
hold on
quiver(displField(1).pos(:,1),displField(1).pos(:,2),displField(1).vec(:,1),displField(1).vec(:,2),'r')
hold off

return;

L=10^(-6);
[fx fy x_out y_out sol_coef]=calcSolFromFwdMapFastBEM(M,u,forceMesh,L,[],[]);

figure(101)
quiver(x_out,y_out,fx,fy,'b')
hold on
quiver(displField(1).pos(:,1),displField(1).pos(:,2),displField(1).vec(:,1),displField(1).vec(:,2),'r')
hold off