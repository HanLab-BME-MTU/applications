%calcForceFromDisplBEM
% calculate the traction force from measured displacement fields using BEM

%readin the displacement field:
E=10000;
L=10^(-4);
load([pwd,'/mech/strain.mat'])

numPoints_u=round(sqrt(length(strain(1).pos(:,1))))
numPoints_f=numPoints_u/2

xmin=min(strain(1).pos(:,1))
xmax=max(strain(1).pos(:,1))
ymin=min(strain(1).pos(:,2))
ymax=max(strain(1).pos(:,2))

[x_mat_f y_mat_f]=meshgrid(linspace(xmin,xmax,numPoints_f) , linspace(ymin,ymax,numPoints_f));
 x_vec_f=reshape(x_mat_f,[],1)
 y_vec_f=reshape(y_mat_f,[],1)


%***************** here starts the BEM-reconstruction *********************
['expected computation time:',num2str(numPoints_u^2*numPoints_f^2*27.6*10^(-3)),'s these are:',num2str(numPoints_u^2*numPoints_f^2*27.6*10^(-3)/3600),'h']    


tic
forceMesh=createMeshAndBasis(x_vec_f,y_vec_f);
toc

%[fx_BEM fy_BEM x_out y_out]=BEM_force_reconstruction(x_mat_u,y_mat_u,ux,uy,forceMesh,E,L);
%in the upper exmaple, the solution will be calculated on the nodes of the force mesh.
tic
[fx_BEM fy_BEM x_out y_out M pos_u u_M] = BEM_force_reconstruction(strain(1).pos(:,1),strain(1).pos(:,2),strain(1).vec(:,1),strain(1).vec(:,2),forceMesh,E,L);
toc

quiver(x_out, y_out, fx_BEM, fy_BEM)