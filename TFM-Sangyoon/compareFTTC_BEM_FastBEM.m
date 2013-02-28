%% compareFTTC_BEM_FastBEM
% It uses an artificial displacements of beads and reconstruct 
% forces using FTTC and see how correct the method is.

% for displacement field u, use the natural measured mesh.
% for the force field f, use the either the same mesh as for u or construct
% a more dense/sparse mesh by your own.

% create mesh for force field using a set of points and applying a delaunay
% triangulation:

%% parameter setup
E=1000;  %Young's modulus
addNoise=1;
percentNoise=10/100;
doSave=0;
savePath = 'xy0150fastBEM_10pNoise.mat';

s=0.5;  %Poisson's ratio, only needed for FTTC

meshPtsFwdSol=2^10;
L=0;
numPoints_u=50;       %(must be even number)
numPoints_f=50;       %(also must be even number)
numPoints_out=50;

xmin =1;
xmax =50;
ymin =1;
ymax =50;
%% Mesh generation and artificial force generation
[x_mat_u, y_mat_u]=meshgrid(linspace(xmin,xmax,numPoints_u) , linspace(ymin,ymax,numPoints_u));
x_vec_u=reshape(x_mat_u,[],1);
y_vec_u=reshape(y_mat_u,[],1);

% I need force distribution with multiple sources including one at the
% boundary to see the boundary effect from FTTC vs. BEM -Sangyoon 013113
force_x1 = assumedForceShifted(1,x_mat_u,y_mat_u,21,24,1.4,-1.4);
force_y1 = assumedForceShifted(2,x_mat_u,y_mat_u,21,24,1.4,-1.4);
force_x2 = assumedForceShifted(1,x_mat_u,y_mat_u,7,5,2,0.1);
force_y2 = assumedForceShifted(2,x_mat_u,y_mat_u,7,5,2,0.1);
force_x3 = assumedForceShifted(1,x_mat_u,y_mat_u,44,47,0.1,-2);
force_y3 = assumedForceShifted(2,x_mat_u,y_mat_u,44,47,0.1,-2);

force_x = force_x1+force_x2+force_x3;
force_y = force_y1+force_y2+force_y3;
figure(1)
quiver(x_mat_u,y_mat_u,force_x,force_y);
drawnow

%% Forward solution
[ux, uy]=fwdSolution(x_mat_u,y_mat_u,E,xmin,xmax,ymin,ymax,...
    @(x,y) assumedForceShifted(1,x,y,21,24,1.4,-1.4)+...
    assumedForceShifted(1,x,y,7,5,2,0.1)+...
    assumedForceShifted(1,x,y,44,47,0.1,-2),...
    @(x,y) assumedForceShifted(2,x,y,21,24,1.4,-1.4)+...
    assumedForceShifted(2,x,y,7,5,2,0.1)+...
    assumedForceShifted(2,x,y,44,47,0.1,-2),'fft',[],meshPtsFwdSol);
% [ux uy]=fwdSolution(x_mat_u,y_mat_u,E,xmin,xmax,ymin,ymax,@(x,y) assumedForce(2,x,y),@(x,y) assumedForce(2,x,y),'fft',[],meshPtsFwdSol);

if addNoise==1
    max_u=max([max(max(ux)) max(max(uy))]); 
    noise_x=normrnd(0,percentNoise*max_u,numPoints_u,numPoints_u);
    noise_y=normrnd(0,percentNoise*max_u,numPoints_u,numPoints_u);
    ux=ux+noise_x;
    uy=uy+noise_y;
end


u_mat(:,:,1)=ux;
u_mat(:,:,2)=uy;
ux_vec=reshape(ux,[],1);
uy_vec=reshape(uy,[],1);
u=vertcat(ux_vec,uy_vec);

figure(2)
quiver(x_mat_u,y_mat_u,ux,uy);
xlim([xmin-1 xmax+1])
ylim([ymin-1 ymax+1])
drawnow
% %% Forward solution with convolution method - don't use it, it's
% extremely slow!! Sangyoon
% [ux, uy]=fwdSolution(x_mat_u,y_mat_u,E,xmin,xmax,ymin,ymax,@(x,y) assumedForce(2,x,y),@(x,y) assumedForce(2,x,y),'conv',[],meshPtsFwdSol);
% 
% u_mat(:,:,1)=ux;
% u_mat(:,:,2)=uy;
% ux_vec=reshape(ux,[],1);
% uy_vec=reshape(uy,[],1);
% u=vertcat(ux_vec,uy_vec);
% 
% figure(3)
% quiver(x_mat_u,y_mat_u,ux,uy);
%% FTTC-reconstruction with no regularization
pix_durch_my=1; %seem to be important only for the energy
grid_mat_u(:,:,1)=x_mat_u;
grid_mat_u(:,:,2)=y_mat_u;
i_max = size(grid_mat_u,1);
j_max = size(grid_mat_u,2);
cluster_size = grid_mat_u(1,1,1) - grid_mat_u(2,2,1);

% [pos,vec,force_FTTC, fnorm_FTTC] = reg_fourier_TFM_used_till_2010_07_16(grid_mat_u,u_mat,E,s, pix_durch_my, cluster_size, i_max, j_max, 0);  %0.000020
[pos,vec,force_FTTC, fnorm_FTTC] = reg_fourier_TFM(grid_mat_u,u_mat,E,s, pix_durch_my, cluster_size, i_max, j_max, 0);  %0.000020

rec_force_FTTC(:,:,1)=reshape(force_FTTC(:,1),i_max,j_max);
rec_force_FTTC(:,:,2)=reshape(force_FTTC(:,2),i_max,j_max);
%% FTTC-reconstruction with regularization 
pix_durch_my=1; %seem to be important only for the energy
grid_mat_u(:,:,1)=x_mat_u;
grid_mat_u(:,:,2)=y_mat_u;
i_max = size(grid_mat_u,1);
j_max = size(grid_mat_u,2);
cluster_size = grid_mat_u(1,1,1) - grid_mat_u(2,2,1);
reg = 0.00001; %0.0001 made too much smoothing

[pos,vec,force_FTTC, fnorm_FTTC] = reg_fourier_TFM(grid_mat_u,u_mat,E,s, pix_durch_my, cluster_size, i_max, j_max, reg);  

rec_force_FTTCreg(:,:,1)=reshape(force_FTTC(:,1),i_max,j_max);
rec_force_FTTCreg(:,:,2)=reshape(force_FTTC(:,2),i_max,j_max);

%% Now BEM reconstruction
display(['expected computation time:',num2str(numPoints_u^2*numPoints_f^2*27.6*10^(-3)),'s these are:',num2str(numPoints_u^2*numPoints_f^2*27.6*10^(-3)/3600),'h']);
[x_mat_f, y_mat_f]=meshgrid(linspace(xmin,xmax,numPoints_f) , linspace(ymin,ymax,numPoints_f));
x_vec_f=reshape(x_mat_f,[],1);
y_vec_f=reshape(y_mat_f,[],1);
[x_out, y_out]=meshgrid(linspace(xmin,xmax,numPoints_out) , linspace(ymin,ymax,numPoints_out));

tic
forceMesh=createMeshAndBasis(x_vec_f,y_vec_f);
toc
tic
L=0; %regularization factor
[fx_BEM, fy_BEM, x_out, y_out, M, pos_u_M, u_M] = ...
    BEM_force_reconstruction(x_mat_u,y_mat_u,ux,uy,forceMesh,E,L,[],[],[],...
    meshPtsFwdSol);
toc

rec_force_BEM(:,:,1)=fx_BEM;
rec_force_BEM(:,:,2)=fy_BEM;
%% BEM-reconstruction with regularization 
L = 0.00001;
[fx_BEMreg, fy_BEMreg, x_out, y_out, M, pos_u_M, u_Mreg] = BEM_force_reconstruction(x_mat_u,y_mat_u,ux,uy,forceMesh,E,L,[],[],[],meshPtsFwdSol,'fwdMap',M);
rec_force_BEMreg(:,:,1)=fx_BEMreg;
rec_force_BEMreg(:,:,2)=fy_BEMreg;
%% Now FastBEM
tic
forceMeshFastBEM=createMeshAndBasisFastBEM(x_vec_f,y_vec_f,true,[],1);
% forceMesh=createMeshAndBasisFastBEM(xvec,yvec,keepBDPts,[],doPlot);
tic
L = 0;
basisClassTablePath = '/home/sh268/Documents/TFM-simulation/basisFunction5050.mat';
[fx_FastBEM fy_FastBEM x_out y_out M_FastBEM pos_u_M_FastBEM u_M_FastBEM] ...
    = BEM_force_reconstruction(x_mat_u,y_mat_u,ux,uy,forceMeshFastBEM,...
    E,L,x_out,y_out,'fast',meshPtsFwdSol,...
    'QR','basisClassTblPath',basisClassTablePath,'imgRows',49,'imgCols',49);
toc
rec_force_FastBEM(:,:,1)=fx_FastBEM;
rec_force_FastBEM(:,:,2)=fy_FastBEM;
%% FastBEM with regularization
tic
L = 0.00001;
[fx_FastBEMreg fy_FastBEMreg x_outreg y_outreg M_FastBEM pos_u_M_FastBEMreg u_M_FastBEMreg]...
    = BEM_force_reconstruction(x_mat_u,y_mat_u,ux,uy,forceMeshFastBEM,...
    E,L,x_out,y_out,'fast',meshPtsFwdSol,...
    'QR','basisClassTblPath',basisClassTablePath,'imgRows',49,'imgCols',49);
toc
rec_force_FastBEMreg(:,:,1)=fx_FastBEMreg;
rec_force_FastBEMreg(:,:,2)=fy_FastBEMreg;
%% FastBEM-reconstruction with scaled regularization 
L = 0.00001; %same regularization factor
[fx_FastBEMregsc fy_FastBEMregsc x_outregsc y_outregsc M_FastBEMsc pos_u_M_FastBEMregsc u_M_FastBEMregsc]...
    = BEM_force_reconstruction(x_mat_u,y_mat_u,ux,uy,forceMeshFastBEM,...
    E,L,x_out,y_out,'fast',meshPtsFwdSol,...
    'QRscaled','basisClassTblPath',basisClassTablePath,'imgRows',49,'imgCols',49);
toc
rec_force_FastBEMregsc(:,:,1)=fx_FastBEMregsc;
rec_force_FastBEMregsc(:,:,2)=fy_FastBEMregsc;
%% FastBEM-reconstruction with scaled regularization (more regularization case)
L = 1; %more regularization factor
[fx_FastBEMregscmore fy_FastBEMregscmore x_outregsc y_outregsc M_FastBEMscmore pos_u_M_FastBEMregscmore u_M_FastBEMregscmore]...
    = BEM_force_reconstruction(x_mat_u,y_mat_u,ux,uy,forceMeshFastBEM,...
    E,L,x_out,y_out,'fast',meshPtsFwdSol,...
    'QRscaled','basisClassTblPath',basisClassTablePath,'imgRows',49,'imgCols',49);
toc
rec_force_FastBEMregscmore(:,:,1)=fx_FastBEMregscmore;
rec_force_FastBEMregscmore(:,:,2)=fy_FastBEMregscmore;
%% Quiver plots of BEM results compared to original, FTTC, regularized FTTC 
figure(50)
forceScale=3*sqrt(max(max(assumedForce(2,x_mat_u,y_mat_u).^2+assumedForce(2,x_mat_u,y_mat_u).^2)));
quiver(x_out,y_out,rec_force_BEM(:,:,1)/forceScale,rec_force_BEM(:,:,2)/forceScale,0,'b');
hold on
quiver(x_out,y_out,rec_force_BEMreg(:,:,1)/forceScale,rec_force_BEMreg(:,:,2)/forceScale,0,'m');
quiver(grid_mat_u(:,:,1),grid_mat_u(:,:,2),rec_force_FTTC(:,:,1)/forceScale,rec_force_FTTC(:,:,2)/forceScale,0,'r');
quiver(grid_mat_u(:,:,1),grid_mat_u(:,:,2),rec_force_FTTCreg(:,:,1)/forceScale,rec_force_FTTCreg(:,:,2)/forceScale,0,'g');
quiver(x_out,y_out,rec_force_FastBEM(:,:,1)/forceScale,rec_force_FastBEM(:,:,2)/forceScale,0,'c');
quiver(x_out,y_out,rec_force_FastBEMreg(:,:,1)/forceScale,rec_force_FastBEMreg(:,:,2)/forceScale,0,'y');
quiver(x_out,y_out,rec_force_FastBEMregsc(:,:,1)/forceScale,rec_force_FastBEMregsc(:,:,2)/forceScale,0,'g');
quiver(x_out,y_out,rec_force_FastBEMregscmore(:,:,1)/forceScale,rec_force_FastBEMregscmore(:,:,2)/forceScale,0,'r');
quiver(x_mat_u,y_mat_u,force_x/forceScale,force_y/forceScale,0,'k');

hold off
xlim([xmin-1 xmax+1])
ylim([ymin-1 ymax+1])
% Make the text of the legend italic and color it brown
hleg = legend('BEM force','Regularized BEM','FTTC Force','Regularized FTTC',...
    'FastBEM','FastBEMreg','FastBEM scaled L=0.00001',...
    'FastBEM scaled L=1','Org Force',...
              'Location','NorthEastOutside');
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
title('reconstructed forces')
% We see that non-regularized force has 1) incorrect directions on some
% points near high forces 2) forces at the boundary of the image. In case
% of regularized reconstructed forces, the orientation looks better, but
% the it underestimates the high forces compared original given forces.
% BEM is appears to be correct even without the regularization. It has no boundary effect.
%% Heat map presentation for each force magnitude
pos = [x_vec_u y_vec_u];
disp_vec = [ux_vec uy_vec];
rec_force_FTTC_vec = [reshape(rec_force_FTTC(:,:,1),[],1) reshape(rec_force_FTTC(:,:,2),[],1)];
rec_force_FTTCreg_vec = [reshape(rec_force_FTTCreg(:,:,1),[],1) reshape(rec_force_FTTCreg(:,:,2),[],1)];
rec_force_org_vec = [reshape(force_x,[],1) reshape(force_y,[],1)];

[grid_mat,tmat_fttc, i_max, j_max] = interp_vec2grid(pos+disp_vec, rec_force_FTTC_vec,[],grid_mat_u); %1:cluster size
tnorm_fttc = (tmat_fttc(:,:,1).^2 + tmat_fttc(:,:,2).^2).^0.5;
[grid_mat,tmat_fttcreg, i_max, j_max] = interp_vec2grid(pos+disp_vec, rec_force_FTTCreg_vec,[],grid_mat_u); %1:cluster size
tnorm_fttcreg = (tmat_fttcreg(:,:,1).^2 + tmat_fttcreg(:,:,2).^2).^0.5;
force_org(:,:,1)=force_x;
force_org(:,:,2)=force_y;
[grid_mat,tmat_org, i_max, j_max] = interp_vec2grid(pos+disp_vec, rec_force_org_vec,[],grid_mat_u); %1:cluster size
tnorm_org = (tmat_org(:,:,1).^2 + tmat_org(:,:,2).^2).^0.5;

% disp_vec = reshape(u_M,[],2);
rec_force_BEM_vec = [reshape(rec_force_BEM(:,:,1),[],1) reshape(rec_force_BEM(:,:,2),[],1)];
[grid_mat,tmat_BEM, i_max, j_max] = interp_vec2grid(pos+disp_vec, rec_force_BEM_vec,[],grid_mat_u); %1:cluster size
tnorm_BEM = (tmat_BEM(:,:,1).^2 + tmat_BEM(:,:,2).^2).^0.5;

% disp_vec = reshape(u_Mreg,[],2);
rec_force_BEMreg_vec = [reshape(rec_force_BEMreg(:,:,1),[],1) reshape(rec_force_BEMreg(:,:,2),[],1)];
[grid_mat,tmat_BEMreg, i_max, j_max] = interp_vec2grid(pos+disp_vec, rec_force_BEMreg_vec,[],grid_mat_u); %1:cluster size
tnorm_BEMreg = (tmat_BEMreg(:,:,1).^2 + tmat_BEMreg(:,:,2).^2).^0.5;

rec_force_FastBEM_vec = [reshape(rec_force_FastBEM(:,:,1),[],1) reshape(rec_force_FastBEM(:,:,2),[],1)];
[grid_mat,tmat_FastBEM, i_max, j_max] = interp_vec2grid(pos+disp_vec, rec_force_FastBEM_vec,[],grid_mat_u); %1:cluster size
tnorm_FastBEM = (tmat_FastBEM(:,:,1).^2 + tmat_FastBEM(:,:,2).^2).^0.5;

rec_force_FastBEMreg_vec = [reshape(rec_force_FastBEMreg(:,:,1),[],1) reshape(rec_force_FastBEMreg(:,:,2),[],1)];
[grid_mat,tmat_FastBEMreg, i_max, j_max] = interp_vec2grid(pos+disp_vec, rec_force_FastBEMreg,[],grid_mat_u); %1:cluster size
tnorm_FastBEMreg = (tmat_FastBEMreg(:,:,1).^2 + tmat_FastBEMreg(:,:,2).^2).^0.5;

rec_force_FastBEMregscmore_vec = [reshape(rec_force_FastBEMregscmore(:,:,1),[],1) reshape(rec_force_FastBEMregscmore(:,:,2),[],1)];
[grid_mat,tmat_FastBEMregscmore, i_max, j_max] = interp_vec2grid(pos+disp_vec, rec_force_FastBEMregscmore_vec,[],grid_mat_u); %1:cluster size
tnorm_FastBEMregscmore = (tmat_FastBEMregscmore(:,:,1).^2 + tmat_FastBEMregscmore(:,:,2).^2).^0.5;

h51=figure(51); 
set(h51, 'Position', [100 300 900 900])
tmin = min(min(tnorm_org));
tmax = max(max(tnorm_org));
colormap('jet');
subplot(3,3,1), surf(grid_mat(:,:,1), grid_mat(:,:,2), tnorm_org,'FaceColor','interp',...
	'EdgeColor','none', 'FaceLighting','phong'), zlim([tmin-1 tmax+1]), ...
    view(0,90), title('Original force')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
xlim([xmin xmax])
ylim([ymin ymax])
subplot(3,3,2), surf(grid_mat(:,:,1), grid_mat(:,:,2), tnorm_fttc), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('FTTC no regularization')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
xlim([xmin xmax])
ylim([ymin ymax])

subplot(3,3,3), surf(grid_mat(:,:,1), grid_mat(:,:,2), tnorm_fttcreg), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('FTTC regularization')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
xlim([xmin xmax])
ylim([ymin ymax])

subplot(3,3,4), surf(grid_mat(:,:,1), grid_mat(:,:,2), tnorm_BEM), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('BEM no regularization')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
xlim([xmin xmax])
ylim([ymin ymax])

subplot(3,3,5), surf(grid_mat(:,:,1), grid_mat(:,:,2), tnorm_BEMreg), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('BEM regularization')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
xlim([xmin xmax])
ylim([ymin ymax])

subplot(3,3,7), surf(grid_mat(:,:,1), grid_mat(:,:,2), tnorm_FastBEM), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('FastBEM no regularization')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
xlim([xmin xmax])
ylim([ymin ymax])

subplot(3,3,8), surf(grid_mat(:,:,1), grid_mat(:,:,2), tnorm_FastBEMreg), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('FastBEM regularization')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
xlim([xmin xmax])
ylim([ymin ymax])

subplot(3,3,9), surf(grid_mat(:,:,1), grid_mat(:,:,2), tnorm_FastBEMregscmore), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('FastBEM scaled regularization')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
xlim([xmin xmax])
ylim([ymin ymax])

subplot(3,3,6), caxis([tmin-0.01 tmax+0.01]), axis off; colorbar('West')
%% Residuals
% magnitude residual
res_tnorm_fttc = tnorm_org-tnorm_fttc;
res_tnorm_fttcreg = tnorm_org-tnorm_fttcreg;
res_tnorm_BEM = tnorm_org-tnorm_BEM;
res_tnorm_BEMreg = tnorm_org-tnorm_BEMreg;
res_tnorm_FastBEM = tnorm_org-tnorm_FastBEM;
res_tnorm_FastBEMreg = tnorm_org-tnorm_FastBEMreg;
res_tnorm_FastBEMregscmore = tnorm_org-tnorm_FastBEMregscmore;

h52=figure(52); 
set(h52, 'Position', [100 300 900 900])
tmin = min(min(res_tnorm_fttcreg));
tmax = max(max(res_tnorm_fttcreg));
tmin = min(tmin,min(min(res_tnorm_fttc)));
tmax = max(tmax,max(max(res_tnorm_fttc)));
colormap('jet');
subplot(3,3,1), surf(grid_mat(:,:,1), grid_mat(:,:,2), res_tnorm_fttc), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('FTTC residual')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
xlim([xmin xmax])
ylim([ymin ymax])

subplot(3,3,2), surf(grid_mat(:,:,1), grid_mat(:,:,2), res_tnorm_fttcreg), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('FTTC regularization res')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
xlim([xmin xmax])
ylim([ymin ymax])

subplot(3,3,4), surf(grid_mat(:,:,1), grid_mat(:,:,2), res_tnorm_BEM), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('BEM residual')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
xlim([xmin xmax])
ylim([ymin ymax])

subplot(3,3,5), surf(grid_mat(:,:,1), grid_mat(:,:,2), res_tnorm_BEMreg), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('BEM regularization res')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
xlim([xmin xmax])
ylim([ymin ymax])

subplot(3,3,7), surf(grid_mat(:,:,1), grid_mat(:,:,2), res_tnorm_FastBEM), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('FastBEM residual')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
xlim([xmin xmax])
ylim([ymin ymax])

subplot(3,3,8), surf(grid_mat(:,:,1), grid_mat(:,:,2), res_tnorm_FastBEMreg), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('FastBEM regularization res')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
xlim([xmin xmax])
ylim([ymin ymax])

subplot(3,3,9), surf(grid_mat(:,:,1), grid_mat(:,:,2), res_tnorm_FastBEMregscmore), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('FastBEM regularization res')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
xlim([xmin xmax])
ylim([ymin ymax])

subplot(3,3,6), caxis([tmin-0.01 tmax+0.01]), axis off; colorbar('West')
%% Total residual
tot_res(1) = sum((rec_force_org_vec(:,1)-rec_force_FTTC_vec(:,1)).^2+...
    (rec_force_org_vec(:,2)-rec_force_FTTC_vec(:,2)).^2)^.5;
tot_res(2) = sum((rec_force_org_vec(:,1)-rec_force_FTTCreg_vec(:,1)).^2+...
    (rec_force_org_vec(:,2)-rec_force_FTTCreg_vec(:,2)).^2)^.5;
tot_res(3) = sum((rec_force_org_vec(:,1)-rec_force_BEM_vec(:,1)).^2+...
    (rec_force_org_vec(:,2)-rec_force_BEM_vec(:,2)).^2)^.5;
tot_res(4) = sum((rec_force_org_vec(:,1)-rec_force_BEMreg_vec(:,1)).^2+...
    (rec_force_org_vec(:,2)-rec_force_BEMreg_vec(:,2)).^2)^.5;
tot_res(5) = sum((rec_force_org_vec(:,1)-rec_force_FastBEM_vec(:,1)).^2+...
    (rec_force_org_vec(:,2)-rec_force_FastBEM_vec(:,2)).^2)^.5;
tot_res(6) = sum((rec_force_org_vec(:,1)-rec_force_FastBEMreg_vec(:,1)).^2+...
    (rec_force_org_vec(:,2)-rec_force_FastBEMreg_vec(:,2)).^2)^.5;
tot_res(7) = sum((rec_force_org_vec(:,1)-rec_force_FastBEMregscmore_vec(:,1)).^2+...
    (rec_force_org_vec(:,2)-rec_force_FastBEMregscmore_vec(:,2)).^2)^.5;
figure(53); bar(tot_res)
title('mean squared deviation FTTC vs. FTTC regularized vs. BEM vs. BEM regularized vs. FastBEM vs. FastBEM regularized vs. FastBEM scaled regularization');
% first one is fttc without regularization, second bar is total residual
% in regularized fttc. Total residual is lower in regularized fttc.

%% Force at cross-sectional line at y=24
tnorm1D_FastBEMreg = tnorm_FastBEMreg(24,:);
tnorm1D_FastBEMregscmore = tnorm_FastBEMregscmore(24,:);
tnorm1D_org = tnorm_org(24,:);
x1D = grid_mat(24,:,1);
figure(54),plot(x1D,tnorm1D_org,'k')
hold on
plot(x1D,tnorm1D_FastBEMreg,'r')
plot(x1D,tnorm1D_FastBEMregscmore,'g')
xlabel('x')
ylabel('Traction stress at y=24')
%% Peak underestimation
[maxorg_c,tmax_r] = max(tnorm_org);
[maxorg(1),tmax_c] = max(maxorg_c); % maxorg(1) = tnorm_org(tmax_r(tmax_c),tmax_c)
maxorg(2) = tnorm_fttc(tmax_r(tmax_c),tmax_c);
maxorg(3) = tnorm_fttcreg(tmax_r(tmax_c),tmax_c);
maxorg(4) = tnorm_BEM(tmax_r(tmax_c),tmax_c);
maxorg(5) = tnorm_BEMreg(tmax_r(tmax_c),tmax_c);
maxorg(6) = tnorm_FastBEM(tmax_r(tmax_c),tmax_c);
maxorg(7) = tnorm_FastBEMreg(tmax_r(tmax_c),tmax_c);
maxorg(8) = tnorm_FastBEMregscmore(tmax_r(tmax_c),tmax_c);

figure(43); bar(maxorg)
title('Peaks in original, fttc, regularized fttc, BEM, regularized BEM, FastBEM and regularized FastBEM');
% peak is overestimated in FTTC, but underestimated a lot in regularized FTTC.

%% Boundary cutting for FastBEM
tnorm_FastBEM(1,:)=0;
tnorm_FastBEM(:,1)=0;
tnorm_FastBEM(end,:)=0;
tnorm_FastBEM(:,end)=0;

tnorm_FastBEMregscmore(1,:)=0;
tnorm_FastBEMregscmore(:,1)=0;
tnorm_FastBEMregscmore(end,:)=0;
tnorm_FastBEMregscmore(:,end)=0;

%% Precision in locating peak force
peakxy_legend = {'      Original     ','        FTTC       ',' Regularized FTTC  ',...
    '        BEM        ','  Regularized BEM  ','     FastBEM       ','Regularized FastBEM','FastBEM Scaled Reg '};
tnorm(:,:,1) = tnorm_org;
tnorm(:,:,2) = tnorm_fttc;
tnorm(:,:,3) = tnorm_fttcreg;
tnorm(:,:,4) = tnorm_BEM;
tnorm(:,:,5) = tnorm_BEMreg;
tnorm(:,:,6) = tnorm_FastBEM;
tnorm(:,:,7) = tnorm_FastBEMreg;
tnorm(:,:,8) = tnorm_FastBEMregscmore;
zeroI = [-4 -4]; %index of zero coordinates
disp('         Methods          Peak value   X       Y')
for jj=1:8
    clear tmax
%     str = sprintf('%s has a peak of %0.2f at x = %d, y = %d.', peakxy_legend{jj}, maxInt(jj), max_r(jj),max_c(jj));
    [~,locMaxI,sigtVal] = findMaxScoreI(tnorm(:,:,jj),zeroI,3,0.5);
    nPeaks = size(locMaxI,1);
    method = cell(nPeaks,1);
    method(ceil(nPeaks/2)) = peakxy_legend(jj);
    for k=1:nPeaks
        tmax(k,1) = tnorm(locMaxI(k,1),locMaxI(k,2),jj);
    end
    peakinfo = [method num2cell(tmax) num2cell(locMaxI)];
    disp(peakinfo);
end
% We see that BEM as well as FTTC couldn't located the x y position at (10,
% 11), and the values are over-estimated.
%% save
save(savePath);