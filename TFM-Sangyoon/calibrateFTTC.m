%% calibrateFTTC 
% It uses an artificial displacements of beads and reconstruct 
% forces using FTTC and see how correct the method is.

% for displacement field u, use the natural measured mesh.
% for the force field f, use the either the same mesh as for u or construct
% a more dense/sparse mesh by your own.

% create mesh for force field using a set of points and applying a delaunay
% triangulation:

%% parameter setup
E=1000;  %Young's modulus
addNoise=0;
percentNoise=10/100;
doSave=0;

s=0.5;  %Poisson's ratio, only needed for FTTC

meshPtsFwdSol=2^10;
L=0;
numPoints_u=100;       %(must be even number)
numPoints_f=100;       %(also must be even number)
numPoints_out=100;

xmin =2;
xmax =98;
ymin =2;
ymax =98;
%% Mesh generation and artificial force generation
[x_mat_u, y_mat_u]=meshgrid(linspace(xmin,xmax,numPoints_u) , linspace(ymin,ymax,numPoints_u));
x_vec_u=reshape(x_mat_u,[],1);
y_vec_u=reshape(y_mat_u,[],1);

% I need force distribution with multiple sources including one at the
% boundary to see the boundary effect from FTTC vs. BEM -Sangyoon 013113
force_x1 = assumedForceShifted(1,x_mat_u,y_mat_u,65,65,7.5,-7.5);
force_y1 = assumedForceShifted(2,x_mat_u,y_mat_u,65,65,7.5,-7.5);
force_x2 = assumedForceShifted(1,x_mat_u,y_mat_u,15,10,10,1);
force_y2 = assumedForceShifted(2,x_mat_u,y_mat_u,15,10,10,1);
force_x3 = assumedForceShifted(1,x_mat_u,y_mat_u,87,90,1,-10);
force_y3 = assumedForceShifted(2,x_mat_u,y_mat_u,87,90,1,-10);

force_x = force_x1+force_x2+force_x3;
force_y = force_y1+force_y2+force_y3;
figure(1)
quiver(x_mat_u,y_mat_u,force_x,force_y);

%% Forward solution
[ux, uy]=fwdSolution(x_mat_u,y_mat_u,E,xmin,xmax,ymin,ymax,...
    @(x,y) assumedForceShifted(1,x,y,65,65,7.5,-7.5)+...
    assumedForceShifted(1,x,y,15,10,10,1)+...
    assumedForceShifted(1,x,y,87,90,1,-10),...
    @(x,y) assumedForceShifted(2,x,y,65,65,7.5,-7.5)+...
    assumedForceShifted(2,x,y,15,10,10,1)+...
    assumedForceShifted(2,x,y,87,90,1,-10),'fft',[],meshPtsFwdSol);
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

%% Evaluation of the force fields 

forceScale=3*sqrt(max(max(assumedForce(2,x_mat_u,y_mat_u).^2+assumedForce(2,x_mat_u,y_mat_u).^2)));

figure(21)
quiver(grid_mat_u(:,:,1),grid_mat_u(:,:,2),rec_force_FTTC(:,:,1)/forceScale,rec_force_FTTC(:,:,2)/forceScale,0,'r');
hold on
quiver(grid_mat_u(:,:,1),grid_mat_u(:,:,2),rec_force_FTTCreg(:,:,1)/forceScale,rec_force_FTTCreg(:,:,2)/forceScale,0,'b');
quiver(x_mat_u,y_mat_u,force_x/forceScale,force_y/forceScale,0,'k');

hold off
xlim([xmin-1 xmax+1])
ylim([ymin-1 ymax+1])
% Make the text of the legend italic and color it brown
hleg = legend('Nonreg Force','Reg Force','Org Force',...
              'Location','NorthEastOutside');
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
title('reconstructed forces')

% We see that non-regularized force has 1) incorrect directions on some
% points near high forces 2) forces at the boundary of the image. In case
% of regularized reconstructed forces, the orientation looks better, but
% the it underestimates the high forces compared original given forces.
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

h22=figure(22); 
set(h22, 'Position', [100 300 1200 300])
tmin = min(min(tnorm_org));
tmax = max(max(tnorm_org));
colormap('jet');
subplot(1,4,1), surf(grid_mat(:,:,1), grid_mat(:,:,2), tnorm_org,'FaceColor','interp',...
	'EdgeColor','none', 'FaceLighting','phong'), zlim([tmin-1 tmax+1]), ...
    view(0,90), title('Original force')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
xlim([xmin xmax])
ylim([ymin ymax])
subplot(1,4,2), surf(grid_mat(:,:,1), grid_mat(:,:,2), tnorm_fttc), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('FTTC no regularization')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
xlim([xmin xmax])
ylim([ymin ymax])

subplot(1,4,3), surf(grid_mat(:,:,1), grid_mat(:,:,2), tnorm_fttcreg), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('FTTC regularization')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
xlim([xmin xmax])
ylim([ymin ymax])

subplot(1,4,4), caxis([tmin-0.01 tmax+0.01]), axis off; colorbar('West')
%% Residuals
% magnitude residual
res_tnorm_fttc = tnorm_org-tnorm_fttc;
res_tnorm_fttcreg = tnorm_org-tnorm_fttcreg;

h41=figure(41); 
set(h41, 'Position', [100 300 900 300])
tmin = min(min(res_tnorm_fttc));
tmax = max(max(res_tnorm_fttc));
colormap('jet');
subplot(1,3,1), surf(grid_mat(:,:,1), grid_mat(:,:,2), res_tnorm_fttc), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('FTTC residual')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
xlim([xmin xmax])
ylim([ymin ymax])

subplot(1,3,2), surf(grid_mat(:,:,1), grid_mat(:,:,2), res_tnorm_fttcreg), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('FTTC regularization res')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
xlim([xmin xmax])
ylim([ymin ymax])

subplot(1,3,3), caxis([tmin-0.01 tmax+0.01]), axis off; colorbar('West')
%% Total residual
tot_res(1) = sum((rec_force_org_vec(:,1)-rec_force_FTTC_vec(:,1)).^2+...
    (rec_force_org_vec(:,2)-rec_force_FTTC_vec(:,2)).^2)^.5;
tot_res(2) = sum((rec_force_org_vec(:,1)-rec_force_FTTCreg_vec(:,1)).^2+...
    (rec_force_org_vec(:,2)-rec_force_FTTCreg_vec(:,2)).^2)^.5;
figure(42); bar(tot_res)
title('mean squared deviation FTTC vs. FTTC regularized');
% first one is fttc without regularization, second bar is total residual
% in regularized fttc. Total residual is lower in regularized fttc.

%% Peak underestimation
[maxorg_c,tmax_r] = max(tnorm_org);
[maxorg(1),tmax_c] = max(maxorg_c); % maxorg(1) = tnorm_org(tmax_r(tmax_c),tmax_c)
maxorg(2) = tnorm_fttc(tmax_r(tmax_c),tmax_c);
maxorg(3) = tnorm_fttcreg(tmax_r(tmax_c),tmax_c);

figure(43); bar(maxorg)
title('Peaks in original, fttc, regularized fttc');
% peak is overestimated in FTTC, but underestimated a lot in regularized FTTC.

%% Precision in locating peak force
% [maxorg_c,tmax_r] = max(tnorm_org);
% [maxorg(1),tmax_c] = max(maxorg_c); % maxorg(1) = tnorm_org(tmax_r(tmax_c),tmax_c)
% max_r(1) = tmax_r(tmax_c); max_c(1) = tmax_c;
% maxInt(1) = maxorg(1);
% [maxfttc_c,tmax_r] = max(tnorm_fttc);
% [maxInt(2),tmax_c] = max(maxfttc_c); % maxorg(1) = tnorm_org(tmax_r(tmax_c),tmax_c)
% max_r(2) = tmax_r(tmax_c); max_c(2) = tmax_c;
% 
% [maxfttcreg_c,tmax_r] = max(tnorm_fttcreg);
% [maxInt(3),tmax_c] = max(maxfttcreg_c); % maxorg(1) = tnorm_org(tmax_r(tmax_c),tmax_c)
% max_r(3) = tmax_r(tmax_c); max_c(3) = tmax_c;
peakxy_legend = {'     Original   ','       FTTC     ','Regularized FTTC'};
tnorm(:,:,1) = tnorm_org;
tnorm(:,:,2) = tnorm_fttc;
tnorm(:,:,3) = tnorm_fttcreg;
zeroI = [-4 -4]; %index of zero coordinates
disp('         Methods          Peak value   X       Y')
for jj=1:3
%     str = sprintf('%s has a peak of %0.2f at x = %d, y = %d.', peakxy_legend{jj}, maxInt(jj), max_r(jj),max_c(jj));
    [~,locMaxI,sigtVal] = findMaxScoreI(tnorm(:,:,jj),zeroI,3);
    nPeaks = size(locMaxI,1);
    method = cell(nPeaks,1);
    method(ceil(nPeaks/2)) = peakxy_legend(jj);
    for k=1:nPeaks
        tmax(k,1) = tnorm(locMaxI(k,1),locMaxI(k,2),jj);
    end
    peakinfo = [method num2cell(tmax) num2cell(locMaxI)];
    disp(peakinfo);
end
% We see that FTTC couldn't located the x y position at (10, 11).

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
[fx_BEM, fy_BEM, x_out, y_out, M, pos_u_M, u_M] = BEM_force_reconstruction(x_mat_u,y_mat_u,ux,uy,forceMesh,E,L,[],[],[],meshPtsFwdSol);
toc

rec_force_BEM(:,:,1)=fx_BEM;
rec_force_BEM(:,:,2)=fy_BEM;
%% BEM-reconstruction with regularization 
L = 0.00001;
[fx_BEMreg, fy_BEMreg, x_out, y_out, M, pos_u_M, u_Mreg] = BEM_force_reconstruction(x_mat_u,y_mat_u,ux,uy,forceMesh,E,L,[],[],[],meshPtsFwdSol);
rec_force_BEMreg(:,:,1)=fx_BEMreg;
rec_force_BEMreg(:,:,2)=fy_BEMreg;
%% BEM results compared to original, FTTC, regularized FTTC 
figure(50)
quiver(x_out,y_out,rec_force_BEM(:,:,1)/forceScale,rec_force_BEM(:,:,2)/forceScale,0,'b');
hold on
quiver(x_out,y_out,rec_force_BEMreg(:,:,1)/forceScale,rec_force_BEMreg(:,:,2)/forceScale,0,'m');
quiver(grid_mat_u(:,:,1),grid_mat_u(:,:,2),rec_force_FTTC(:,:,1)/forceScale,rec_force_FTTC(:,:,2)/forceScale,0,'r');
quiver(grid_mat_u(:,:,1),grid_mat_u(:,:,2),rec_force_FTTCreg(:,:,1)/forceScale,rec_force_FTTCreg(:,:,2)/forceScale,0,'g');
quiver(x_mat_u,y_mat_u,force_x/forceScale,force_y/forceScale,0,'k');

hold off
xlim([xmin-1 xmax+1])
ylim([ymin-1 ymax+1])
% Make the text of the legend italic and color it brown
hleg = legend('BEM force','Regularized BEM','FTTC Force','Regularized FTTC','Org Force',...
              'Location','NorthEastOutside');
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
title('reconstructed forces')
% BEM is appears to be correct even without the regularization. It has no boundary effect.
%% Heat map presentation for each force magnitude
pos = [x_vec_u y_vec_u];
disp_vec = reshape(u_M,[],2);
rec_force_BEM_vec = [reshape(rec_force_BEM(:,:,1),[],1) reshape(rec_force_BEM(:,:,2),[],1)];
[grid_mat,tmat_BEM, i_max, j_max] = interp_vec2grid(pos+disp_vec, rec_force_BEM_vec,[],grid_mat_u); %1:cluster size
tnorm_BEM = (tmat_BEM(:,:,1).^2 + tmat_BEM(:,:,2).^2).^0.5;

% disp_vec = reshape(u_Mreg,[],2);
rec_force_BEMreg_vec = [reshape(rec_force_BEMreg(:,:,1),[],1) reshape(rec_force_BEMreg(:,:,2),[],1)];
[grid_mat,tmat_BEMreg, i_max, j_max] = interp_vec2grid(pos+disp_vec, rec_force_BEMreg_vec,[],grid_mat_u); %1:cluster size
tnorm_BEMreg = (tmat_BEMreg(:,:,1).^2 + tmat_BEMreg(:,:,2).^2).^0.5;

h51=figure(51); 
set(h51, 'Position', [100 300 900 600])
tmin = min(min(tnorm_org));
tmax = max(max(tnorm_org));
colormap('jet');
subplot(2,3,1), surf(grid_mat(:,:,1), grid_mat(:,:,2), tnorm_org,'FaceColor','interp',...
	'EdgeColor','none', 'FaceLighting','phong'), zlim([tmin-1 tmax+1]), ...
    view(0,90), title('Original force')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
xlim([xmin xmax])
ylim([ymin ymax])
subplot(2,3,2), surf(grid_mat(:,:,1), grid_mat(:,:,2), tnorm_fttc), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('FTTC no regularization')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
xlim([xmin xmax])
ylim([ymin ymax])

subplot(2,3,3), surf(grid_mat(:,:,1), grid_mat(:,:,2), tnorm_fttcreg), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('FTTC regularization')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
xlim([xmin xmax])
ylim([ymin ymax])

subplot(2,3,4), surf(grid_mat(:,:,1), grid_mat(:,:,2), tnorm_BEM), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('BEM no regularization')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
xlim([xmin xmax])
ylim([ymin ymax])

subplot(2,3,5), surf(grid_mat(:,:,1), grid_mat(:,:,2), tnorm_BEMreg), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('BEM regularization')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
xlim([xmin xmax])
ylim([ymin ymax])

subplot(2,3,6), caxis([tmin-0.01 tmax+0.01]), axis off; colorbar('West')
%% Residuals
% magnitude residual
res_tnorm_BEM = tnorm_org-tnorm_BEM;
res_tnorm_BEMreg = tnorm_org-tnorm_BEMreg;

h52=figure(52); 
set(h52, 'Position', [100 300 600 600])
tmin = min(min(res_tnorm_fttc));
tmax = max(max(res_tnorm_fttc));
colormap('jet');
subplot(2,3,1), surf(grid_mat(:,:,1), grid_mat(:,:,2), res_tnorm_fttc), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('FTTC residual')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
xlim([xmin xmax])
ylim([ymin ymax])

subplot(2,3,2), surf(grid_mat(:,:,1), grid_mat(:,:,2), res_tnorm_fttcreg), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('FTTC regularization res')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
xlim([xmin xmax])
ylim([ymin ymax])

subplot(2,3,4), surf(grid_mat(:,:,1), grid_mat(:,:,2), res_tnorm_BEM), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('BEM residual')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
xlim([xmin xmax])
ylim([ymin ymax])

subplot(2,3,5), surf(grid_mat(:,:,1), grid_mat(:,:,2), res_tnorm_BEMreg), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('BEM regularization res')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
xlim([xmin xmax])
ylim([ymin ymax])

subplot(2,3,6), caxis([tmin-0.01 tmax+0.01]), axis off; colorbar('West')
%% Total residual
tot_res(3) = sum((rec_force_org_vec(:,1)-rec_force_BEM_vec(:,1)).^2+...
    (rec_force_org_vec(:,2)-rec_force_BEM_vec(:,2)).^2)^.5;
tot_res(4) = sum((rec_force_org_vec(:,1)-rec_force_BEMreg_vec(:,1)).^2+...
    (rec_force_org_vec(:,2)-rec_force_BEMreg_vec(:,2)).^2)^.5;
figure(53); bar(tot_res)
title('mean squared deviation FTTC vs. FTTC regularized vs. BEM vs. BEM regularized');
% first one is fttc without regularization, second bar is total residual
% in regularized fttc. Total residual is lower in regularized fttc.

%% Peak underestimation
maxorg(4) = tnorm_BEM(tmax_r(tmax_c),tmax_c);
maxorg(5) = tnorm_BEMreg(tmax_r(tmax_c),tmax_c);

figure(43); bar(maxorg)
title('Peaks in original, fttc, regularized fttc, BEM and regularized BEM');
% peak is overestimated in FTTC, but underestimated a lot in regularized FTTC.

%% Precision in locating peak force
% [maxorg_c,tmax_r] = max(tnorm_org);
% [maxorg(1),tmax_c] = max(maxorg_c); % maxorg(1) = tnorm_org(tmax_r(tmax_c),tmax_c)
% max_r(1) = tmax_r(tmax_c); max_c(1) = tmax_c;
% maxInt(1) = maxorg(1);
% [maxfttc_c,tmax_r] = max(tnorm_fttc);
% [maxInt(2),tmax_c] = max(maxfttc_c); % maxorg(1) = tnorm_org(tmax_r(tmax_c),tmax_c)
% max_r(2) = tmax_r(tmax_c); max_c(2) = tmax_c;
% 
% [maxfttcreg_c,tmax_r] = max(tnorm_fttcreg);
% [maxInt(3),tmax_c] = max(maxfttcreg_c); % maxorg(1) = tnorm_org(tmax_r(tmax_c),tmax_c)
% max_r(3) = tmax_r(tmax_c); max_c(3) = tmax_c;
peakxy_legend = {'     Original   ','       FTTC     ','Regularized FTTC','       BEM      ',' Regularized BEM'};
tnorm(:,:,1) = tnorm_org;
tnorm(:,:,2) = tnorm_fttc;
tnorm(:,:,3) = tnorm_fttcreg;
tnorm(:,:,4) = tnorm_BEM;
tnorm(:,:,5) = tnorm_BEMreg;
zeroI = [-4 -4]; %index of zero coordinates
disp('         Methods          Peak value   X       Y')
for jj=1:5
    clear tmax
%     str = sprintf('%s has a peak of %0.2f at x = %d, y = %d.', peakxy_legend{jj}, maxInt(jj), max_r(jj),max_c(jj));
    [~,locMaxI,sigtVal] = findMaxScoreI(tnorm(:,:,jj),zeroI,5,0.6);
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