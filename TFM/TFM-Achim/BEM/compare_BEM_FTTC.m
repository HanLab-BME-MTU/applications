% for displacement field u, use the natural measured mesh.
% for the force field f, use the either the same mesh as for u or construct
% a more dense/sparse mesh by your own.

% create mesh for force field using a set of points and applying a delaunay
% triangulation:

E=100;  %Young's modulus
FTTC=1;
%FTTC_old=0;
BEM=0;
FastBEM=0;
addNoise=0;
percentNoise=10/100;
doSave=0;

if FTTC==1
    s=0.5;  %Poisson's ratio, only needed for FTTC
end

meshPtsFwdSol=2^10;
L=0
numPoints_u=20       %(muss gerade Anzahl sein)
numPoints_f=20       %(muss gerade Anzahl sein)
numPoints_out=20

xmin =2;
xmax =12;
ymin =2;
ymax =12;

[x_mat_u y_mat_u]=meshgrid(linspace(xmin,xmax,numPoints_u) , linspace(ymin,ymax,numPoints_u));
x_vec_u=reshape(x_mat_u,[],1);
y_vec_u=reshape(y_mat_u,[],1);

figure(1)
quiver(x_mat_u,y_mat_u,assumedForce(1,x_mat_u,y_mat_u),assumedForce(2,x_mat_u,y_mat_u));

[ux uy]=fwdSolution(x_mat_u,y_mat_u,E,xmin,xmax,ymin,ymax,@(x,y) assumedForce(1,x,y),@(x,y) assumedForce(2,x,y),'fft',[],meshPtsFwdSol);

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

%***************** here starts the FastBEM-reconstruction *********************
if FastBEM==1;
    
    % ['expected computation time:',num2str(numPoints_u^2*numPoints_f^2*27.6*10^(-3)),'s these are:',num2str(numPoints_u^2*numPoints_f^2*27.6*10^(-3)/3600),'h']
    
    [x_mat_f y_mat_f]=meshgrid(linspace(xmin,xmax,numPoints_f) , linspace(ymin,ymax,numPoints_f));
    x_vec_f=reshape(x_mat_f,[],1);
    y_vec_f=reshape(y_mat_f,[],1);

    [x_out y_out]=meshgrid(linspace(xmin,xmax,numPoints_out) , linspace(ymin,ymax,numPoints_out));

    tic
    forceMeshFastBEM=createMeshAndBasisFastBEM(x_vec_f,y_vec_f);
    toc

    %[fx_BEM fy_BEM x_out y_out]=BEM_force_reconstruction(x_mat_u,y_mat_u,ux,uy,forceMesh,E,L);
    %in the upper exmaple, the solution will be calculated on the nodes of the force mesh.
    tic
    [fx_FastBEM fy_FastBEM x_out y_out M_FastBEM pos_u_M_FastBEM u_M_FastBEM] = BEM_force_reconstruction(x_mat_u,y_mat_u,ux,uy,forceMeshFastBEM,E,L,x_out,y_out,'fast',meshPtsFwdSol);
    toc

    rec_force_FastBEM(:,:,1)=fx_FastBEM;
    rec_force_FastBEM(:,:,2)=fy_FastBEM;
end


%***************** here starts the BEM-reconstruction *********************
if BEM==1;
    
    % ['expected computation time:',num2str(numPoints_u^2*numPoints_f^2*27.6*10^(-3)),'s these are:',num2str(numPoints_u^2*numPoints_f^2*27.6*10^(-3)/3600),'h']
    
    [x_mat_f y_mat_f]=meshgrid(linspace(xmin,xmax,numPoints_f) , linspace(ymin,ymax,numPoints_f));
    x_vec_f=reshape(x_mat_f,[],1);
    y_vec_f=reshape(y_mat_f,[],1);

    [x_out y_out]=meshgrid(linspace(xmin,xmax,numPoints_out) , linspace(ymin,ymax,numPoints_out));

    tic
    forceMesh=createMeshAndBasis(x_vec_f,y_vec_f);
    toc

    %[fx_BEM fy_BEM x_out y_out]=BEM_force_reconstruction(x_mat_u,y_mat_u,ux,uy,forceMesh,E,L);
    %in the upper exmaple, the solution will be calculated on the nodes of the force mesh.
    tic
    [fx_BEM fy_BEM x_out y_out M pos_u_M u_M] = BEM_force_reconstruction(x_mat_u,y_mat_u,ux,uy,forceMesh,E,L,x_out,y_out,[],meshPtsFwdSol);
    toc

    rec_force_BEM(:,:,1)=fx_BEM;
    rec_force_BEM(:,:,2)=fy_BEM;
end


%***************** here starts the FTTC-reconstruction ********************
if FTTC==1
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
end

%***************** here starts the FTTC-reconstruction Old %***************
% if FTTC_old==1
%     pix_durch_my=1; %seem to be important only for the energy
%     grid_mat_u(:,:,1)=x_mat_u;
%     grid_mat_u(:,:,2)=y_mat_u;
%     i_max = size(grid_mat_u,1);
%     j_max = size(grid_mat_u,2);
%     cluster_size = grid_mat_u(1,1,1) - grid_mat_u(2,2,1);
%     
%     % [pos,vec,force_FTTC, fnorm_FTTC] = reg_fourier_TFM_used_till_2010_07_16(grid_mat_u,u_mat,E,s, pix_durch_my, cluster_size, i_max, j_max, 0);  %0.000020
%     [pos,vec,force_FTTC_old, fnorm_FTTC_old] = reg_fourier_TFM_used_till_2010_07_16(grid_mat_u,u_mat,E,s, pix_durch_my, cluster_size, i_max, j_max, 0);  %0.000020
%     
%     rec_force_FTTC_old(:,:,1)=reshape(force_FTTC_old(:,1),i_max,j_max);
%     rec_force_FTTC_old(:,:,2)=reshape(force_FTTC_old(:,2),i_max,j_max);
% end

%***************** here starts the evaluation of the force fields *********

forceScale=3*sqrt(max(max(assumedForce(1,x_mat_u,y_mat_u).^2+assumedForce(2,x_mat_u,y_mat_u).^2)));

figure(41)
if FastBEM==1    
    quiver(x_out,y_out,rec_force_FastBEM(:,:,1)/forceScale,rec_force_FastBEM(:,:,2)/forceScale,0,'m');
    hold on
end
if BEM==1
    quiver(x_out,y_out,rec_force_BEM(:,:,1)/forceScale,rec_force_BEM(:,:,2)/forceScale,0,'b');
    hold on
end
if FTTC==1
    quiver(grid_mat_u(:,:,1),grid_mat_u(:,:,2),rec_force_FTTC(:,:,1)/forceScale,rec_force_FTTC(:,:,2)/forceScale,0,'r');
    hold on
end
% if FTTC_old==1
%     quiver(grid_mat_u(:,:,1),grid_mat_u(:,:,2),rec_force_FTTC_old(:,:,1)/forceScale,rec_force_FTTC_old(:,:,2)/forceScale,0,'g');
%     hold on
% end
quiver(x_mat_u,y_mat_u,assumedForce(1,x_mat_u,y_mat_u)/forceScale,assumedForce(2,x_mat_u,y_mat_u)/forceScale,0,'k');
hold off
xlim([xmin xmax])
ylim([ymin ymax])
title('reconstructed forces')

return;


figure(42)
if BEM==1
    res_force_BEM(:,:,1)=assumedForce(1,x_out,y_out)-rec_force_BEM(:,:,1);
    res_force_BEM(:,:,2)=assumedForce(2,x_out,y_out)-rec_force_BEM(:,:,2);
    quiver(x_out,y_out,res_force_BEM(:,:,1),res_force_BEM(:,:,2));
end
if FTTC==1
    hold on
    res_force_FTTC(:,:,1)=assumedForce(1,grid_mat_u(:,:,1),grid_mat_u(:,:,2))-rec_force_FTTC(:,:,1);
    res_force_FTTC(:,:,2)=assumedForce(2,grid_mat_u(:,:,1),grid_mat_u(:,:,2))-rec_force_FTTC(:,:,2);
    quiver(grid_mat_u(:,:,1),grid_mat_u(:,:,2),res_force_FTTC(:,:,1),res_force_FTTC(:,:,2),'r');
    hold off
end
xlim([xmin xmax])
ylim([ymin ymax])
title('residuals')


figure(43)
if BEM==1
    surf(x_out,y_out,sqrt(res_force_BEM(:,:,1).^2+res_force_BEM(:,:,2).^2));
end
if FTTC==1
    hold on
    surf(grid_mat_u(:,:,1),grid_mat_u(:,:,2),sqrt(res_force_FTTC(:,:,1).^2+res_force_FTTC(:,:,2).^2));
    hold off
end 
xlim([xmin xmax])
ylim([ymin ymax])
title('root squared deviation')


figure(44)
if BEM==1
    rel_res_force_BEM(:,:,1)=(assumedForce(1,x_out,y_out)-rec_force_BEM(:,:,1))/max(max(sqrt(assumedForce(1,x_out,y_out).^2+assumedForce(2,x_out,y_out).^2)));
    rel_res_force_BEM(:,:,2)=(assumedForce(2,x_out,y_out)-rec_force_BEM(:,:,2))/max(max(sqrt(assumedForce(1,x_out,y_out).^2+assumedForce(2,x_out,y_out).^2)));
    quiver(x_out,y_out,rel_res_force_BEM(:,:,1),rel_res_force_BEM(:,:,2));
end
if FTTC==1
    hold on
    rel_res_force_FTTC(:,:,1)=(assumedForce(1,grid_mat_u(:,:,1),grid_mat_u(:,:,2))-rec_force_FTTC(:,:,1))/max(max(sqrt(assumedForce(1,grid_mat_u(:,:,1),grid_mat_u(:,:,2)).^2+assumedForce(2,grid_mat_u(:,:,1),grid_mat_u(:,:,2)).^2)));
    rel_res_force_FTTC(:,:,2)=(assumedForce(2,grid_mat_u(:,:,1),grid_mat_u(:,:,2))-rec_force_FTTC(:,:,2))/max(max(sqrt(assumedForce(1,grid_mat_u(:,:,1),grid_mat_u(:,:,2)).^2+assumedForce(2,grid_mat_u(:,:,1),grid_mat_u(:,:,2)).^2)));
    quiver(grid_mat_u(:,:,1),grid_mat_u(:,:,2),rel_res_force_FTTC(:,:,1),rel_res_force_FTTC(:,:,2),'r');
    hold off
end
xlim([xmin xmax])
ylim([ymin ymax])
title('relative residuals')


figure(45)
if BEM==1
    surf(x_out,y_out,sqrt(rel_res_force_BEM(:,:,1).^2+rel_res_force_BEM(:,:,2).^2));
end
if FTTC==1
    hold on
    surf(grid_mat_u(:,:,1),grid_mat_u(:,:,2),sqrt(rel_res_force_FTTC(:,:,1).^2+rel_res_force_FTTC(:,:,2).^2));
    hold off
end
xlim([xmin xmax])
ylim([ymin ymax])
title('root squared relative deviation')

if BEM==1
    'mean squared deviation BEM:'
    sum(sum(res_force_BEM(:,:,1).^2+res_force_BEM(:,:,2).^2))/numPoints_u^2
end

if FTTC==1
    'mean squared deviation FTTC:'
    sum(sum(res_force_FTTC(:,:,1).^2+res_force_FTTC(:,:,2).^2))/numPoints_u^2
end

if BEM==1 && doSave==1
    filename = ['~/matlab/WS-M/WS_BEM-rec_pu_',int2str(numPoints_u),'_pf_',int2str(numPoints_f),'_po_',int2str(numPoints_out),'_noise_',num2str(addNoise),'_L_',num2str(L),'.mat'];
    save(filename)
end

%[x_fine y_fine]=meshgrid(linspace(xmin,xmax,1*numPoints_u) , linspace(ymin,ymax,1*numPoints_u));
%[fx_FTTC_fine fy_FTTC_fine]=calcSolFromFwdMap(M,u_M,forceMesh,L,x_fine,y_fine);
%quiver(x_fine,y_fine,fx_FTTC_fine,fy_FTTC_fine,'g');


%Calculate the total force applied to the surface:
if BEM==1
    [fx_for_quad_BEM fy_for_quad_BEM]=calcSolFromFwdMap(M,u_M,forceMesh,L,forceMesh.p(:,1),forceMesh.p(:,2));

    fx_intp_BEM = TriScatteredInterp(forceMesh.p(:,1),forceMesh.p(:,2),fx_for_quad_BEM,'linear');
    fy_intp_BEM = TriScatteredInterp(forceMesh.p(:,1),forceMesh.p(:,2),fy_for_quad_BEM,'linear');

    %fx_for_quad=@(x,y) calcSolFromFwdMap(M,u_M,forceMesh,L,x,y)
    int_fx_BEM = quad2d(@(x,y) fx_intp_BEM(x,y),xmin,xmax,ymin,ymax)
    int_fx_ass = quad2d(@(x,y) assumedForce(1,x,y),xmin,xmax,ymin,ymax)

    int_fy_BEM = quad2d(@(x,y) fy_intp_BEM(x,y),xmin,xmax,ymin,ymax)
    int_fy_ass = quad2d(@(x,y) assumedForce(2,x,y),xmin,xmax,ymin,ymax)
end

if FTTC==1

    fx_intp_FTTC = TriScatteredInterp(pos(:,1),pos(:,2),force_FTTC(:,1),'linear');
    fy_intp_FTTC = TriScatteredInterp(pos(:,1),pos(:,2),force_FTTC(:,2),'linear');

    %fx_for_quad=@(x,y) calcSolFromFwdMap(M,u_M,forceMesh,L,x,y)
    int_fx_FTTC = quad2d(@(x,y) fx_intp_FTTC(x,y),xmin,xmax,ymin,ymax)
    int_fx_ass = quad2d(@(x,y) assumedForce(1,x,y),xmin,xmax,ymin,ymax)

    int_fy_FTTC = quad2d(@(x,y) fy_intp_FTTC(x,y),xmin,xmax,ymin,ymax)
    int_fy_ass = quad2d(@(x,y) assumedForce(2,x,y),xmin,xmax,ymin,ymax)
end

% plot the L-curve:
% plotLcurve(M,sol_mats,u_M,forceMesh)

return;

%now recalculate the deviations imposing L>0
L_reg=1*10^(-5)

[fx fy]=calcSolFromFwdMap(M,u_M,forceMesh,L_reg);

figure(50)
quiver(x_vec_u,y_vec_u,fx,fy)

figure(420)
res_force_BEM_reg(:,1)=assumedForce(1,x_vec_u,y_vec_u)-fx;
res_force_BEM_reg(:,2)=assumedForce(2,x_vec_u,y_vec_u)-fy;
quiver(x_vec_u,y_vec_u,res_force_BEM_reg(:,1),res_force_BEM_reg(:,2));
xlim([xmin xmax])
ylim([ymin ymax])
title('residuals')

'mean squared deviation BEM regularized:'
sum(res_force_BEM_reg(:,1).^2+res_force_BEM_reg(:,2).^2)/numPoints_u^2

%L_reg=3*10^(-4)
[pos,vec,force_FTTC_reg, fnorm_FTTC_reg] = reg_fourier_TFM(grid_mat_u,u_mat,E,s, pix_durch_my, cluster_size, i_max, j_max, L_reg);  %0.000020
    
rec_force_FTTC_reg(:,:,1)=reshape(force_FTTC_reg(:,1),i_max,j_max);
rec_force_FTTC_reg(:,:,2)=reshape(force_FTTC_reg(:,2),i_max,j_max);

figure(51)
quiver(x_vec_u,y_vec_u,force_FTTC_reg(:,1),force_FTTC_reg(:,2));

figure(421)
res_force_FTTC_reg(:,1)=assumedForce(1,x_vec_u,y_vec_u)-force_FTTC_reg(:,1);
res_force_FTTC_reg(:,2)=assumedForce(2,x_vec_u,y_vec_u)-force_FTTC_reg(:,2);
quiver(x_vec_u,y_vec_u,res_force_FTTC_reg(:,1),res_force_FTTC_reg(:,2));
xlim([xmin xmax]);
ylim([ymin ymax]);
title('residuals');

'mean squared deviation FTTC regularized:'
sum(res_force_FTTC_reg(:,1).^2+res_force_FTTC_reg(:,2).^2)/numPoints_u^2