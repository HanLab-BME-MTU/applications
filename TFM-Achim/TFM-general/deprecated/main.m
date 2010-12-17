xmin =-5;
xmax = 5;
ymin =-5;
ymax = 5;

numPoints=10; %(muss gerade Anzahl sein)

E=100;

[x_mat y_mat]=meshgrid(linspace(xmin,xmax,numPoints) , linspace(xmin,xmax,numPoints));

figure(1)
quiver(x_mat,y_mat,assumedForce(1,x_mat,y_mat),assumedForce(2,x_mat,y_mat));

[ux uy]=fwdSolution(x_mat,y_mat,E,xmin,xmax,ymin,ymax,@(x,y) assumedForce(1,x,y),@(x,y) assumedForce(2,x,y));

figure(2)
quiver(x_mat,y_mat,ux,uy);

u=cat(3,ux,uy);
grid_mat=cat(3,x_mat,y_mat);

s=0.5;
pix_durch_my=1; %seem to be important only for the energy
i_max = size(grid_mat,1);
j_max = size(grid_mat,2);
cluster_size = grid_mat(1,1,1) - grid_mat(2,2,1);


[pos,vec,force, fnorm,energie,f] = reg_fourier_TFM(grid_mat,u,E,s, pix_durch_my, cluster_size, i_max, j_max, 0.000020);  %0.000020


% Use Benedikts version for the creating grid_mat:
x_vec=reshape(x_mat,[],1);
y_vec=reshape(y_mat,[],1);
ux_vec=reshape(ux,[],1);
uy_vec=reshape(uy,[],1);

%[grid_mat_Bene,u_Bene, i_max_Bene, j_max_Bene] = interp_vec2grid([x_vec, y_vec], [ux_vec, uy_vec], cluster_size, [])

%[pos,vec,force, fnorm,energie,f] = reg_fourier_TFM(grid_mat_Bene,u_Bene,E,s, pix_durch_my, cluster_size, i_max_Bene, j_max_Bene, 0);  %0.000020


figure(40)
quiver(pos(:,1),pos(:,2),force(:,1),force(:,2))
title('reconstructed force')

figure(41)
rec_force(:,:,1)=reshape(force(:,1),numPoints,numPoints);
rec_force(:,:,2)=reshape(force(:,2),numPoints,numPoints);
quiver(x_mat,y_mat,rec_force(:,:,1),rec_force(:,:,2));

figure(42)
res_force(:,:,1)=assumedForce(1,x_mat,y_mat)-rec_force(:,:,1);
res_force(:,:,2)=assumedForce(2,x_mat,y_mat)-rec_force(:,:,2);
quiver(x_mat,y_mat,res_force(:,:,1),res_force(:,:,2));

figure(43)
surf(x_mat,y_mat,sqrt(res_force(:,:,1).^2+res_force(:,:,2).^2));
title('root squared deviation')


figure(44)
rel_res_force(:,:,1)=(assumedForce(1,x_mat,y_mat)-rec_force(:,:,1))/max(max(sqrt(assumedForce(1,x_mat,y_mat).^2+assumedForce(2,x_mat,y_mat).^2)));
rel_res_force(:,:,2)=(assumedForce(2,x_mat,y_mat)-rec_force(:,:,2))/max(max(sqrt(assumedForce(1,x_mat,y_mat).^2+assumedForce(2,x_mat,y_mat).^2)));
quiver(x_mat,y_mat,rel_res_force(:,:,1),rel_res_force(:,:,2));

figure(45)
surf(x_mat,y_mat,sqrt(rel_res_force(:,:,1).^2+rel_res_force(:,:,2).^2));
title('root squared relative deviation')

sum(sum(res_force(:,:,1).^2+res_force(:,:,2).^2))/numPoints^2

