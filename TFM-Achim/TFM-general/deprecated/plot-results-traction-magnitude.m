for j=1:56
  for i=1:48
    Zmatrix(i,j)=TFM_results(1).traction_magnitude((i-1)*56+j)
  end
end
surf(TFM_results(1).pos(1:56,1),TFM_results(1).pos(1:56:end,2),Zmatrix)


% or use Benedikts functions:
im_direct='/home/ab236/.gvfs/orchestra on files.med.harvard.edu/home/ab236/matlab/data/CONTROL/110906/cell5/c488_im';
make_traction_images(im_direct,TFM_results,TFM_settings); %the input "colorscale" is optional!

figure(1)
% plot the strain field after loading it to the workspace:
quiver(strain(1).pos(:,1),strain(1).pos(:,2),strain(1).vec(:,1),strain(1).vec(:,2))

%interpolate to regular grid using cubic splines (see Benedikts function
%interp_vec2grid):
cluster_size=11;
[grid_mat,u, i_max,j_max] = interp_vec2grid(strain(1).pos, strain(1).vec,cluster_size)

figure(2) 
surf(grid_mat(:,:,1),grid_mat(:,:,2),u(:,:,1));

figure(3) 
surf(grid_mat(:,:,1),grid_mat(:,:,2),u(:,:,2));

u_vec_window=[];
for k=1:5
    u_vec_window=[u_vec_window 
                   u(:,k,1)
                   u(:,k,2)];
end
%das stimmt noch nicht!
s_window_size=5;
[u_noise noise]=wiener2(u(:,1:5,1),[s_window_size,s_window_size]);
%noise=1

s_window_size=3;
%Smooth it with a Wiener Filter:
[u_wiener(:,:,1), noise_1] = wiener2(u(:,:,1),[s_window_size,s_window_size],noise);
[u_wiener(:,:,2), noise_2] = wiener2(u(:,:,2),[s_window_size,s_window_size],noise);


figure(20)
surf(grid_mat(:,:,1),grid_mat(:,:,2),u_wiener(:,:,1));

figure(30)
surf(grid_mat(:,:,1),grid_mat(:,:,2),u_wiener(:,:,2));

% Now calculte the results:
%function  [pos,vec,force, fnorm,energie,f] = reg_fourier_TFM(grid_mat,u,E,s, pix_durch_my, cluster_size, i_max, j_max, L)
E=10;
s=0.5;
pix_durch_my=0.67
i_max=56;
j_max=48;
[pos,vec,force, fnorm,energie,f] = reg_fourier_TFM(grid_mat,u_wiener,E,s, pix_durch_my, cluster_size, i_max, j_max, 0);

figure(40)
quiver(pos(:,1),pos(:,2),force(:,1),force(:,2))


[grid_mat_new,tmat, i_max, j_max] = interp_vec2grid(pos+vec, force, cluster_size); 
tnorm = (tmat(:,:,1).^2 + tmat(:,:,2).^2).^0.5;

figure(41)
surf(grid_mat_new(:,:,1), grid_mat_new(:,:,2), tnorm), shading interp,

figure(42)
colormap default; surf(grid_mat_new(:,:,1), grid_mat_new(:,:,2), tnorm),view(0,90), shading interp, axis equal;
set(gca, 'DataAspectRatio', [1,1,10],'YDir','reverse');
if nargin >3
   set(gca,'CLim',colorscale);
end
colorbar;


[X,Y]=meshgrid( [min(min(grid_mat_new(:,:,1))):1:max(max(grid_mat_new(:,:,1)))] , [min(min(grid_mat_new(:,:,2))):1:max(max(grid_mat_new(:,:,2)))] );
my_X=X';
my_Y=Y';

fnorm_int=interp2(grid_mat_new(:,:,1),grid_mat_new(:,:,2),tnorm,my_X,my_Y);%,'cubic')

figure(43)
imagec(tnorm)