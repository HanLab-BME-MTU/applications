% This program (regularized fourier transform traction force
% reconstruction) was produced at the University of Heidelberg, BIOMS
% group of Ulich Schwarz. It calculates traction from a gel displacement
% field.
%
% Benedikt Sabass 20-5-2007
function  [pos,vec,force, fnorm,energie,f] = reg_fourier_TFM(grid_mat,u,E,s, pix_durch_my, cluster_size, i_max, j_max, L)

    nN_pro_pix_fakt = 1/(10^3*pix_durch_my^2);
    nN_pro_my_fakt = 1/(10^3);
    
    V = 2*(1+s)/E;
    
    kx_vec = 2*pi/i_max/cluster_size.*[0:(i_max/2-1) (-i_max/2:-1)];
    ky_vec = 2*pi/j_max/cluster_size.*[0:(j_max/2-1) (-j_max/2:-1)];
    kx = repmat(kx_vec',1,j_max);
    ky = repmat(ky_vec,i_max,1);
    
    kx(1,1) = 1;
    ky(1,1) = 1;
    
    Ginv_xx =(kx.^2+ky.^2).^(-1/2).*V.*(kx.^2.*L+ky.^2.*L+V.^2).^(-1).*(kx.^2.* ...
              L+ky.^2.*L+((-1)+s).^2.*V.^2).^(-1).*(kx.^4.*(L+(-1).*L.*s)+ ...
              kx.^2.*((-1).*ky.^2.*L.*((-2)+s)+(-1).*((-1)+s).*V.^2)+ky.^2.*( ...
              ky.^2.*L+((-1)+s).^2.*V.^2));
    Ginv_yy = (kx.^2+ky.^2).^(-1/2).*V.*(kx.^2.*L+ky.^2.*L+V.^2).^(-1).*(kx.^2.* ...
              L+ky.^2.*L+((-1)+s).^2.*V.^2).^(-1).*(kx.^4.*L+(-1).*ky.^2.*((-1)+ ...
              s).*(ky.^2.*L+V.^2)+kx.^2.*((-1).*ky.^2.*L.*((-2)+s)+((-1)+s).^2.* ...
              V.^2));
    Ginv_xy = (-1).*kx.*ky.*(kx.^2+ky.^2).^(-1/2).*s.*V.*(kx.^2.*L+ky.^2.*L+ ...
              V.^2).^(-1).*(kx.^2.*L+ky.^2.*L+((-1)+s).*V.^2).*(kx.^2.*L+ky.^2.* ...
              L+((-1)+s).^2.*V.^2).^(-1);


    Ginv_xx(1,1) = 0;
    Ginv_yy(1,1) = 0;
    Ginv_xy(1,1) = 0;

    Ginv_xy(i_max/2+1,:) = 0;
    Ginv_xy(:,j_max/2+1) = 0;

    Ftu(:,:,1) = fft2(u(:,:,1));
    Ftu(:,:,2) = fft2(u(:,:,2));

    Ftf(:,:,1) = Ginv_xx.*Ftu(:,:,1) + Ginv_xy.*Ftu(:,:,2);
    Ftf(:,:,2) = Ginv_xy.*Ftu(:,:,1) + Ginv_yy.*Ftu(:,:,2);

    f(:,:,1) = real(ifft2(Ftf(:,:,1)));
    f(:,:,2) = real(ifft2(Ftf(:,:,2)));
    
    pos(:,1) = reshape(grid_mat(:,:,1),i_max*j_max,1);
    pos(:,2) = reshape(grid_mat(:,:,2),i_max*j_max,1);

    vec(:,1) = reshape(u(:,:,1),i_max*j_max,1);
    vec(:,2) = reshape(u(:,:,2),i_max*j_max,1);

    force(:,1) = reshape(f(:,:,1),i_max*j_max,1);
    force(:,2) = reshape(f(:,:,2),i_max*j_max,1);     
   
    fnorm = (force(:,1).^2 + force(:,2).^2).^0.5;
    energie = 1/2*sum(sum(u(2:end-1,2:end-1,1).*f(2:end-1,2:end-1,1) + u(2:end-1,2:end-1,2).*f(2:end-1,2:end-1,2)))*(cluster_size)^2*pix_durch_my^3/10^6; 
end