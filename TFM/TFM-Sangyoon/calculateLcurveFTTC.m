function  [rho,eta,reg_corner,alphas] = calculateLcurveFTTC(grid_mat,u,E,s, cluster_size, i_max, j_max, L, LcurveFactor)
% this function calculateLcurveFTTC calculates L-curve by obtaining residual norm and self-norm
% Input : grid_mat, u, cluster_size have to be in the same units, namely
%         pixels. If they are given in different units e.g. meters, then a
%         different scaling factor for the elastic energy has to be
%         applied! The force value remains, however, unaffected, see below.
% Output: The output force is actually a surface stress with the same units
%         as the input E! In particular, the unit of the output force is
%         independent of the units of the input grid_mat,u and cluster_size
%         The reason for this is essentially that the elastic stress is
%         only dependent on the non-dimensional strain which is given by
%         spatial derivatives of the displacements, that is du/dx. If u and
%         dx (essentially cluster_size) are in the same units, then the
%         resulting force has the same dimension as the input E.
%     LcurveFactor = 10;
    alphas=10.^(log10(L)-5:1.25/LcurveFactor:log10(L)+2.5);
    rho=zeros(length(alphas),1);
    eta=zeros(length(alphas),1);
    for i=1:length(alphas);
        disp(['testing L = ' num2str(alphas(i)) '... '])
        curL=alphas(i);
        V = 2*(1+s)/E;

        kx_vec = 2*pi/i_max/cluster_size.*[0:(i_max/2-1) (-i_max/2:-1)];
        ky_vec = 2*pi/j_max/cluster_size.*[0:(j_max/2-1) (-j_max/2:-1)];
        kx = repmat(kx_vec',1,j_max);
        ky = repmat(ky_vec,i_max,1);
        if nargin<10
            Rx=ones(size(kx));
            Ry=ones(size(ky));
        end

        kx(1,1) = 1;
        ky(1,1) = 1;

        X = i_max*cluster_size/2;
        Y = j_max*cluster_size/2; 

        g0x = pi.^(-1).*V.*((-1).*Y.*log((-1).*X+sqrt(X.^2+Y.^2))+Y.*log( ...
          X+sqrt(X.^2+Y.^2))+((-1)+s).*X.*(log((-1).*Y+sqrt(X.^2+Y.^2) ...
          )+(-1).*log(Y+sqrt(X.^2+Y.^2))));
        g0y = pi.^(-1).*V.*(((-1)+s).*Y.*(log((-1).*X+sqrt(X.^2+Y.^2))+( ...
          -1).*log(X+sqrt(X.^2+Y.^2)))+X.*((-1).*log((-1).*Y+sqrt( ...
          X.^2+Y.^2))+log(Y+sqrt(X.^2+Y.^2))));
        % kx, ky: wave vector
        k = (kx.^2+ky.^2).^(-1/2);
        Ginv_xx =k.*V.*(kx.^2.*curL.*Rx+ky.^2.*curL.*Ry+V.^2).^(-1).*(kx.^2.* ...
                  curL.*Rx+ky.^2.*curL.*Ry+((-1)+s).^2.*V.^2).^(-1).*(kx.^4.*(curL.*Rx+(-1).*curL.*s.*Rx)+ ...
                  kx.^2.*((-1).*ky.^2.*curL.*Ry.*((-2)+s)+(-1).*((-1)+s).*V.^2)+ky.^2.*( ...
                  ky.^2.*curL.*Ry+((-1)+s).^2.*V.^2));
        Ginv_yy = k.*V.*(kx.^2.*curL+ky.^2.*curL.*Ry+V.^2).^(-1).*(kx.^2.* ...
                  curL.*Rx+ky.^2.*curL.*Ry+((-1)+s).^2.*V.^2).^(-1).*(kx.^4.*curL+(-1).*ky.^2.*((-1)+ ...
                  s).*(ky.^2.*curL.*Rx+V.^2)+kx.^2.*((-1).*ky.^2.*curL.*Ry.*((-2)+s)+((-1)+s).^2.* ...
                  V.^2));
        Ginv_xy = (-1).*kx.*ky.*k.*s.*V.*(kx.^2.*curL.*Rx+ky.^2.*curL.*Ry+ ...
                  V.^2).^(-1).*(kx.^2.*curL.*Rx+ky.^2.*curL.*Ry+((-1)+s).*V.^2).*(kx.^2.*curL.*Rx+ky.^2.* ...
                  curL.*Ry+((-1)+s).^2.*V.^2).^(-1);


        Ginv_xx(1,1) = 1/g0x;
        Ginv_yy(1,1) = 1/g0y;
        Ginv_xy(1,1) = 0;

        Ginv_xy(i_max/2+1,:) = 0;
        Ginv_xy(:,j_max/2+1) = 0;

        Ftu(:,:,1) = fft2(u(:,:,1));
        Ftu(:,:,2) = fft2(u(:,:,2));

        Ftf(:,:,1) = Ginv_xx.*Ftu(:,:,1) + Ginv_xy.*Ftu(:,:,2);
        Ftf(:,:,2) = Ginv_xy.*Ftu(:,:,1) + Ginv_yy.*Ftu(:,:,2);

        f(:,:,1) = ifft2(Ftf(:,:,1),'symmetric');
        f(:,:,2) = ifft2(Ftf(:,:,2),'symmetric');

        pos(:,1) = reshape(grid_mat(:,:,1),i_max*j_max,1);
        pos(:,2) = reshape(grid_mat(:,:,2),i_max*j_max,1);

        vec(:,1) = reshape(u(:,:,1),i_max*j_max,1);
        vec(:,2) = reshape(u(:,:,2),i_max*j_max,1);

        force(:,1) = reshape(f(:,:,1),i_max*j_max,1);
        force(:,2) = reshape(f(:,:,2),i_max*j_max,1);     

        fnorm = (force(:,1).^2 + force(:,2).^2).^0.5;
        eta(i) = sum(fnorm(:));
        % residual norm calculation
        G_xx=V*((1-s).*k.^2+s*ky.^2)./(k.^3);
        G_xy=V*s*kx.*ky./(k.^3);
        G_yy=V*((1-s).*k.^2+s*kx.^2)./(k.^3);

%         Ginv_xx =k*(-kx.^2+ky.^2.*(-1+s))/(V*(-1+s));
% 
%         Ginv_xx =k.*(kx.^2+(1-s)*ky.^2)./(V*(1-s));
%         G_xxfor =V.*(1-s)./(k.*(k.^2-ky.^2.*s));
        
%         G_xx(1,1) = g0x;
%         G_yy(1,1) = g0y;
        G_xx(1,1) = 0;
        G_yy(1,1) = 0;
        G_xy(1,1) = 0;
        G_xy(i_max/2+1,:) = 0;
        G_xy(:,j_max/2+1) = 0;

        Ftu_estm(:,:,1) = G_xx.*Ftf(:,:,1)+G_xy.*Ftf(:,:,2);
        Ftu_estm(:,:,2) = G_xy.*Ftf(:,:,1)+G_yy.*Ftf(:,:,2);
%         ue(:,:,1) = ifft2(Ftu_estm(:,:,1),'symmetric');
%         ue(:,:,2) = ifft2(Ftu_estm(:,:,2),'symmetric');
        ue(:,:,1) = real(ifft2(Ftu_estm(:,:,1),'symmetric'));
        ue(:,:,2) = real(ifft2(Ftu_estm(:,:,2),'symmetric'));
%         
%         figure, quiver(grid_mat(:,:,1),grid_mat(:,:,2),u(:,:,1),u(:,:,2),0)     
%         hold on
%         quiver(grid_mat(:,:,1),grid_mat(:,:,2),ue(:,:,1),ue(:,:,2),0,'r')
        
%         [ue(:,:,1),ue(:,:,2)]=fwdSolution(grid_mat(:,:,1),grid_mat(:,:,2),E,...
%             grid_mat(1,1,1),grid_mat(end,end,1),grid_mat(1,1,2),grid_mat(end,end,2),...
%             f(:,:,1),f(:,:,2),'fft','Intp',2^10,[],0.5,false,false);
%         [ue(:,:,1),ue(:,:,2)]=fwdSolution(grid_mat(:,:,1),grid_mat(:,:,2),E,...
%             grid_mat(1,1,1),grid_mat(end,end,1),grid_mat(1,1,2),grid_mat(end,end,2),...
%             f(:,:,1),f(:,:,2),'fft','noIntp',[],[],0.5,false,true);
        residueU = ((ue(:,:,1)-u(:,:,1)).^2+(ue(:,:,2)-u(:,:,2)).^2).^0.5;
        rho(i) = sum(residueU(:));
    end
    % Find the L-corner
%     save(LcurveDataPath,'rho','eta','alphas','L','msparse','-v7.3'); % saving before selection.
%     [reg_corner,ireg_corner,~]=regParamSelecetionLcurve(rho,eta,alphas,L);
    [reg_corner,ireg_corner,~]=regParamSelecetionLcurve(rho,eta,alphas,L);
    
end