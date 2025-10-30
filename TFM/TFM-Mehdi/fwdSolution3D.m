function [ux,uy,uz,x_grid,y_grid,z_grid]=fwdSolution3D(x0,y0,z0,E,xmin,xmax,ymin,ymax,zmin,zmax,force_x,force_y,force_z,method,opt,v)
% This forward solution is only valid for a Poisson's ratio v=0.5 if not
% specified.
% Input: No matter what the dimension of x0 and y0 is (pix, or um), the
%        dimension of the surface stresses (force_x, force_y) must have the
%        same dimension as the elastic modulus E, usually Pa.

% Output: The calculated ux and uy have the same dimension as the input
%         x0, y0.
% Achim Besser 2012
% Updated with allowing Poisson's ratio other than 0.5.
% Sangyoon Han, August 2014
% if nargin <18
%     v=0.5;
% %     refine = true;
%     useSameSampling = false;
% elseif nargin <19
% %     refine = true;
%     useSameSampling = false;
% elseif nargin <20
%     useSameSampling = false;
% end
if strcmpi(method,'conv')
    tic;
    disp('Calulate the convolution explicitely in free triangulated mesh')
    [nRow,nCol,nz]=size(x0);

    ux = zeros(nRow,1);
    uy = zeros(nRow,1);
    uz = zeros(nRow,1);
    for i=1:nRow
        for jj=1:nCol
            for kk=1:nz
            integrandx = @(x,y,z) boussinesqGreens3D(1,1,x0(i,jj,kk)-x,y0(i,jj,kk)-y,z0(i,jj,kk)-z,E,v).*force_x(i,jj,kk) + boussinesqGreens3D(1,2,x0(i,jj,kk)-x,y0(i,jj,kk)-y,z0(i,jj,kk)-z,E,v).*force_y(i,jj,kk) + boussinesqGreens3D(1,3,x0(i,jj,kk)-x,y0(i,jj,kk)-y,z0(i,jj,kk)-z,E,v).*force_z(i,jj,kk);
            integrandy = @(x,y,z) boussinesqGreens3D(2,1,x0(i,jj,kk)-x,y0(i,jj,kk)-y,z0(i,jj,kk)-z,E,v).*force_x(i,jj,kk) + boussinesqGreens3D(2,2,x0(i,jj,kk)-x,y0(i,jj,kk)-y,z0(i,jj,kk)-z,E,v).*force_y(i,jj,kk) + boussinesqGreens3D(2,3,x0(i,jj,kk)-x,y0(i,jj,kk)-y,z0(i,jj,kk)-z,E,v).*force_z(i,jj,kk);
            integrandz = @(x,y,z) boussinesqGreens3D(3,1,x0(i,jj,kk)-x,y0(i,jj,kk)-y,z0(i,jj,kk)-z,E,v).*force_x(i,jj,kk) + boussinesqGreens3D(3,2,x0(i,jj,kk)-x,y0(i,jj,kk)-y,z0(i,jj,kk)-z,E,v).*force_y(i,jj,kk) + boussinesqGreens3D(3,3,x0(i,jj,kk)-x,y0(i,jj,kk)-y,z0(i,jj,kk)-z,E,v).*force_z(i,jj,kk);

            ux(i,jj,kk) = integral3(integrandx,xmin,xmax,ymin,ymax,zmin,zmax);
            uy(i,jj,kk) = integral3(integrandy,xmin,xmax,ymin,ymax,zmin,zmax);
            uz(i,jj,kk) = integral3(integrandz,xmin,xmax,ymin,ymax,zmin,zmax);% RelTol sucks! 'RelTol',5e-13);
            end
         end
    end

 

elseif strcmpi(method,'fft')
    disp('Use fast convolution')    

    %***************************************************************
    % Here starts the calculation using the fast fourier transform *
    %***************************************************************
    
    % Number of points to calculate the force field and the Greensfunction.
    % Since below Nx_G and Ny_G are odd, the Greensfunctions will be
    % evaluated at zero. Since the Greensfunctions at x=0 diverges, this
    % will in general cause a problem. For this I have found a work around.
    % The Greensfunction will be evaluated as is and the divergent value
    % at x=0 will be set to zero, this part of the support will be
    % integrated seperately. In order to do this, I assume that for dense
    % sampling, the force field doesn't very strongly over one gridsize.
    % Assuming it to be constant around x=0, allows to integrate the
    % Greensfunction around a domain x=+-r and y=+-r. This yields a
    % correction term which is of particular importance for sparse
    % sampling, meaning that Nx_F is small. This alogorithm performs very
    % well and has been cross-validated with the results obtained using the
    % 'conv' option. If you want to repeat the test, use the few lines of
    % code at the very end of this function.
    
    %tic;
    
    % This determines the sampling of the force field:
%     if (nargin < 16 || isempty(meshPtsFwdSol)) && ~useSameSampling
%         disp('Use meshPtsFwdSol=2^10. This value should be given with the function call!!!');
%         meshPtsFwdSol=2^10;
%     end
    
    
        Nx_F=2^7; % 2^10 is the densest sampling possible.
        Ny_F=Nx_F;
        Nz_F=Nx_F;
    
    
    % To account for dx*dy in the convolution integral one has to finally
    % rescale the result by the following scaling factor:
     xRange=(max(max(max(x0)))-min(min(min(x0))));
     yRange=(max(max(max(y0)))-min(min(min(y0))));
     zRange=(max(max(max(z0)))-min(min(min(z0))));


    scalingFactor=(xRange*yRange*zRange)/(Nx_F*Ny_F*Nz_F);
%     scalingFactor=1;
    
    % To cover the whole support of the force field, the domain over which
    % the Greensfunctions have to be calculated need to be at least of size:
    % (2Nx-1)x(2Ny-1).

    Nx_G=2*Nx_F-1;
    Ny_G=2*Ny_F-1;
    Nz_G=2*Nz_F-1;
    
    % Subsequently, these have to be padded with zeros according to:
    Nx_pad=Nx_F+Nx_G-1;
    Ny_pad=Ny_F+Ny_G-1;
    Nz_pad=Nz_F+Nz_G-1;
    
    % These might not be a power of 2, make sure that they are:
    Nx_pad=pow2(nextpow2(Nx_pad));
    Ny_pad=pow2(nextpow2(Ny_pad));
    Nz_pad=pow2(nextpow2(Nz_pad));


    % First determine the boundaries of the mesh:
    leftUpperCorner =[min(min(min(x0))) min(min(min(y0))) min(min(min(z0)))];
    rightLowerCorner=[max(max(max(x0))) max(max(max(y0))) max(max(max(z0)))];

    % create a regular mesh with Nx*Ny meshpoints where the force field is
    % calculated. This need not to be a power of 2 yet:
    xvec_F=linspace(leftUpperCorner(1),rightLowerCorner(1),Nx_F);
    yvec_F=linspace(leftUpperCorner(2),rightLowerCorner(2),Ny_F);
    zvec_F=linspace(leftUpperCorner(3),rightLowerCorner(3),Nz_F);

    [xgrid_F,ygrid_F,zgrid_F]=meshgrid(xvec_F,yvec_F,zvec_F);
    % create a mesh centered at zero with Nx_G*Ny_G meshpoints, where the
    % Greensfunctions are calculated.
    xvec_G=linspace(-xRange,xRange,Nx_G);
    yvec_G=linspace(-yRange,yRange,Ny_G);
    zvec_G=linspace(-zRange,zRange,Nz_G);

%     [xgrid_G,ygrid_G]=meshgrid(yvec_G,xvec_G);
    [xgrid_G,ygrid_G,zgrid_G]=meshgrid(xvec_G,yvec_G,zvec_G);
      
    %calculate the force values at the grid_F positions:
%       if useSameSampling
%           discrete_Force_x=force_x(xgrid_F,ygrid_F); %this has only to be calculated over the support xmin,xmax,ymin,ymax rest is zero
%           discrete_Force_y=force_y(xgrid_F,ygrid_F);
%           discrete_Force_z=force_z(xgrid_F,ygrid_F);
%       else
    % Make force_x and force_y a function handle if it is a matrix
%           if ismatrix(force_x) && ismatrix(force_y) && ismatrix(force_y) &&
          if ~isa(force_x,'TriScatteredInterp') && ~isa(force_x,'function_handle')
        %     [xmat,ymat]=meshgrid(xmin:xmax,ymin:ymax);
        %     xvec=xmat(:);
        %     yvec=ymat(:);
            xvec=x0(:);
            yvec=y0(:);
            zvec=z0(:);
            force_x_vec=force_x(:);
            force_y_vec=force_y(:);
            force_z_vec=force_z(:);
            force_x = scatteredInterpolant(xvec,yvec,zvec,force_x_vec);
            force_y = scatteredInterpolant(xvec,yvec,zvec,force_y_vec);
            force_z = scatteredInterpolant(xvec,yvec,zvec,force_z_vec);

         end
        discrete_Force_x=force_x(xgrid_F,ygrid_F,zgrid_F); %this has only to be calculated over the support xmin,xmax,ymin,ymax rest is zero
        discrete_Force_y=force_y(xgrid_F,ygrid_F,zgrid_F);
        discrete_Force_z=force_z(xgrid_F,ygrid_F,zgrid_F);

%         taking care of nans
         checkVec=isnan(discrete_Force_x);
         discrete_Force_x(checkVec)=0;
         checkVec=isnan(discrete_Force_y);
         discrete_Force_y(checkVec)=0;
         checkVec=isnan(discrete_Force_z);
         discrete_Force_z(checkVec)=0;
         
%       end
    
    % Calculate the Greens-function values at the grid_G positions. This can
    % be improved since the Greensfunction never change for a given grid
    % size. When the Basis functions are calculated this has to be done
    % only once (as well as the FFT for these fields!!!):
    discrete_boussinesqGreens3D11=boussinesqGreens3D(1,1,xgrid_G,ygrid_G,zgrid_G);
    discrete_boussinesqGreens3D12=boussinesqGreens3D(1,2,xgrid_G,ygrid_G,zgrid_G);
    discrete_boussinesqGreens3D13=boussinesqGreens3D(1,3,xgrid_G,ygrid_G,zgrid_G);
    discrete_boussinesqGreens3D21=boussinesqGreens3D(2,1,xgrid_G,ygrid_G,zgrid_G); %#ok
    discrete_boussinesqGreens3D22=boussinesqGreens3D(2,2,xgrid_G,ygrid_G,zgrid_G);
    discrete_boussinesqGreens3D23=boussinesqGreens3D(2,3,xgrid_G,ygrid_G,zgrid_G);
    discrete_boussinesqGreens3D31=boussinesqGreens3D(3,1,xgrid_G,ygrid_G,zgrid_G);
    discrete_boussinesqGreens3D32=boussinesqGreens3D(3,2,xgrid_G,ygrid_G,zgrid_G);
    discrete_boussinesqGreens3D33=boussinesqGreens3D(3,3,xgrid_G,ygrid_G,zgrid_G);
    
    % Pad the calculated fields with zero to the next power larger than 
    % (2*N-1), see above. For this setup, the FFT is fastest.
%     discrete_Force_x=padarray(discrete_Force_x_unPadded,[Nx_pad-Nx_F Ny_pad-Ny_F],0,'post');%'symmetric','post');
% %     discrete_Force_y=padarray(discrete_Force_y_unPadded,[Nx_pad-Nx_F Ny_pad-Ny_F],0,'post');
    discrete_Force_x=padarray(discrete_Force_x,[Nz_pad-Nz_F Ny_pad-Ny_F Nx_pad-Nx_F],0,'post');%'symmetric','post');
     discrete_Force_y=padarray(discrete_Force_y,[Nz_pad-Nz_F Ny_pad-Ny_F Nx_pad-Nx_F],0,'post');
     discrete_Force_z=padarray(discrete_Force_z,[Nz_pad-Nz_F Ny_pad-Ny_F Nx_pad-Nx_F],0,'post');
   

    
%     discrete_boussinesqGreens11=padarray(discrete_boussinesqGreens11,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
%     discrete_boussinesqGreens12=padarray(discrete_boussinesqGreens12,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
%    %discrete_boussinesqGreens22=padarray(discrete_boussinesqGreens22,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post'); 
%     discrete_boussinesqGreens22=padarray(discrete_boussinesqGreens22,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
    
     discrete_boussinesqGreens3D11=padarray(discrete_boussinesqGreens3D11,[Nz_pad-Nz_G Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
     discrete_boussinesqGreens3D12=padarray(discrete_boussinesqGreens3D12,[Nz_pad-Nz_G Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
     discrete_boussinesqGreens3D13=padarray(discrete_boussinesqGreens3D13,[Nz_pad-Nz_G Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
%      discrete_boussinesqGreens3D21=padarray(discrete_boussinesqGreens3D21,[Nx_pad-Nx_G Ny_pad-Ny_G Nz_pad-Nz_G],0,'post');
    discrete_boussinesqGreens3D22=padarray(discrete_boussinesqGreens3D22,[Nz_pad-Nz_G Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
     discrete_boussinesqGreens3D23=padarray(discrete_boussinesqGreens3D23,[Nz_pad-Nz_G Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
     discrete_boussinesqGreens3D31=padarray(discrete_boussinesqGreens3D31,[Nz_pad-Nz_G Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
     discrete_boussinesqGreens3D32=padarray(discrete_boussinesqGreens3D32,[Nz_pad-Nz_G Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
     discrete_boussinesqGreens3D33=padarray(discrete_boussinesqGreens3D33,[Nz_pad-Nz_G Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
    
    % Now calculate the fourier transforms:
    dFT_Force_x=fftn(discrete_Force_x);
    clear discrete_Force_x;
    dFT_Force_y=fftn(discrete_Force_y);
    clear discrete_Force_y;
    dFT_Force_z=fftn(discrete_Force_z);
    clear discrete_Force_z;
    
    % This has to be calculated only once for all basis functions!
%     dFT_boussinesqGreens11=fft2(discrete_boussinesqGreens11);
%     clear discrete_boussinesqGreens11;
%     dFT_boussinesqGreens12=fft2(discrete_boussinesqGreens12);
%     clear discrete_boussinesqGreens12;
%     dFT_boussinesqGreens21=dFT_boussinesqGreens12;
%     % nothing to clear here!
%     dFT_boussinesqGreens22=fft2(discrete_boussinesqGreens22);
%     clear discrete_boussinesqGreens22;
    dFT_boussinesqGreens3D11=fftn(discrete_boussinesqGreens3D11);
    dFT_boussinesqGreens3D12=fftn(discrete_boussinesqGreens3D12);
    dFT_boussinesqGreens3D13=fftn(discrete_boussinesqGreens3D13);
    dFT_boussinesqGreens3D21=dFT_boussinesqGreens3D12;
    dFT_boussinesqGreens3D22=fftn(discrete_boussinesqGreens3D22);
    dFT_boussinesqGreens3D23=fftn(discrete_boussinesqGreens3D23);
    dFT_boussinesqGreens3D31=fftn(discrete_boussinesqGreens3D31);
    dFT_boussinesqGreens3D32=fftn(discrete_boussinesqGreens3D32);
    dFT_boussinesqGreens3D33=fftn(discrete_boussinesqGreens3D33);
    
    % Now calculate the solution:                
    ux_grid=ifftn(dFT_boussinesqGreens3D11.*dFT_Force_x+dFT_boussinesqGreens3D12.*dFT_Force_y+dFT_boussinesqGreens3D13.*dFT_Force_z);
    uy_grid=ifftn(dFT_boussinesqGreens3D21.*dFT_Force_x+dFT_boussinesqGreens3D22.*dFT_Force_y+dFT_boussinesqGreens3D23.*dFT_Force_z);
    uz_grid=ifftn(dFT_boussinesqGreens3D31.*dFT_Force_x+dFT_boussinesqGreens3D32.*dFT_Force_y+dFT_boussinesqGreens3D33.*dFT_Force_z); 
    
    % Plot the solution:
%     figure(10)
%     imshow(ux_grid,[])
%     
%     figure(11)
%     surf(uy_grid)
    
    % Now extract the essential part from the solution. It is located in
    % the center of the padded field.    
    % I really don't understand why to cut it out like this, but it works!
    startIndex_x=abs(Nx_G-Nx_F)+1; % Or is it just: stzartIndex_x=Nx_F
    startIndex_y=abs(Ny_G-Ny_F)+1;
    startIndex_z=abs(Nz_G-Nz_F)+1;

    endIndex_x=startIndex_x+Nx_F-1;
    endIndex_y=startIndex_y+Ny_F-1;
    endIndex_z=startIndex_z+Nz_F-1;

    

    % Remove imaginary part caused by round off errors:
%      ux_grid=real(ux_grid(1:128,1:128,1:128));
%      uy_grid=real(uy_grid(1:128,1:128,1:128));
%      uy_grid=real(uy_grid(32:95,32:95,32:95));
%      uy_grid=real(uy_grid(startIndex_x:endIndex_x,startIndex_y:endIndex_y));
     ux_grid=real(ux_grid(startIndex_y:endIndex_y,startIndex_x:endIndex_x,startIndex_z:endIndex_z));
     uy_grid=real(uy_grid(startIndex_y:endIndex_y,startIndex_x:endIndex_x,startIndex_z:endIndex_z));
     uz_grid=real(uz_grid(startIndex_y:endIndex_y,startIndex_x:endIndex_x,startIndex_z:endIndex_z));


%!!! This could be improved by using the analytical solution for the Fourie
%!!! Transform of the Greensfunction!
    % Add the solution for G(0,0). This is a correction term which becomes
    % irrelevant for very dense sampling. But for small Nx_F it is REALLY
    % essential!
    % Set the Poisson's ratio to 0.5:
    v=0.5;
%     dx=abs(xvec_G(2)-xvec_G(1));
%     dy=abs(yvec_G(2)-yvec_G(1));
%     dz=abs(zvec_G(2)-zvec_G(1));
    
%     int_x2_over_r3=2*dy*log((dy^2+2*dx*(dx+sqrt(dx^2+dy^2)))/(dy^2));    
%     int_y2_over_r3=2*dx*log((dx^2+2*dy*(dy+sqrt(dx^2+dy^2)))/(dx^2));    
%     int_1_over_r  =int_x2_over_r3 + int_y2_over_r3;
        
%     corrTerm_11=(1+v)/(pi*E)*((1-v)*int_1_over_r+v*int_x2_over_r3);
%     corrTerm_22=(1+v)/(pi*E)*((1-v)*int_1_over_r+v*int_y2_over_r3);
    r=1;
    corrTermXY =(1+v)*(2-v)/(2*pi*E)*(4*r*log(3+2*sqrt(2)));
    corrTermZ  =2*(1-v)/(2-v)*corrTermXY;
    
    ux_grid=ux_grid+force_x(xgrid_F,ygrid_F,zgrid_F)*corrTermXY;
    uy_grid=uy_grid+force_y(xgrid_F,ygrid_F,zgrid_F)*corrTermXY;
    uz_grid=uz_grid+force_z(xgrid_F,ygrid_F,zgrid_F)*corrTermZ;    
    
   
    
    
    
    % Recursive call to fwdSolution:   
%     refine =true;
%     if ~isempty(xmin) && ~isempty(xmax) && ~isempty(ymin) && ~isempty(ymax) && ~isempty(zmin) && ~isempty(zmax) && refine
%         % and the support of the force is much small than ROI for the
%         % displacement, this check should be included otherwise one calculate
%         % the same thing twice
%         % extend the size of the region a little bit, here by a factor of 1.2,
%         % but this is arbitrary.
%         % The range of the support is:
%         xsupp=xmax-xmin;
%         ysupp=ymax-ymin;
%         zsupp=zmax-zmin;
%         
%         % Now expand the x- and y-range:
%         xminExp=xmin-xsupp/10;
%         xmaxExp=xmax+xsupp/10;
%         yminExp=ymin-ysupp/10;
%         ymaxExp=ymax+ysupp/10;
%         zminExp=zmin-zsupp/10;
%         zmaxExp=zmax+zsupp/10;
%         
%         [ux_fine uy_fine uz_fine x_grid_fine y_grid_fine z_grid_fine]=fwdSolution([xminExp xmaxExp],[yminExp ymaxExp],[zminExp zmaxExp],E,[],[],[],[],force_x,force_y,force_z,method,'noIntp',meshPtsFwdSol);
%         
%         % Later on we want to have ux and uy defined on a regular
%         % grid. for this reason we now interpolate the fine solution onto
%         % the regular xgrid_F and ygrid_F. This grid is usually so fine
%         % that it already corresponds to subpixel sampling: (e.g. 2^10=1024
%         % positions along each image dimension)
%         
%         % find the positions that are within
%         
%         iux_fine = griddata(x_grid_fine,y_grid_fine,z_grid_fine,ux_fine,xgrid_F,ygrid_F,zgrid_F,'linear');%'*cubic'
%         iuy_fine = griddata(x_grid_fine,y_grid_fine,uy_fine,xgrid_F,ygrid_F,'linear');%'*linear'
%         iuz_fine = griddata(x_grid_fine,y_grid_fine,z_grid_fine,ux_fine,xgrid_F,ygrid_F,zgrid_F,'linear');%'*cubic'
% 
%         % those contain a lot of NaNs since most of the points are outside
%         % of [xminExp xmaxExp] and [yminExp ymaxExp]. This will give us two
%         % masks that we can now use to update ux_grid and uy_grid:        
%         valMat=~isnan(iux_fine);
%         
%         diff=2*abs(uy_grid(valMat)-iuy_fine(valMat))./abs(uy_grid(valMat)+iuy_fine(valMat));
%         display(['max. correction: ',num2str(max(diff(:)))]);
%         
%         % Update the values in the coarse solution:
%         ux_grid(valMat)=iux_fine(valMat);
%         uy_grid(valMat)=iuy_fine(valMat);
%     end
    
        
    % Now interpolate the displacement field from the regular grid to the irregular
    % measured grid. Since we have the force field defined on a regular grid
    % we can use the fast *option. 'linear' is about a factor of two faster than
    % 'cubic'. Hard to tell if cubic performs better than linear.    
     if nargin==18 && strcmp(opt,'noIntp')
        ux =scalingFactor*(ux_grid);%'*cubic'
        uy =scalingFactor*(uy_grid);%'*linear'
        uz =scalingFactor*(uz_grid);
        x_grid=xgrid_F;
        y_grid=ygrid_F;
        z_grid=zgrid_F;
     else
        ux =scalingFactor*interp3(xgrid_F,ygrid_F,zgrid_F,ux_grid,x0,y0,z0);%'*cubic'
        uy =scalingFactor*interp3(xgrid_F,ygrid_F,zgrid_F,uy_grid,x0,y0,z0);%'*linear'
        uz =scalingFactor*interp3(xgrid_F,ygrid_F,zgrid_F,uz_grid,x0,y0,z0);
        x_grid=x0;
        y_grid=y0;
        z_grid=z0;
     end
    
end
    
    %toc;
    
    %x0_vec=reshape(x0,[],1);
    %y0_vec=reshape(y0,[],1);

    %x_vec=reshape(xgrid_F,[],1);
    %y_vec=reshape(ygrid_F,[],1);
    %ux_vec=reshape(ux_grid,[],1);
    %uy_vec=reshape(uy_grid,[],1);


    %tic;
    %[~, ~, ux] = griddata(x_vec,y_vec,ux_vec,x0,y0,'cubic');
    %[~, ~, uy] = griddata(x_vec,y_vec,uy_vec,x0,y0,'cubic');
    %toc;
% elseif strcmpi(method,'fftn_finite')
%     %display('Use fast convolution')    
% 
%     % This determines the sampling of the force field:
%     if nargin < 16 || isempty(meshPtsFwdSol)
%         display('Use meshPtsFwdSol=2^10. This value should be given with the function call!!!');
%         meshPtsFwdSol=2^10;
%     end
%         
%     Nx_F=meshPtsFwdSol; % 2^10 is the densest sampling possible.
%     Ny_F=Nx_F;
%     
%     % To account for dx*dy in the convolution integral one has to finally
%     % rescale the result by the following scaling factor:
%     xRange=(max(max(x0))-min(min(x0)));
%     yRange=(max(max(y0))-min(min(y0)));
%     scalingFactor=(xRange*yRange)/(Nx_F*Ny_F);
%     
%     % To cover the whole support of the force field, the domain over which
%     % the Greensfunctions have to be calculated need to be at least of size:
%     % (2Nx-1)x(2Ny-1).
% 
%     Nx_G=2*Nx_F-1;
%     Ny_G=2*Ny_F-1;    
%     
%     % Subsequently, these have to be padded with zeros according to:
%     Nx_pad=Nx_F+Nx_G-1;
%     Ny_pad=Ny_F+Ny_G-1;
%     
%     % These might not be a power of 2, make sure that they are:
%     Nx_pad=pow2(nextpow2(Nx_pad));
%     Ny_pad=pow2(nextpow2(Ny_pad));
% 
%     % First determine the boundaries of the mesh:
%     leftUpperCorner =[min(min(x0)) min(min(y0))];
%     rightLowerCorner=[max(max(x0)) max(max(y0))];
% 
%     % create a regular mesh with Nx*Ny meshpoints where the force field is
%     % calculated. This need not to be a power of 2 yet:
%     xvec_F=linspace(leftUpperCorner(1),rightLowerCorner(1),Nx_F);
%     yvec_F=linspace(leftUpperCorner(2),rightLowerCorner(2),Ny_F);
%     [xgrid_F,ygrid_F]=meshgrid(xvec_F,yvec_F);
%     
%     % create a mesh centered at zero with Nx_G*Ny_G meshpoints, where the
%     % Greensfunctions are calculated.
%     xvec_G=linspace(-xRange,xRange,Nx_G);
%     yvec_G=linspace(-yRange,yRange,Ny_G);
%     
%     [xgrid_G,ygrid_G]=meshgrid(xvec_G,yvec_G);
%       
%     %calculate the force values at the grid_F positions:
%     discrete_Force_x_unPadded=force_x(xgrid_F,ygrid_F); %this has only to be calculated over the support xmin,xmax,ymin,ymax rest is zero
%     discrete_Force_y_unPadded=force_y(xgrid_F,ygrid_F);
% 
%     % This part if what's different from 'fft' method that uses
%     % boussinesque greens function
%     discrete_boussinesqGreens11=finiteThicknessGreens(1,1,xgrid_G,ygrid_G,E,h);
%     discrete_boussinesqGreens12=finiteThicknessGreens(1,2,xgrid_G,ygrid_G,E,h);
%    %discrete_boussinesqGreens21=discrete_boussinesqGreens12;
%     discrete_boussinesqGreens22=finiteThicknessGreens(2,2,xgrid_G,ygrid_G,E,h);
%     
%     % Pad the calculated fields with zero to the next power larger than 
%     % (2*N-1), see above. For this setup, the FFT is fastest.
%     discrete_Force_x=padarray(discrete_Force_x_unPadded,[Nx_pad-Nx_F Ny_pad-Ny_F],0,'post');%'symmetric','post');
%     discrete_Force_y=padarray(discrete_Force_y_unPadded,[Nx_pad-Nx_F Ny_pad-Ny_F],0,'post');
%     
%     discrete_boussinesqGreens11=padarray(discrete_boussinesqGreens11,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
%     discrete_boussinesqGreens12=padarray(discrete_boussinesqGreens12,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
%    %discrete_boussinesqGreens22=padarray(discrete_boussinesqGreens22,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post'); 
%     discrete_boussinesqGreens22=padarray(discrete_boussinesqGreens22,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
%     
%     % Now calculate the fourier transforms:
%     dFT_Force_x=fft2(discrete_Force_x);
%     clear discrete_Force_x;
%     dFT_Force_y=fft2(discrete_Force_y);
%     clear discrete_Force_y;
%     
%     % This has to be calculated only once for all basis functions!
%     dFT_boussinesqGreens11=fft2(discrete_boussinesqGreens11);
%     clear discrete_boussinesqGreens11;
%     dFT_boussinesqGreens12=fft2(discrete_boussinesqGreens12);
%     clear discrete_boussinesqGreens12;
%     dFT_boussinesqGreens21=dFT_boussinesqGreens12;
%     % nothing to clear here!
%     dFT_boussinesqGreens22=fft2(discrete_boussinesqGreens22);
%     clear discrete_boussinesqGreens22;
%     
%     % Now calculate the solution:                
%     ux_grid=ifft2(dFT_boussinesqGreens11.*dFT_Force_x+dFT_boussinesqGreens12.*dFT_Force_y);
%     clear dFT_boussinesqGreens11 dFT_boussinesqGreens12;
%     uy_grid=ifft2(dFT_boussinesqGreens21.*dFT_Force_x+dFT_boussinesqGreens22.*dFT_Force_y);
%     clear dFT_boussinesqGreens21 dFT_Force_x dFT_boussinesqGreens22 dFT_Force_y;
%     
%     % Now extract the essential part from the solution. It is located in
%     % the center of the padded field.    
%     % I really don't understand why to cut it out like this, but it works!
%     startIndex_x=abs(Nx_G-Nx_F)+1; % Or is it just: startIndex_x=Nx_F
%     startIndex_y=abs(Ny_G-Ny_F)+1;
%     
%     endIndex_x=startIndex_x+Nx_F-1;
%     endIndex_y=startIndex_y+Ny_F-1;
%     
% 
%     % Remove imaginary part caused by round off errors:
%     ux_grid=real(ux_grid(startIndex_x:endIndex_x,startIndex_y:endIndex_y));
%     uy_grid=real(uy_grid(startIndex_x:endIndex_x,startIndex_y:endIndex_y));
% 
% %     figure
% %     imshow(ux_grid,[])
%     % scale the solution appropriately!
%     ux_grid=scalingFactor*ux_grid;
%     uy_grid=scalingFactor*uy_grid;
%     
%     % Now interpolate the displacement field from the regular grid to the irregular
%     % measured grid. Since we have the force field defined on a regular grid
%     % we can use the fast *option. 'linear' is about a factor of two faster than
%     % 'cubic'. Hard to tell if cubic performs better than linear.    
%     if nargin>10 && strcmp(opt,'noIntp')
%         ux = ux_grid;%'*cubic'
%         uy = uy_grid;%'*linear'
%         x_grid=xgrid_F;
%         y_grid=ygrid_F;
%     else
%         ux = interp2(xgrid_F,ygrid_F,ux_grid,x0,y0,'*cubic');%'*cubic'
%         uy = interp2(xgrid_F,ygrid_F,uy_grid,x0,y0,'*cubic');%'*linear'
%         x_grid=x0;
%         y_grid=y0;
%     end
%     
%     %toc;
%     
%     %x0_vec=reshape(x0,[],1);
%     %y0_vec=reshape(y0,[],1);
% 
%     %x_vec=reshape(xgrid_F,[],1);
%     %y_vec=reshape(ygrid_F,[],1);
%     %ux_vec=reshape(ux_grid,[],1);
%     %uy_vec=reshape(uy_grid,[],1);
% 
% 
%     %tic;
%     %[~, ~, ux] = griddata(x_vec,y_vec,ux_vec,x0,y0,'cubic');
%     %[~, ~, uy] = griddata(x_vec,y_vec,uy_vec,x0,y0,'cubic');
%     %toc;
% 
% 
% 
% 
% % to test the example:
% 
% % strange!!!:
% %x0_vec=linspace(-5,20,30);
% %y0_vec=linspace(1,6,15);
% 
% x0_vec=linspace(-50,200,51);
% y0_vec=linspace(0,600,61);
% 
% [x0 y0]=meshgrid(x0_vec,y0_vec);
% 
% E=10000;
% meshPtsFwdSol=2^10;
% 
% xmin=min(x0_vec);
% xmax=max(x0_vec);
% ymin=min(y0_vec);
% ymax=max(y0_vec);
% 
% xmin=0;
% xmax=2;
% ymin=1;
% ymax=3;
% 
% force_x=@(x,y) assumedForce(1,x,y);
% force_y=@(x,y) assumedForce(2,x,y);
% 
% [ux_conv uy_conv]=fwdSolution(x0,y0,E,xmin,xmax,ymin,ymax,force_x,force_y,'conv');
% 
% figure(1)
% quiver(x0,y0,ux_conv,uy_conv);
% 
% 
% [ux_fft uy_fft]=fwdSolution(x0,y0,E,xmin,xmax,ymin,ymax,force_x,force_y,'fft',[],meshPtsFwdSol);
% 
% % figure(11)
% % quiver(x0,y0,ux_fft,uy_fft,'r');
% 
% %compare the two results:
% scalePlot=0.3*sqrt(max(max(ux_conv.^2+uy_conv.^2)));
% figure(40)
% quiver(x0,y0,ux_fft/scalePlot,uy_fft/scalePlot,0,'r');
% hold on
% quiver(x0,y0,ux_conv/scalePlot,uy_conv/scalePlot,0,'g');
% hold off
% 
% %figure(5)
% %quiver(x0,y0,assumedForce(1,x0,y0),assumedForce(2,x0,y0),0,'g');
% 
% corr_x=2*abs(ux_fft-ux_conv)./abs(ux_fft+ux_conv);
% corr_y=2*abs(uy_fft-uy_conv)./abs(uy_fft+uy_conv);
% 
% figure(3)
% surf(x0,y0,corr_x)
% 
% figure(4)
% surf(x0,y0,corr_y)
% 
% display(['mean rel. deviation in %: ',num2str(mean(corr_y(:)))]);
% display([' max rel. deviation in %: ',num2str(max(corr_y(:)))]);
% %uncorr_x=2*abs(ux_fft-ux_conv)./abs(ux_fft+ux_conv)
% %uncorr_y=2*abs(uy_fft-uy_conv)./abs(uy_fft+uy_conv)
% 
% % uncorr_y-corr_y
% 
% display('This should be 0')
% 
% 
% 
% % in case of:
% % x0_vec=linspace(-10,10,25);
% % y0_vec=linspace(-2,2,15);
% % Important note: for a sampling of 2^8 in the fast solution,
% % mean(corr_y(:)) decreases for increasing numerical precision in the
% % direct integration of the fwd solution (up to 'AbsTol' 10^-(8)). Thus, it seems 
% % like as if 2^8 sampling in the Fourier solution is almost as precise as what
% % can be achieved by direct numerical intergeation (up to 'AbsTol' 10^-(10))... which is amazing!!!
% % Fourier sampling: 2^(8)
% % AbsTol:       mean(corr_y(:))
% % 10^(-10)      0.0067
% % 10^(- 9)      0.0070
% % 10^(- 8)      0.0074  (here precision: fast sol ~= direct conv)
% % 10^(- 7)      0.0210
% % 10^(- 6)      0.0442
% % 10^(- 5)      0.1079
% % 10^(- 4)      0.1179
% 
% % Fourier sampling: 2^(10)
% % AbsTol:       mean(corr_y(:))
% % 10^(-10)      0.0011
% % 10^(- 9)      0.0013  (here precision: fast sol ~= direct conv)
% % 10^(- 8)      0.0023
% % 10^(- 7)      0.0162
% % 10^(- 6)      0.0391
% % 10^(- 5)      0.1024
% % 10^(- 4)      0.1123
% 
% % Note also that the best precision is where the displacement is predicted
% % to be highest, that is at the force center! The mean rel. difference
% % between sampling 2^(8) and 2^(10) is 0.0056 which roughly indicates the
% % accuracy of the 2^(8) sampling, which is inline with the results above.

















