function [ux uy uz x_grid y_grid z0]=fwdSolution3D(x0,y0,z0,E,v,xmin,xmax,ymin,ymax,force_x,force_y,force_z,method,opt)
% This calculates the forward solution for the general problem of a force
% acting on the top surface, returning the displacement at the single plane
% z0. x0, y0 are matrices but z0 is allowed to be a single value only!

if isempty(force_z)
    force_z = @(x,y) zeros(size(x));
end

if nargin<14 || strcmpi(method,'conv')
    display('This part of the code has not been tested yet!')
    tic;
    display('Calulate the convolution explicitely')
    [nRow nCol]=size(x0);

    for i=1:nRow
        for j=1:nCol  
            integrandx = @(x,y) boussinesqGreens3D(1,1,x0(i,j)-x,y0(i,j)-y,z0,E,v).*force_x(x,y) + boussinesqGreens3D(1,2,x0(i,j)-x,y0(i,j)-y,z0,E,v).*force_y(x,y) + boussinesqGreens3D(1,3,x0(i,j)-x,y0(i,j)-y,z0,E,v).*force_z(x,y);
            integrandy = @(x,y) boussinesqGreens3D(2,1,x0(i,j)-x,y0(i,j)-y,z0,E,v).*force_x(x,y) + boussinesqGreens3D(2,2,x0(i,j)-x,y0(i,j)-y,z0,E,v).*force_y(x,y) + boussinesqGreens3D(2,3,x0(i,j)-x,y0(i,j)-y,z0,E,v).*force_z(x,y);
            integrandz = @(x,y) boussinesqGreens3D(3,1,x0(i,j)-x,y0(i,j)-y,z0,E,v).*force_x(x,y) + boussinesqGreens3D(3,2,x0(i,j)-x,y0(i,j)-y,z0,E,v).*force_y(x,y) + boussinesqGreens3D(3,3,x0(i,j)-x,y0(i,j)-y,z0,E,v).*force_z(x,y);

            ux(i,j) = quad2d(integrandx,xmin,xmax,ymin,ymax,'MaxFunEvals',10^4,'AbsTol',5e-9);
            uy(i,j) = quad2d(integrandy,xmin,xmax,ymin,ymax,'MaxFunEvals',10^4,'AbsTol',5e-9);%'AbsTol',5e-11
            uz(i,j) = quad2d(integrandz,xmin,xmax,ymin,ymax,'MaxFunEvals',10^4,'AbsTol',5e-9);
        end
    end
    toc;
elseif strcmpi(method,'fft')
    %display('Use fast convolution')    

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
    % correction term which is of particular cimportance for sparce
    % sampling, meaning that Nx_F is small. This alogorithm performs very
    % well and has been cross-validated with the results obtained using the
    % 'conv' option. If you want to repeat the test, use the few lines of
    % code at the very end of this function.
    
    tic;
    
    % This determines the sampling of the force field:
    Nx_F=2^8; % 2^10 is the most dense sampling possible.
    Ny_F=Nx_F;
    
    % To account for dx*dy in the convolution integral one has to finally
    % rescale the result by the following scaling factor:
    xRange=(max(max(x0))-min(min(x0)));
    yRange=(max(max(y0))-min(min(y0)));
    scalingFactor=(xRange*yRange)/(Nx_F*Ny_F);
    
    % To cover the whole support of the force field, the domain over which
    % the Greensfunctions have to be calculated need to be at least of size:
    % (2Nx-1)x(2Ny-1).

    Nx_G=2*Nx_F-1;
    Ny_G=2*Ny_F-1;    
    
    % Subsequently, these have to be padded with zeros according to:
    Nx_pad=Nx_F+Nx_G-1;
    Ny_pad=Ny_F+Ny_G-1;
    
    % These might not be a power of 2, make sure that they are:
    Nx_pad=pow2(nextpow2(Nx_pad));
    Ny_pad=pow2(nextpow2(Ny_pad));

    % First determine the boundaries of the mesh:
    leftUpperCorner =[min(min(x0)) min(min(y0))];
    rightLowerCorner=[max(max(x0)) max(max(y0))];

    % create a regular mesh with Nx*Ny meshpoints where the force field is
    % calculated. This need not to be a power of 2 yet:
    xvec_F=linspace(leftUpperCorner(1),rightLowerCorner(1),Nx_F);
    yvec_F=linspace(leftUpperCorner(2),rightLowerCorner(2),Ny_F);
    [xgrid_F ygrid_F]=meshgrid(xvec_F,yvec_F);
    
    % create a mesh centered at zero with Nx_G*Ny_G meshpoints, where the
    % Greensfunctions are calculated.
    xvec_G=linspace(-xRange,xRange,Nx_G);
    yvec_G=linspace(-yRange,yRange,Ny_G);
    
    [xgrid_G ygrid_G]=meshgrid(xvec_G,yvec_G);
      
    %calculate the force values at the grid_F positions:
    discrete_Force_x=force_x(xgrid_F,ygrid_F); %this has only to be calculated over the support xmin,xmax,ymin,ymax rest is zero
    discrete_Force_y=force_y(xgrid_F,ygrid_F);
    discrete_Force_z=force_z(xgrid_F,ygrid_F);

    %calculate the Greens-function values at the grid_G positions. This can
    %be improved since the Greensfunction never change for a given grid
    %size. So when the Basis functions are calculated this has to be done
    %only once (as well as the FFT for these fields!!!):
    
    % it should be checked first, if all Greensfunctions are really needed!
    discrete_boussinesqGreens3D11=boussinesqGreens3D(1,1,xgrid_G,ygrid_G,z0,E,v);
    discrete_boussinesqGreens3D12=boussinesqGreens3D(1,2,xgrid_G,ygrid_G,z0,E,v);
    discrete_boussinesqGreens3D13=boussinesqGreens3D(1,3,xgrid_G,ygrid_G,z0,E,v);
    discrete_boussinesqGreens3D21=discrete_boussinesqGreens3D12;
    discrete_boussinesqGreens3D22=boussinesqGreens3D(2,2,xgrid_G,ygrid_G,z0,E,v);
    discrete_boussinesqGreens3D23=boussinesqGreens3D(2,3,xgrid_G,ygrid_G,z0,E,v);
    discrete_boussinesqGreens3D31=boussinesqGreens3D(3,1,xgrid_G,ygrid_G,z0,E,v);
    discrete_boussinesqGreens3D32=boussinesqGreens3D(3,2,xgrid_G,ygrid_G,z0,E,v);
    discrete_boussinesqGreens3D33=boussinesqGreens3D(3,3,xgrid_G,ygrid_G,z0,E,v);
       
    % Pad the calculated fields with zero to the next power larger than 
    % (2*N-1), see above. For this setup, the FFT is fastest.
    discrete_Force_x=padarray(discrete_Force_x,[Nx_pad-Nx_F Ny_pad-Ny_F],0,'post');%'symmetric','post');
    discrete_Force_y=padarray(discrete_Force_y,[Nx_pad-Nx_F Ny_pad-Ny_F],0,'post');
    discrete_Force_z=padarray(discrete_Force_z,[Nx_pad-Nx_F Ny_pad-Ny_F],0,'post');
    
    discrete_boussinesqGreens3D11=padarray(discrete_boussinesqGreens3D11,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
    discrete_boussinesqGreens3D12=padarray(discrete_boussinesqGreens3D12,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
    discrete_boussinesqGreens3D13=padarray(discrete_boussinesqGreens3D13,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
    discrete_boussinesqGreens3D21=discrete_boussinesqGreens3D12;
    discrete_boussinesqGreens3D22=padarray(discrete_boussinesqGreens3D22,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
    discrete_boussinesqGreens3D23=padarray(discrete_boussinesqGreens3D23,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
    discrete_boussinesqGreens3D31=padarray(discrete_boussinesqGreens3D31,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
    discrete_boussinesqGreens3D32=padarray(discrete_boussinesqGreens3D32,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
    discrete_boussinesqGreens3D33=padarray(discrete_boussinesqGreens3D33,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
    
    %Now calculate the fourier transforms:
    dFT_Force_x=fft2(discrete_Force_x);
    dFT_Force_y=fft2(discrete_Force_y);
    dFT_Force_z=fft2(discrete_Force_z);
    
    %This has to be calculated only once for all basis functions!
    dFT_boussinesqGreens3D11=fft2(discrete_boussinesqGreens3D11);
    dFT_boussinesqGreens3D12=fft2(discrete_boussinesqGreens3D12);
    dFT_boussinesqGreens3D13=fft2(discrete_boussinesqGreens3D13);
    dFT_boussinesqGreens3D21=dFT_boussinesqGreens3D12;
    dFT_boussinesqGreens3D22=fft2(discrete_boussinesqGreens3D22);
    dFT_boussinesqGreens3D23=fft2(discrete_boussinesqGreens3D23);
    dFT_boussinesqGreens3D31=fft2(discrete_boussinesqGreens3D31);
    dFT_boussinesqGreens3D32=fft2(discrete_boussinesqGreens3D32);
    dFT_boussinesqGreens3D33=fft2(discrete_boussinesqGreens3D33);
    
    % Now calculate the solution:                
    ux_grid=ifft2(dFT_boussinesqGreens3D11.*dFT_Force_x+dFT_boussinesqGreens3D12.*dFT_Force_y+dFT_boussinesqGreens3D13.*dFT_Force_z);
    uy_grid=ifft2(dFT_boussinesqGreens3D21.*dFT_Force_x+dFT_boussinesqGreens3D22.*dFT_Force_y+dFT_boussinesqGreens3D23.*dFT_Force_z);
    uz_grid=ifft2(dFT_boussinesqGreens3D31.*dFT_Force_x+dFT_boussinesqGreens3D32.*dFT_Force_y+dFT_boussinesqGreens3D33.*dFT_Force_z);

    % Plot the solution:
%     figure(10)
%     surf(ux_grid)
%     
%     figure(11)
%     surf(uy_grid)
    
    % Now extract the essential part from the solution. It is located in
    % the center of the padded field.    
    % I really don't understand why to cut it out like this, but it works!
    startIndex_x=abs(Nx_G-Nx_F)+1; % Or is it just: startIndex_x=Nx_F
    startIndex_y=abs(Ny_G-Ny_F)+1;
    
    endIndex_x=startIndex_x+Nx_F-1;
    endIndex_y=startIndex_y+Ny_F-1;
    

    % Remove imaginary part caused by round off errors:
    ux_grid=real(ux_grid(startIndex_x:endIndex_x,startIndex_y:endIndex_y));
    uy_grid=real(uy_grid(startIndex_x:endIndex_x,startIndex_y:endIndex_y));
    uz_grid=real(uz_grid(startIndex_x:endIndex_x,startIndex_y:endIndex_y));
    
    % Add the solution for G(0,0). This is a correction term which becomes
    % irrelevant for very dense sampling. But for small Nx_F it is REALLY
    % essential!
%!!! This assumes that the grid in x and y are equally spaced!
    
    display('The correction term has to be recalculated!')
    r=1; % what is r???
    corrTermXY =(1+v)*(2-v)/(2*pi*E)*(4*r*log(3+2*sqrt(2)));
    corrTermZ  =2*(1-v)/(2-v)*corrTermXY;
    
    % There is only one integral for each component to be calculated, all
    % other intergrals are zero, either because z=0 or because of symmetry:
    ux_grid=ux_grid+force_x(xgrid_F,ygrid_F)*corrTermXY;
    uy_grid=uy_grid+force_y(xgrid_F,ygrid_F)*corrTermXY;
    uz_grid=uz_grid+force_z(xgrid_F,ygrid_F)*corrTermZ;

    % Now interpolate the displacement field from the regular grid to the irregular
    % measured grid. Since we have the force field defined on a regular grid
    % we can use the fast *option. 'linear' is about a factor of two faster than
    % 'cubic'. Hard to tell if cubic performs better than linear.
    
    if nargin==14 && strcmp(opt,'noIntp')
        ux =scalingFactor*(ux_grid);%'*cubic'
        uy =scalingFactor*(uy_grid);%'*linear'
        uz =scalingFactor*(uz_grid);
        x_grid=xgrid_F;
        y_grid=ygrid_F;
    else
        ux =scalingFactor*interp2(xgrid_F,ygrid_F,ux_grid,x0,y0,'*cubic');%'*cubic'
        uy =scalingFactor*interp2(xgrid_F,ygrid_F,uy_grid,x0,y0,'*cubic');%'*linear'
        uz =scalingFactor*interp2(xgrid_F,ygrid_F,uz_grid,x0,y0,'*cubic');
        x_grid=x0;
        y_grid=y0;
    end
    
    toc;
    
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
end














