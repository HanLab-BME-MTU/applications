function [ux uy x_grid y_grid]=fwdSolution_New(x0,y0,E,xmin,xmax,ymin,ymax,force_x,force_y,method,opt)
% This forward solution is only valid for a Poisson's ratio v=0.5!
% Input: No matter what the dimension of x0 and y0 is (pix, or um), the
%        dimension of the surface stresses (force_x, force_y) must have the
%        same dimension as the elastic modulus E, usually Pa.

% Output: The calculated ux and uy have the same dimension as the input
%         x0, y0.

if nargin<10 || strcmpi(method,'conv')
    tic;
    display('Calulate the convolution explicitely')
    [nRow nCol]=size(x0);

    for i=1:nRow
        for j=1:nCol  
            integrandx = @(x,y) boussinesqGreens(1,1,x0(i,j)-x,y0(i,j)-y,E).*force_x(x,y) + boussinesqGreens(1,2,x0(i,j)-x,y0(i,j)-y,E).*force_y(x,y);
            integrandy = @(x,y) boussinesqGreens(2,1,x0(i,j)-x,y0(i,j)-y,E).*force_x(x,y) + boussinesqGreens(2,2,x0(i,j)-x,y0(i,j)-y,E).*force_y(x,y);

            ux(i,j) = quad2d(integrandx,xmin,xmax,ymin,ymax,'MaxFunEvals',10^5,'AbsTol',5e-7);
            uy(i,j) = quad2d(integrandy,xmin,xmax,ymin,ymax,'MaxFunEvals',10^5,'AbsTol',5e-7);%'AbsTol',5e-11
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
    
    %tic;
    
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
    discrete_Force_x_unPadded=force_x(xgrid_F,ygrid_F); %this has only to be calculated over the support xmin,xmax,ymin,ymax rest is zero
    discrete_Force_y_unPadded=force_y(xgrid_F,ygrid_F);

    % Calculate the Greens-function values at the grid_G positions. This can
    % be improved since the Greensfunction never change for a given grid
    % size. When the Basis functions are calculated this has to be done
    % only once (as well as the FFT for these fields!!!):
    discrete_boussinesqGreens11=boussinesqGreens(1,1,xgrid_G,ygrid_G,E);
    discrete_boussinesqGreens12=boussinesqGreens(1,2,xgrid_G,ygrid_G,E);
   %discrete_boussinesqGreens21=discrete_boussinesqGreens12;
    discrete_boussinesqGreens22=boussinesqGreens(2,2,xgrid_G,ygrid_G,E);
    
    % Pad the calculated fields with zero to the next power larger than 
    % (2*N-1), see above. For this setup, the FFT is fastest.
    discrete_Force_x=padarray(discrete_Force_x_unPadded,[Nx_pad-Nx_F Ny_pad-Ny_F],0,'post');%'symmetric','post');
    discrete_Force_y=padarray(discrete_Force_y_unPadded,[Nx_pad-Nx_F Ny_pad-Ny_F],0,'post');
    
    discrete_boussinesqGreens11=padarray(discrete_boussinesqGreens11,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
    discrete_boussinesqGreens12=padarray(discrete_boussinesqGreens12,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
   %discrete_boussinesqGreens22=padarray(discrete_boussinesqGreens22,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post'); 
    discrete_boussinesqGreens22=padarray(discrete_boussinesqGreens22,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
    
    % Now calculate the fourier transforms:
    dFT_Force_x=fft2(discrete_Force_x);
    clear discrete_Force_x;
    dFT_Force_y=fft2(discrete_Force_y);
    clear discrete_Force_y;
    
    % This has to be calculated only once for all basis functions!
    dFT_boussinesqGreens11=fft2(discrete_boussinesqGreens11);
    clear discrete_boussinesqGreens11;
    dFT_boussinesqGreens12=fft2(discrete_boussinesqGreens12);
    clear discrete_boussinesqGreens12;
    dFT_boussinesqGreens21=dFT_boussinesqGreens12;
    % nothing to clear here!
    dFT_boussinesqGreens22=fft2(discrete_boussinesqGreens22);
    clear discrete_boussinesqGreens22;
    
    dx=xRange;%/Nx_pad;
    dy=yRange;%/Ny_pad;
    %kx_vec = 1/dx.*(-Nx_pad/2:Nx_pad/2);
    kx_vec = (2*pi/Nx_pad/dx)*[(-Nx_pad/2:-1) (1:Nx_pad/2)];
    ky_vec = (2*pi/Nx_pad/dy)*[(-Ny_pad/2:-1) (1:Ny_pad/2)];
    [kx_grid ky_grid]=meshgrid(kx_vec,ky_vec);
    
    figure(1)
    surf(real(abs(fftshift(dFT_boussinesqGreens22))));
    
    figure(2)
    surf(real(FTBoussinesqGreens(2,2,kx_grid,ky_grid,E)));
    % kx_vec = 2*pi/i_max/cluster_size.*[0:(i_max/2-1) (-i_max/2:-1)];
    % ky_vec = 2*pi/j_max/cluster_size.*[0:(j_max/2-1)
    % (-j_max/2:-1)];
    
    figure(3); surf(abs(real(fftshift(dFT_boussinesqGreens22)))./real(FTBoussinesqGreens(2,2,kx_grid,ky_grid,E)));
    
    % Now calculate the k-vector!
    
    % Now calculate the solution:                
    ux_grid=ifft2(dFT_boussinesqGreens11.*dFT_Force_x+dFT_boussinesqGreens12.*dFT_Force_y);
    clear dFT_boussinesqGreens11 dFT_boussinesqGreens12;
    uy_grid=ifft2(dFT_boussinesqGreens21.*dFT_Force_x+dFT_boussinesqGreens22.*dFT_Force_y);
    clear dFT_boussinesqGreens21 dFT_Force_x dFT_boussinesqGreens22 dFT_Force_y;

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

%!!! This could be improved by using the analytical solution for the Fourie
%!!! Transform of the Greensfunction!
    % Add the solution for G(0,0). This is a correction term which becomes
    % irrelevant for very dense sampling. But for small Nx_F it is REALLY
    % essential!
    % Set the Poisson's ratio to 0.5:
    v=0.5;
    r=1; % what is r???
    corrTerm=(1+v)*(2-v)/(2*pi*E)*(4*r*log(3+2*sqrt(2)));
    
    ux_grid=ux_grid+discrete_Force_x_unPadded*corrTerm;
    clear discrete_Force_x_unPadded;
    uy_grid=uy_grid+discrete_Force_y_unPadded*corrTerm;
    clear discrete_Force_y_unPadded;

    % Now interpolate the displacement field from the regular grid to the irregular
    % measured grid. Since we have the force field defined on a regular grid
    % we can use the fast *option. 'linear' is about a factor of two faster than
    % 'cubic'. Hard to tell if cubic performs better than linear.
    
    if nargin==11 && strcmp(opt,'noIntp')
        ux =scalingFactor*(ux_grid);%'*cubic'
        uy =scalingFactor*(uy_grid);%'*linear'
        x_grid=xgrid_F;
        y_grid=ygrid_F;
    else
        ux =scalingFactor*interp2(xgrid_F,ygrid_F,ux_grid,x0,y0,'*cubic');%'*cubic'
        uy =scalingFactor*interp2(xgrid_F,ygrid_F,uy_grid,x0,y0,'*cubic');%'*linear'
        x_grid=x0;
        y_grid=y0;
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
end

return;

% to test the example:

x0_vec=linspace(-10,10,15);
y0_vec=linspace(-2,2,10);

[x0 y0]=meshgrid(x0_vec,y0_vec);

E=10000;

xmin=min(x0_vec);
xmax=max(x0_vec);
ymin=min(y0_vec);
ymax=max(y0_vec);

force_x=@(x,y) assumedForce(1,x,y);
force_y=@(x,y) assumedForce(2,x,y);

[ux uy]=fwdSolution(x0,y0,E,xmin,xmax,ymin,ymax,force_x,force_y,'conv');

figure(1)
quiver(x0,y0,ux,uy);


[ux_conv uy_conv]=fwdSolution_New(x0,y0,E,xmin,xmax,ymin,ymax,force_x,force_y,'fft');

% figure(11)
% quiver(x0,y0,ux_conv,uy_conv,'r');

%compare the two results:
scalePlot=15*sqrt(max(max(ux.^2+ux.^2)))
figure(40)
quiver(x0,y0,ux_conv/scalePlot,uy_conv/scalePlot,0,'r');
hold on
quiver(x0,y0,ux/scalePlot,uy/scalePlot,0,'g');
hold off

%figure(5)
%quiver(x0,y0,assumedForce(1,x0,y0),assumedForce(2,x0,y0),0,'g');

2*abs(ux_conv-ux)./abs(ux_conv+ux)
2*abs(uy_conv-uy)./abs(uy_conv+uy)
display('This should be 1/2 which is not fulfilled for a non-centered grid! at small displacements?!')














