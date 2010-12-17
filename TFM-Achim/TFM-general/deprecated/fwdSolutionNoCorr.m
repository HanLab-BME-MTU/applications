function [ux uy]=fwdSolutionNoCorr(x0,y0,E,xmin,xmax,ymin,ymax,force_x,force_y,method)

if nargin<10 || strcmpi(method,'conv')
    tic;
    display('Calulate the convolution explicitely')
    [nRow nCol]=size(x0);

    for i=1:nRow
        for j=1:nCol  
            integrandx = @(x,y) boussinesqGreens(1,1,x0(i,j)-x,y0(i,j)-y,E).*force_x(x,y) + boussinesqGreens(1,2,x0(i,j)-x,y0(i,j)-y,E).*force_y(x,y);
            integrandy = @(x,y) boussinesqGreens(2,1,x0(i,j)-x,y0(i,j)-y,E).*force_x(x,y) + boussinesqGreens(2,2,x0(i,j)-x,y0(i,j)-y,E).*force_y(x,y);

            ux(i,j) = quad2d(integrandx,xmin,xmax,ymin,ymax,'MaxFunEvals',10^4,'AbsTol',5e-9);
            uy(i,j) = quad2d(integrandy,xmin,xmax,ymin,ymax,'MaxFunEvals',10^4,'AbsTol',5e-9);%'AbsTol',5e-11
        end
    end
    toc;
elseif strcmpi(method,'fft')
    display('Use fast convolution')
    
    display('There is still a bug in this code! Debug for non-centered grids! (non-equal spacing seems ok, though?!)')
    

    %***************************************************************
    % Here starts the calculation using the fast fourier transform *
    %***************************************************************
    
    % Number of points to calculate the force field and the Greensfunction.
    % Both Nx and Ny have to be even, otherwise there will be problems with
    % the divergent Greens function at x=0. If Nx and Ny is still chosen to
    % be odd, the Greens function will be cut off at the center (kind of
    % dirty, so better chose Nx and Ny to be even). If Nx and Ny are even
    % though they can not be exactly centered, then the non-zero part is
    % cropped from the convolution result. Shift is about 1/2 of the grid
    % and becomes negligible if Nx and Ny are large anyways. Odd values 
    % for Nx and Ny don't improve the result anyways (no shift required but
    % there is the problem with the cut-off).
    tic;
    
    % This determines the sampling of the force field:
    Nx_F=2^5;
    Ny_F=Nx_F;
    
    % To account for dx*dy in the convolution integral one has to finally
    % rescale the result by the following scaling factor:
    xRange=(max(max(x0))-min(min(x0)));
    yRange=(max(max(y0))-min(min(y0)));
    scalingFactor=(xRange*yRange)/(Nx_F*Ny_F);
    
    % To cover the whole support of the force field, the domain over which
    % the Greensfunctions have to be calculated need to be at least of size:
    % (2Nx-1)x(2Ny-1). Since the number of points have to be even, as
    % discussed above, we do: (2Nx)x(2Ny):
%!!!For this odd/even Problem I have found a solution! See my notebook!

    Nx_G=2*Nx_F;
    Ny_G=2*Ny_F;    
    
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
    
    % create a mesh centered at zero with Nx*Ny meshpoints, where the
    % Greensfunctions are calculated.
    xvec_G=linspace(-xRange,xRange,Nx_G);
    yvec_G=linspace(-yRange,yRange,Ny_G);
    
    [xgrid_G ygrid_G]=meshgrid(xvec_G,yvec_G);
      
    %calculate the force values at the grid_F positions:
    toc
    display('the next one is the interesting time: ')
    tic
    discrete_Force_x=force_x(xgrid_F,ygrid_F); %this has only to be calculated over the support xmin,xmax,ymin,ymax rest is zero
    discrete_Force_y=force_y(xgrid_F,ygrid_F);
    toc
    tic
    %calculate the Greens-function values at the grid_G positions. This can
    %be improved since the Greensfunction never change for a given grid
    %size. So when the Basis functions are calculated this has to be done
    %only once (as well as the FFT for these fields!!!):
    discrete_boussinesqGreens11=boussinesqGreens(1,1,xgrid_G,ygrid_G,E);
    discrete_boussinesqGreens12=boussinesqGreens(1,2,xgrid_G,ygrid_G,E);
   %discrete_boussinesqGreens21=discrete_boussinesqGreens12;
    discrete_boussinesqGreens22=boussinesqGreens(2,2,xgrid_G,ygrid_G,E);
    
    % Pad the calculated fields with zero to the next power larger than 
    % (2*N-1), see above. For this setup, the FFT is fastest. In stead of
    % 0 one should be able to also use symmetric but then, u_x is
    % essentially zero, is there an error in the code?!
    discrete_Force_x=padarray(discrete_Force_x,[Nx_pad-Nx_F Ny_pad-Ny_F],0,'post');%'symmetric','post');
    discrete_Force_y=padarray(discrete_Force_y,[Nx_pad-Nx_F Ny_pad-Ny_F],0,'post');
    
    discrete_boussinesqGreens11=padarray(discrete_boussinesqGreens11,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
    discrete_boussinesqGreens12=padarray(discrete_boussinesqGreens12,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
   %discrete_boussinesqGreens22=padarray(discrete_boussinesqGreens22,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post'); 
    discrete_boussinesqGreens22=padarray(discrete_boussinesqGreens22,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
    
    %Now calculate the fourier transforms:
    dFT_Force_x=fft2(discrete_Force_x);
    dFT_Force_y=fft2(discrete_Force_y);
    
    %This has to be calculated only once for all basis functions!
    dFT_boussinesqGreens11=fft2(discrete_boussinesqGreens11);
    dFT_boussinesqGreens12=fft2(discrete_boussinesqGreens12);
    dFT_boussinesqGreens21=dFT_boussinesqGreens12;
    dFT_boussinesqGreens22=fft2(discrete_boussinesqGreens22);
    
    % Now calculate the solution:                
    ux_grid=ifft2(dFT_boussinesqGreens11.*dFT_Force_x+dFT_boussinesqGreens12.*dFT_Force_y);
    uy_grid=ifft2(dFT_boussinesqGreens21.*dFT_Force_x+dFT_boussinesqGreens22.*dFT_Force_y);

    % Plot the solution:
%     figure(10)
%     surf(ux_grid)
%     
%     figure(11)
%     surf(uy_grid)
    
    % Now extract the essential part from the solution. It is located in
    % the center of the padded field. Here, Nx and Ny have to be even!
    % Error otherwise. Also remove imaginary parts:
    
%     startIndex_x=floor((Nx_pad-Nx_F)/2);
%     if mod(Nx_F,2)==0
%         startIndex_x=startIndex_x+1;
%     end
%     endIndex_x  =startIndex_x+Nx_F-1;
%     
%     startIndex_y=floor((Ny_pad-Ny_F)/2);%+1;
%     if mod(Ny_F,2)==0
%         startIndex_y=startIndex_y+1;
%     end
%     endIndex_y  =startIndex_y+Ny_F-1;
    
    % I really don't understand why to cut it out like this, but it works!
    startIndex_x=abs(Nx_G-Nx_F);
    startIndex_y=abs(Ny_G-Ny_F);
    
    endIndex_x=startIndex_x+Nx_F-1;
    endIndex_y=startIndex_y+Ny_F-1;
    

    %Remove imaginary part caused by round off errors:
    ux_grid=real(ux_grid(startIndex_x:endIndex_x,startIndex_y:endIndex_y));
    uy_grid=real(uy_grid(startIndex_x:endIndex_x,startIndex_y:endIndex_y));

    % Now interpolate the displacement field from the regular grid to the irregular
    % measured grid. Since we have the force field defined on a regular grid
    % we can use the fast *option. 'linear' is about a factor of two faster than
    % 'cubic'. Hard to tell if cubic performs better than linear.
      
    ux =scalingFactor*interp2(xgrid_F,ygrid_F,ux_grid,x0,y0,'*cubic');%'*cubic'
    uy =scalingFactor*interp2(xgrid_F,ygrid_F,uy_grid,x0,y0,'*cubic');%'*linear' 
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

return;

% to test the example:

x0_vec=linspace(-10,2,25);
y0_vec=linspace(-5,2,15);

[x0 y0]=meshgrid(x0_vec,y0_vec);

E=10000;

xmin=min(x0_vec);
xmax=max(x0_vec);
ymin=min(y0_vec);
ymax=max(y0_vec);

force_x=@(x,y) assumedForce(1,x,y);
force_y=@(x,y) assumedForce(2,x,y);

[ux uy]=fwdSolution(x0,y0,E,xmin,xmax,ymin,ymax,force_x,force_y);

%figure(1)
%quiver(x0,y0,ux,uy);


[ux_conv uy_conv]=fwdSolutionElaborate(x0,y0,E,xmin,xmax,ymin,ymax,force_x,force_y,'fft');

% figure(11)
% quiver(x0,y0,ux_conv,uy_conv,'r');

%compare the two results:
figure(40)
quiver(x0,y0,ux_conv,uy_conv,'r');
hold on
quiver(x0,y0,ux,uy,'g');
hold off

figure(5)
quiver(x0,y0,assumedForce(1,x0,y0),assumedForce(2,x0,y0),'g');

ux_conv./(ux_conv+ux)
uy_conv./(uy_conv+uy)
display('This should be 1/2 which is not fulfilled for a non-centered grid! at small displacements?!')














