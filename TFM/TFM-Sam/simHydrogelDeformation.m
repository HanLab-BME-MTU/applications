function [bead_x,bead_y,bead_ux,bead_uy,Av]=simHydrogelDeformation(fx,fy,d,nPoints,E,v,dataPath)
% SIMHYDROGELDEFORMATION generates a 2D gaussian reference bead
% image and then uses the user defined force, area, young's modulus, and 
% poisson's ratio to deform the image using FEA
%
%   Inputs:
%       fx =        X component of force magnitude (Pa)
%       fy =        Y component of force magnitude (Pa)
%       d =         Force application area diameter (circular)
%       nPoints =   Number of Gaussian beads to simulate
%       E =         Young's Modulus (Pa) - usually 1000
%       v =         Poisson's Ratio - usually 0.49
%       dataPath =  File directory for bead images
%
%   Outputs:
%       bead_x =    x coordinates of bead centers
%       bead_y =    y coordinates of bead centers
%       bead_ux =   x component of bead displacement
%       bead_uy =   y component of bead displacement
%       Av =        vector containing gaussian bead amplitudes
%
%   Example:
%       [bead_x,bead_y,bead_ux,bead_uy,Av]=simHydrogelDeformation(0,1000,40,5000,1000,0.49,'C:\Users\Sam Haarman\Documents\MATLAB\Han Lab\Bead Image Creation')

% //Data path configuration ***********************************************
imgPath=[dataPath filesep 'Beads'];
refPath=[dataPath filesep 'Reference'];
orgPath=[dataPath filesep 'Original'];
if ~exist(refPath,'dir') || ~exist(orgPath,'dir')
    mkdir(imgPath);
    mkdir(refPath);
    mkdir(orgPath);
end

chrRef='img1ref_';
chrBead='img3bead_';
%Change these character vectors to rename files for next run
chr1='1kPa5kbeads'; chr13D='1kPa5kbeads3D';
refstring=[chrRef chr1 '.tif'];
imgstring=[chrBead chr1 '.tif'];

% //Generate reference bead image *****************************************
meshPtsFwdSol=2^8;
xmax=meshPtsFwdSol; xmin=1;
ymax=meshPtsFwdSol; ymin=1;
sigma = 1.68;
%2D
[refimg,bead_x, bead_y, ~, Av] = simGaussianBeads(xmax,ymax, sigma, ...
        'npoints', nPoints, 'Border', 'periodic','A',0.3+rand(1,nPoints));
%3D
[refimg3D, beadcenters3D, ~, Av3D] = simGaussianSpots3D([256 256 15], ...
        sigma,'npoints', nPoints, 'Border', 'periodic', 'A',0.3+rand(1,nPoints));
    
%noise addition
refimg = refimg+0.05*rand(ymax,xmax)*max(refimg(:));
%add 3D noise addition

%saving and writing reference image
save(append(chrRef,chr1),'refimg');
save(append(chrRef,chr13D),'refimg3D')
imwrite(uint16(refimg*2^16/max(max(refimg))),[refPath filesep refstring]);
    
% //Generate force field
forceType = 'groupForce';
[x_mat_u, y_mat_u]=meshgrid(xmin:1:xmax,ymin:1:ymax);
force_x = assumedForceAniso2D(1,x_mat_u,y_mat_u,(xmax/2),ymax/2,0,fx,d,d,forceType);
force_y = assumedForceAniso2D(2,x_mat_u,y_mat_u,(xmax/2),ymax/2,0,fy,d,d,forceType);
force_z = [];


% //Generate model container **********************************************
structModel=createpde('structural','static-solid');

% //Define geometry *******************************************************
%reduce to micron size geometry, 72 nm/pix
gm=multicuboid(256,256,15);
structModel.Geometry=gm;

% //Specify material properties *******************************************
structuralProperties(structModel,'YoungsModulus',E,'PoissonsRatio',v);

% //Apply boundary constraints ********************************************
structuralBC(structModel,'Face',1,'Constraint','fixed');

% //Apply force distribution on face 2 ************************************
if isempty(force_z)
force_z=zeros(256); %can be replaced by real data if a normal force is desired
end
[xLoc,yLoc]=ndgrid(-128:127,-128:127);
forceInterpX=griddedInterpolant(xLoc,yLoc,force_x);
forceInterpY=griddedInterpolant(xLoc,yLoc,force_y);
forceInterpZ=griddedInterpolant(xLoc,yLoc,force_z);
%define function handle
hydrogelForce = @(location,state)[forceInterpX(location.x,location.y); ...
                                  forceInterpY(location.x,location.y); ...
                                  forceInterpZ(location.x,location.y)];
%pass function handle and define BCs
structuralBoundaryLoad(structModel,'Face',2,'SurfaceTraction',hydrogelForce,'Vectorize','on');

% //Generate mesh for model ***********************************************
generateMesh(structModel,'Hmax',10);

% //Solving structural model **********************************************
structModelResults=solve(structModel);

% //Visualizing results ***************************************************
figure(2)
pdeplot3D(structModel,'ColorMapData',structModelResults.Displacement.Magnitude, ...
    'Deformation',structModelResults.Displacement)
figure(4)
pdeplot3D(structModel,'ColorMapData',structModelResults.Displacement.uz, ...
    'Deformation',structModelResults.Displacement)

% //Shifting bead locations to apply interpDisp at those locations ********
%2D
bead_xshifted=bead_x-128.5;     %bead_x from simGaussianBeads.m
bead_yshifted=bead_y-128.5;     %bead_y from simGaussianBeads.m
bead_zshifted=15*ones(nPoints,1); %bead_z can be included in future
%3D
bead_x3Dshifted=beadcenters3D(:,1)-128.5; %beadcenters3D from simGaussianSpots3D.m
bead_y3Dshifted=beadcenters3D(:,2)-128.5; %beadcenters3D from simGaussianSpots3D.m
bead_z3Dshifted=beadcenters3D(:,3)-0.5;   %beadcenters3D from simGaussianSpots3D.m

% //Interpolating displacements at bead locations using interpDisp ********
%2D
interpDisp=interpolateDisplacement(structModelResults,bead_xshifted,bead_yshifted,bead_zshifted);
%3D
interpDisp3D=interpolateDisplacement(structModelResults,bead_x3Dshifted,bead_y3Dshifted,bead_z3Dshifted);

% //Reshaping to match bead location vectors ******************************
%2D
bead_ux = reshape(interpDisp.ux,size(bead_x));
bead_uy = reshape(interpDisp.uy,size(bead_y));
magsFEA = sqrt((bead_ux.^2)+(bead_uy.^2));
%3D
bead_ux3D = reshape(interpDisp3D.ux,size(beadcenters3D(:,1)));
bead_uy3D = reshape(interpDisp3D.uy,size(beadcenters3D(:,2)));
bead_uz3D = reshape(interpDisp3D.uz,size(beadcenters3D(:,3)));
beadcenters3Ddisp = [bead_ux3D bead_uy3D bead_uz3D];
newbeadcenters3D = beadcenters3D + beadcenters3Ddisp;
xmax3D = ceil(max(newbeadcenters3D(:,1)));
ymax3D = ceil(max(newbeadcenters3D(:,2)));
zmax3D = ceil(max(newbeadcenters3D(:,3)));

% //Plotting interpolated displacement ************************************
figure(3)
q2=quiver(bead_xshifted,bead_yshifted,bead_ux,bead_uy);
xlim([-128 128]); ylim([-128 128]);
xlabel('X'); ylabel('Y'); title('Ground Truth Displacement Field - 1 kPa');

%\\ Colormapping **********************************************************
%Get the current colormap
currentColormap = colormap(gca);
%Now determine the color to make each arrow using a colormap
[~, ~, ind] = histcounts(magsFEA, size(currentColormap, 1));

%Now map this to a colormap to get RGB
cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
cmap(:,:,4) = 255;
cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);

%We repeat each color 3 times (using 1:3 below) because each arrow has 3 vertices
set(q2.Head, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:3,:,:), [], 4).');   %'

%We repeat each color 2 times (using 1:2 below) because each tail has 2 vertices
set(q2.Tail, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:2,:,:), [], 4).');

set(gca,'Color','k');

%\\ Generating Deformed Bead Image ****************************************
%2D
beadimg = simGaussianBeads(xmax,ymax, sigma, ...
    'x',bead_x + bead_ux,'y',bead_y + bead_uy,'A',Av,'Border', 'padded');
save(append(chrBead,chr1),'beadimg');
imwrite(uint16(beadimg*2^16/max(max(beadimg))),[imgPath filesep imgstring]);
%3D
beadimg3D = simGaussianSpots3D([xmax3D ymax3D zmax3D], sigma, ...
    'X', newbeadcenters3D,'npoints', nPoints, 'A', Av3D, 'Border', 'padded', 'Verbose' ,'on');
save(append(chrBead,chr13D),'beadimg3D');