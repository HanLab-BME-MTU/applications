function [bead_x,bead_y,bead_ux,bead_uy,Av]=simHydrogelDeformation_CoarseMesh(fx,fy,d,nPoints,E,v,dataPath,multiForce)
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
%       multiForce = true or false, utilize multiple force application areas
%       
%
%   Outputs:
%       bead_x =    x coordinates of bead centers
%       bead_y =    y coordinates of bead centers
%       bead_ux =   x component of bead displacement
%       bead_uy =   y component of bead displacement
%       Av =        vector containing gaussian bead amplitudes
%
%   Example:
%       [bead_x,bead_y,bead_ux,bead_uy,Av]=simHydrogelDeformation_CoarseMesh(0,1000,40,4000,1000,0.49,'/home/sehaarma/Documents/MATLAB/Bead Image Creation',false)

%% //Data path configuration ***********************************************
close all
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
char = 'singleforce_coarsemesh_';
chr1=strcat(char,num2str(E),'E_',num2str(nPoints),'beads_(',num2str(fx),',',num2str(fy),')Pa'); %chr13D='newBeads3D(1kpa,0.5k,1kE)';
refstring=[chrRef chr1 '.tif'];
imgstring=[chrBead chr1 '.tif'];

%% //Generate reference bead image *****************************************
meshPtsFwdSol = 2^9; %number of pix/pts in mesh
beadDiameter = 40; %20nm/bead
pixelSize = 72; %20 nm/pix
thickness = meshPtsFwdSol/2; %thickness of hydrogel for FEM purposes (pixels)
%geometry & image lengths
numPix_x = meshPtsFwdSol; numPix_y = meshPtsFwdSol; %pixel length of image & fem geometry
xmin = 1; ymin = 1;
sigma = 1.68; %stdev of gaussian function
%2D
[refimg,bead_x, bead_y, ~, Av] = simGaussianBeads(numPix_x,numPix_y, sigma, ...
        'npoints', nPoints, 'Border', 'periodic','A',0.4+rand(1,nPoints),...
        'pixelSize',pixelSize,'beadDiameter',beadDiameter);
%3D
% [refimg3D, beadcenters3D, ~, Av3D] = simGaussianSpots3D([meshPtsFwdSol meshPtsFwdSol thickness], ...
%         sigma,'npoints', nPoints, 'Border', 'periodic', 'A',0.5+rand(1,nPoints));
    
%noise addition
refimg = refimg+0.05*rand(numPix_y,numPix_x)*max(refimg(:));
%add 3D noise addition in future

%writing reference image
imwrite(uint16(refimg*2^16/max(max(refimg))),[refPath filesep refstring]);
%add saving 3D refimg
    
%% //Generate force field
if multiForce
    [force_x,force_y] = generateMultiForce();
else
    forceType = 'groupForce';
    [x_mat_u, y_mat_u]=meshgrid(xmin:1:numPix_x,ymin:1:numPix_y);
    force_x = assumedForceAniso2D(1,x_mat_u,y_mat_u,numPix_x/2,numPix_y/2,fx,0,d,d,forceType);
    force_y = assumedForceAniso2D(2,x_mat_u,y_mat_u,numPix_x/2,numPix_y/2,0,fy,d,d,forceType);
end
force_z = [];

%% //Generate model container **********************************************
structModel=createpde('structural','static-solid');

%% //Define geometry *******************************************************
%we use pixels to define fem geometry to ensure consistent units
gm=multicuboid(numPix_x,numPix_y,thickness);
structModel.Geometry = gm;

% //Specify material properties *******************************************
structuralProperties(structModel,'YoungsModulus',E,'PoissonsRatio',v);

% //Apply boundary constraints ********************************************
structuralBC(structModel,'Face',1,'Constraint','fixed'); %New face ID for the bottom.

%% //Apply force distribution on face 2 ************************************
if isempty(force_z)
    force_z=zeros(numPix_x); %can be replaced by real data if a normal force is desired
end
[xLoc,yLoc]=ndgrid(-numPix_x/2:(numPix_x/2)-1,-numPix_x/2:(numPix_x/2)-1);
forceInterpX=griddedInterpolant(xLoc,yLoc,force_x);
forceInterpY=griddedInterpolant(xLoc,yLoc,force_y);
forceInterpZ=griddedInterpolant(xLoc,yLoc,force_z);
%define function handle
hydrogelForce = @(location,state)[forceInterpX(location.x,location.y); ...
                                  forceInterpY(location.x,location.y); ...
                                  forceInterpZ(location.x,location.y);];
%pass function handle and define BCs
structuralBoundaryLoad(structModel,'Face',2,'SurfaceTraction',hydrogelForce,'Vectorize','on'); %New face ID (F8)

%% //Generate mesh for model ***********************************************
generateMesh(structModel,'Hmax',25, 'Hmin',1, 'Hgrad', 1.2);
figure,pdeplot3D(structModel)
%% //Solving structural model **********************************************
structModelResults=solve(structModel);
%outputs
%   displacement = pix
%   stress = kg/(sec^2*pix) = Pa * ((10^9)/72)

%% //Visualizing results ***************************************************
figure
pdeplot3D(structModel,'ColorMapData',structModelResults.Displacement.Magnitude, ...
    'Deformation',structModelResults.Displacement,'DeformationScaleFactor',1)
figure
pdeplot3D(structModel,'ColorMapData',structModelResults.Displacement.uz, ...
    'Deformation',structModelResults.Displacement,'DeformationScaleFactor',1)

%% //Shifting bead locations to apply interpDisp at those locations ********
%2D
%Shifting bead locations because the outputs from simGaussianBeads are not
%centered around the origin which is required for our FEM model.
bead_xshifted=bead_x-(numPix_x/2)+0.5;     %bead_x from simGaussianBeads.m
bead_yshifted=bead_y-(numPix_x/2)+0.5;     %bead_y from simGaussianBeads.m
bead_zshifted=thickness*ones(nPoints,1);   %bead_z can be included in future
%3D
% bead_x3Dshifted=beadcenters3D(:,1)-(meshPtsFwdSol/2)+0.5; %beadcenters3D from simGaussianSpots3D.m
% bead_y3Dshifted=beadcenters3D(:,2)-(meshPtsFwdSol/2)+0.5; %beadcenters3D from simGaussianSpots3D.m
% bead_z3Dshifted=beadcenters3D(:,3)-0.5;   %beadcenters3D from simGaussianSpots3D.m

%% //Interpolating displacements at bead locations using interpDisp ********
%2D
interpDisp=interpolateDisplacement(structModelResults,bead_xshifted,bead_yshifted,bead_zshifted);
%3D
% interpDisp3D=interpolateDisplacement(structModelResults,bead_x3Dshifted,bead_y3Dshifted,bead_z3Dshifted);

%% //Reshaping to match bead location vectors ******************************
%2D
bead_ux = reshape(interpDisp.ux,size(bead_x));
bead_uy = reshape(interpDisp.uy,size(bead_y));
%3D
% bead_ux3D = reshape(interpDisp3D.ux,size(beadcenters3D(:,1)));
% bead_uy3D = reshape(interpDisp3D.uy,size(beadcenters3D(:,2)));
% bead_uz3D = reshape(interpDisp3D.uz,size(beadcenters3D(:,3)));
% beadcenters3Ddisp = [bead_ux3D bead_uy3D bead_uz3D];
% newbeadcenters3D = beadcenters3D + beadcenters3Ddisp;
% xmax3D = ceil(max(newbeadcenters3D(:,1)));
% ymax3D = ceil(max(newbeadcenters3D(:,2)));
% zmax3D = ceil(max(newbeadcenters3D(:,3)));

%% //Plotting interpolated displacement ************************************
figure
q2=quiver(bead_xshifted,bead_yshifted,bead_ux,bead_uy);
xlim([-(meshPtsFwdSol/2) (meshPtsFwdSol/2)]); ylim([-(meshPtsFwdSol/2) (meshPtsFwdSol/2)]);
xlabel('X'); ylabel('Y'); title('Ground Truth Displacement Field - 1 kPa');

%% //Colormapping *********************************************************
%Get the current colormap
% currentColormap = colormap(gca);
% %Now determine the color to make each arrow using a colormap
% [~, ~, ind] = histcounts(magsFEA, size(currentColormap, 1));
% 
% %Now map this to a colormap to get RGB
% cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
% cmap(:,:,4) = 255;
% cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);
% 
% %We repeat each color 3 times (using 1:3 below) because each arrow has 3 vertices
% set(q2.Head, ...
%     'ColorBinding', 'interpolated', ...
%     'ColorData', reshape(cmap(1:3,:,:), [], 4).');   %'
% 
% %We repeat each color 2 times (using 1:2 below) because each tail has 2 vertices
% set(q2.Tail, ...
%     'ColorBinding', 'interpolated', ...
%     'ColorData', reshape(cmap(1:2,:,:), [], 4).');
% 
% set(gca,'Color','k');

%% //Generating Deformed Bead Image ****************************************
%2D
nanElements = sum(isnan(bead_ux));
bead_ux(isnan(bead_ux)) = 0;
bead_uy(isnan(bead_uy)) = 0;
beadimg = simGaussianBeads(numPix_x,numPix_y, sigma, ...
    'x',bead_x + bead_ux,'y',bead_y + bead_uy,'A',Av,...
    'pixelSize',20,'beadDiameter',20,'Border','periodic');
imwrite(uint16(beadimg*2^16/max(max(beadimg))),[imgPath filesep imgstring]);
%3D
% beadimg3D = simGaussianSpots3D([xmax3D ymax3D zmax3D], sigma, ...
%     'X', newbeadcenters3D,'npoints', nPoints, 'A', Av3D,'Border', 'periodic');
%add saving 3D beadimg

%% //Saving outputs for comparison after TFM processing ********************
disp(strcat('Removed ',num2str(nanElements),'NaN elements.'))
save(strcat('outputs_',chr1,'.mat'),'bead_ux','bead_uy','force_x','force_y'); %chr1 defined initially by user, matches image file naming
