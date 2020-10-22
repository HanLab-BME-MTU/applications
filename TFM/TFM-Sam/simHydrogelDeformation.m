function [bead_x,bead_y,bead_ux,bead_uy,Av]=simHydrogelDeformation(fx,fy,d,nPoints,E,v,dataPath,multiForce)
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
%       [bead_x,bead_y,bead_ux,bead_uy,Av]=simHydrogelDeformation(0,1000,40,4000,1000,0.49,'/home/sehaarma/Documents/MATLAB/Bead Image Creation',false)

%% //Data path configuration ***********************************************
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
chr1='singleforce_finemesh_1kE_4kbeads'; %chr13D='newBeads3D(1kpa,0.5k,1kE)';
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
% gm=multicuboid(numPix_x,numPix_y,thickness);
% We want to now insert a sphere to ensure the refined mesh near force
% application area
% Create a 3-D rectangular mesh grid first (with sparse spacing)
[xg, yg, zg] = meshgrid(linspace(-numPix_x/2,(numPix_x/2)-1,meshPtsFwdSol/2^4),...
                        linspace(-numPix_x/2,(numPix_x/2)-1,meshPtsFwdSol/2^4),...
                        -thickness:0); %/2^4 defines sparsity
Pcube = [xg(:) yg(:), zg(:)];
% Extract the grid points located outside of the spherical region with
% force diameter 
Pcavitycube = Pcube(vecnorm(Pcube') > d,:);
% Create points on the sphere
[x1,y1,z1] = sphere(24);
Psphere = d*[x1(:) y1(:) z1(:)];
% Get rid of the upper half sphere
PsphereHalf = Psphere(z1(:) <= 0,:); 
PsphereHalf = unique(PsphereHalf,'rows');
% Combine the coordinates of the rectangular grid (without the points
% inside the sphere) and the surface coordinates of the unit sphere. 
Pcombined = [Pcavitycube;PsphereHalf];
% Create an alphaShape object representing the cube with the spherical
% cavity. 
shpCubeWithSphericalCavity = alphaShape(Pcombined(:,1), ...
                                        Pcombined(:,2), ...
                                        Pcombined(:,3));
% Increasing alpha to fill unconnected mesh in the middle
shpCubeWithSphericalCavity.Alpha = 20;
% cavityAlone = alphaShape(Pcavitycube(:,1), ...
%                         Pcavitycube(:,2), ...
%                         Pcavitycube(:,3));       
% sphereAlone = alphaShape(PsphereHalf(:,1), ...
%                         PsphereHalf(:,2), ...
%                         PsphereHalf(:,3));       
% figure
% plot(shpCubeWithSphericalCavity,'FaceAlpha',0.4)
% title('alphaShape: Cube with Half-Spherical Cavity')
% plot(cavityAlone,'FaceAlpha',0.4)
% plot(sphereAlone,'FaceAlpha',0.4)

% Recover the triangulation that defines the domain of the alphaShape object.
[tri,loc] = alphaTriangulation(shpCubeWithSphericalCavity);

% Create a PDE model.
modelCube = createpde;

% Create a geometry from the mesh and import the geometry and the mesh into the model.
[gCube,mshCube] = geometryFromMesh(modelCube,loc',tri');

% % Plot the resulting geometry.
% figure
% pdegplot(modelCube,'FaceAlpha',0.5,'CellLabels','on','FaceLabels','on')
% title('PDEModel: Cube with Half-Spherical Cavity')
% Identified that the half sphere is F7.
%% Solid Half-Sphere Nested in Cube
% Create tetrahedral elements to form a solid sphere by using the spherical
% shell and adding a new node at the center. First, obtain the spherical
% shell by extracting facets of the spherical boundary.  
sphereFacets = boundaryFacets(mshCube,'Face',7);
sphereNodes = findNodes(mshCube,'region','Face',7);

% Add a new node at the center.
newNodeID = size(mshCube.Nodes,2) + 1;

% Construct the tetrahedral elements by using each of the three nodes on the spherical boundary facets and the new node at the origin.
sphereTets =  [sphereFacets; newNodeID*ones(1,size(sphereFacets,2))];

% Create a model that combines the cube with the spherical cavity and a sphere.
model = createpde;

% Create a vector that maps all mshCube elements to cell 1, and all elements of the solid sphere to cell 2.
e2c = [ones(1,size(mshCube.Elements,2)), 2*ones(1,size(sphereTets,2))];

% Add a new node at the center [0;0;0] to the nodes of the cube with the cavity.
combinedNodes = [mshCube.Nodes,[0;0;0]];

% Combine the element connectivity matrices.
combinedElements = [mshCube.Elements,sphereTets];

% Create a two-cell geometry from the mesh.
[g,msh] = geometryFromMesh(model,combinedNodes,combinedElements,e2c);
 
% figure
% pdegplot(model,'FaceAlpha',0.5,'CellLabels','on','FaceLabels','on')
% title('Solid Sphere in Cube')


structModel.Geometry=g;

% //Specify material properties *******************************************
structuralProperties(structModel,'YoungsModulus',E,'PoissonsRatio',v);

% //Apply boundary constraints ********************************************
structuralBC(structModel,'Face',5,'Constraint','fixed'); %New face ID for the bottom.

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
structuralBoundaryLoad(structModel,'Face',8,'SurfaceTraction',hydrogelForce,'Vectorize','on'); %New face ID (F8)

%% //Generate mesh for model ***********************************************
generateMesh(structModel,'Hmax',40, 'Hmin',1, 'Hgrad', 1.2);
% pdeplot3D(structModel)
%% //Solving structural model **********************************************
structModelResults=solve(structModel);
%outputs
%   displacement = pix
%   stress = kg/(sec^2*pix) = Pa * ((10^9)/72)

%% //Visualizing results ***************************************************
figure(2)
pdeplot3D(structModel,'ColorMapData',structModelResults.Displacement.Magnitude, ...
    'Deformation',structModelResults.Displacement,'DeformationScaleFactor',1)
figure(4)
pdeplot3D(structModel,'ColorMapData',structModelResults.Displacement.uz, ...
    'Deformation',structModelResults.Displacement,'DeformationScaleFactor',1)

%% //Shifting bead locations to apply interpDisp at those locations ********
%2D
%Shifting bead locations because the outputs from simGaussianBeads are not
%centered around the origin which is required for our FEM model.
bead_xshifted=bead_x-(numPix_x/2)+0.5;     %bead_x from simGaussianBeads.m
bead_yshifted=bead_y-(numPix_x/2)+0.5;     %bead_y from simGaussianBeads.m
bead_zshifted=0*ones(nPoints,1);   %bead_z can be included in future
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
figure(3)
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
