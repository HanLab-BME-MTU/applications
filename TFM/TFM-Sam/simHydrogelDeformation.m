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
%       [bead_x,bead_y,bead_ux,bead_uy,Av]=simHydrogelDeformation(0,500,40,8000,1000,0.49,'/home/sehaarma/Documents/MATLAB/Bead Image Creation',false)
close all

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
if multiForce
    chr1=['multiforce_finemesh_',num2str(E),'E_',num2str(nPoints),'beads']; %chr13D='newBeads3D(1kpa,0.5k,1kE)';
else
    chr1=['singleforce_finemesh_',num2str(E),'E_',num2str(nPoints),'beads_(',num2str(fx),',',num2str(fy),')Pa'];
end
refstring=[chrRef chr1 '.tif'];
imgstring=[chrBead chr1 '.tif'];

% //Generate reference bead image *****************************************
meshPtsFwdSol = 2^9; %number of points in mesh
beadDiameter = 40; %20nm/bead
pixelSize = 72; %20 nm/pix
thickness = meshPtsFwdSol/4; %thickness of hydrogel for FEM purposes (pixels)
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
    
% //Generate force field **************************************************
if multiForce %multiple force applications
    [force_x,force_y] = generateMultiForce();
else %single force application
    forceType = 'groupForce';
    [x_mat_u, y_mat_u]=meshgrid(xmin:1:numPix_x,ymin:1:numPix_y);
    force_x = assumedForceAniso2D(1,x_mat_u,y_mat_u,numPix_x/2,numPix_y/2,fx,0,d,d,forceType);
    force_y = assumedForceAniso2D(2,x_mat_u,y_mat_u,numPix_x/2,numPix_y/2,0,fy,d,d,forceType);
end
force_z = [];


% //Generate model container **********************************************
structModel=createpde('structural','static-solid');

% //Define geometry *******************************************************
if ~multiForce
    
    d = 64; %diameter of fine mesh section
    
    %extrude method
    halfSide = 512/2;
    bound = [3; 4; -halfSide; halfSide; halfSide; -halfSide; halfSide; halfSide; -halfSide; -halfSide]; %rectangle of substrate size
    %create circles of decreasing radius in order 
    r = d/2;
    circ = [1; 0; 0; r; 0; 0; 0; 0; 0; 0]; %circle of radius r
    circ1 = [1; 0; 0; r/2; 0; 0; 0; 0; 0; 0]; %circle of radius r/2
    circ2 = [1; 0; 0; r/4; 0; 0; 0; 0; 0; 0]; %circle of radius r/4
    circ3 = [1; 0; 0; r/6; 0; 0; 0; 0; 0; 0]; %circle of radius r/6
    gd = [bound,circ,circ1,circ2,circ3];
    ns = char('bound','circ','circ1','circ2','circ3');
    ns = ns';
    sf = 'bound + circ + circ1 + circ2 + circ3';    
    [dl, ~] = decsg(gd,sf,ns);
        
    pdem = createpde;
    g = geometryFromEdges(pdem,dl);
    
    %convert analytical model geometry to discrete geometry
    facets = facetAnalyticGeometry(pdem,g,0);
    gm = analyticToDiscrete(facets);
    pdem.Geometry = gm; %reassign discrete geometry to model

    %extrude 2D geometry into 3D of defined thickness
    g = extrude(gm,thickness);
    
elseif multiForce

    imScale = 10;
    %initialize ellipse locations and angles
    fpos = [10*imScale,10*imScale; 12*imScale,11*imScale; 15*imScale,10*imScale; 11*imScale,14*imScale; % 1st cluster
            15*imScale,26*imScale; 17*imScale,31*imScale; 28*imScale,27*imScale; 21*imScale,20*imScale; % 2nd cluster
            27*imScale,39*imScale; 31*imScale,42*imScale; 37*imScale,43*imScale; 31*imScale,36*imScale; 38*imScale,40*imScale; 40*imScale,33*imScale; % 3rd cluster
            33*imScale,9*imScale; 39*imScale,18*imScale; 43*imScale,24*imScale; 38*imScale,12*imScale; 42*imScale,17*imScale]; % 4th cluster

    deg =  [4.76364; 5.07961; 4.95326; 4.83658;
            90-19.0577; -20.2657; -(90+18.2325); -1.4321;
            -(180+50.1944); -(57.2648); 90-74.4759; -(180+56.3099); -(90+70.71); -(90+70.71);
            -48.0128; -53.9726; -(90-48.8141); -(90-55.008); -75.9638];

    %conduct first decomposition outside for loop to include bounding box
    bound = [2; 4; 0; 512; 512; 0; 512; 512; 0; 0;];
    circ = [4 fpos(1,1) -fpos(1,2)+500 8 16 deg(1) 0 0 0 0]';
    gd = [bound,circ];
    ns = char('bound','circ');
    ns = ns';
    sf = 'bound + circ';
    [dl, ~] = decsg(gd,sf,ns);

    %define ellipse locations
    x = fpos(:,1); y = -fpos(:,2) + 500;
    %define ellipse axes lengths
    s1 = 7.5*ones(length(fpos),1); s2 = 15*ones(length(fpos),1);
    %define ellipses and conduct decomposition
    for i = 2:length(x)
        gd = [4, x(i), y(i), s1(i), s2(i), deg(i), 0, 0, 0, 0]';
        [dlt, ~] = decsg(gd);
        dl = [dl dlt];
    end

    %augment dl matrix to ensure correct face IDs
    dl(7,:) = ones(length(dl),1); %set outside face to one for all ellipses
    m = 9; n = 12; %initialize matrix locations at 1st edge of 2nd ellipse
    for j = 3:20
        dl(6,m:n) = j * ones(4,1); %incrementally set inside face ID for each ellipse
        m = m + 4; n = n + 4; %update matrix location for next ellipse
    end

    %assign model geometry
    pdem = createpde;
    g = geometryFromEdges(pdem,dl);

    %convert analytical model geometry to discrete geometry
    facets = facetAnalyticGeometry(pdem,g,0);
    gm = analyticToDiscrete(facets);
    pdem.Geometry = gm; %reassign discrete geometry to model

    %extrude 2D geometry into 3D of defined thickness
    g = extrude(gm,thickness);
    
end

structModel.Geometry=g;

pdegplot(structModel)

% //Specify material properties *******************************************
structuralProperties(structModel,'YoungsModulus',E,'PoissonsRatio',v);

% //Apply boundary constraints *******************************************
if multiForce
    structuralBC(structModel,'Face',1,'Constraint','fixed'); %New face ID for the bottom.
elseif ~multiForce
    structuralBC(structModel,'Face',1:5,'Constraint','fixed'); %New face ID for the bottom.
end

% //Apply force distribution on face 2 ************************************
if ~multiForce
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
elseif multiForce
    if isempty(force_z)
        force_z=zeros(numPix_x); %can be replaced by real data if a normal force is desired
    end
    [xLoc,yLoc]=ndgrid(1:numPix_x,1:numPix_y);
    forceInterpX=griddedInterpolant(xLoc,yLoc,force_x);
    forceInterpY=griddedInterpolant(xLoc,yLoc,force_y);
    forceInterpZ=griddedInterpolant(xLoc,yLoc,force_z);
    %define function handle
    hydrogelForce = @(location,state)[forceInterpX(location.x,location.y); ...
                                      forceInterpY(location.x,location.y); ...
                                      forceInterpZ(location.x,location.y);];
end

% //Generate mesh for model ***********************************************
generateMesh(structModel,'Hmax',40, 'Hmin',1, 'Hgrad', 1.2);
% pdeplot3D(structModel)

%pass function handle and define BCs
if multiForce
    structuralBoundaryLoad(structModel,'Face',21,'SurfaceTraction',hydrogelForce,'Vectorize','on'); %New face ID (F8)
elseif ~multiForce
    structuralBoundaryLoad(structModel,'Face',6:10,'SurfaceTraction',hydrogelForce,'Vectorize','on'); %New face ID (F8)
end

% //Solving structural model **********************************************
structModelResults=solve(structModel);

% //Visualizing results ***************************************************
figure(2)
pdeplot3D(structModel,'ColorMapData',structModelResults.Displacement.Magnitude, ...
    'Deformation',structModelResults.Displacement,'DeformationScaleFactor',1)
figure(4)
pdeplot3D(structModel,'ColorMapData',structModelResults.Displacement.Magnitude,'Mesh','on')

% //Shifting bead locations to apply interpDisp at those locations ********
%2D
%Shifting bead locations because the outputs from simGaussianBeads are not
%centered around the origin which is required for our FEM model.
if ~multiForce
    % top of the geometry is located at z = thickness
    bead_zshifted=thickness*ones(nPoints,1);   %interpolation location in z direction
    bead_xshifted=bead_x-(numPix_x/2)+0.5;     %bead_x from simGaussianBeads.m
    bead_yshifted=bead_y-(numPix_x/2)+0.5;     %bead_y from simGaussianBeads.m
elseif multiForce
    % top of the geometry is located at z = thickness
    bead_zshifted=thickness*ones(nPoints,1);  
    bead_xshifted=bead_x;     %bead_x from simGaussianBeads.m
    bead_yshifted=bead_y;     %bead_y from simGaussianBeads.m
end
    
%3D
% bead_x3Dshifted=beadcenters3D(:,1)-(meshPtsFwdSol/2)+0.5; %beadcenters3D from simGaussianSpots3D.m
% bead_y3Dshifted=beadcenters3D(:,2)-(meshPtsFwdSol/2)+0.5; %beadcenters3D from simGaussianSpots3D.m
% bead_z3Dshifted=beadcenters3D(:,3)-0.5;   %beadcenters3D from simGaussianSpots3D.m

% //Interpolating displacements at bead locations using interpDisp ********
%2D
interpDisp=interpolateDisplacement(structModelResults,bead_xshifted,bead_yshifted,bead_zshifted);
%3D
% interpDisp3D=interpolateDisplacement(structModelResults,bead_x3Dshifted,bead_y3Dshifted,bead_z3Dshifted);

% //Reshaping to match bead location vectors ******************************
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

% //Plotting interpolated displacement ************************************
if ~multiForce
    figure, grid on;
    q2=quiver(bead_xshifted,bead_yshifted,bead_ux,bead_uy);
    xlim([-(meshPtsFwdSol/2) (meshPtsFwdSol/2)]); ylim([-(meshPtsFwdSol/2) (meshPtsFwdSol/2)]);
    xlabel('X'); ylabel('Y'); title('Ground Truth Displacement Field - 1 kPa');
elseif multiForce
    figure, grid on;
    q2=quiver(bead_xshifted,bead_yshifted,bead_ux,bead_uy);
    xlim([0 numPix_x]); ylim([0 numPix_y]);
    xlabel('X'); ylabel('Y'); title('Ground Truth Displacement Field - 1 kPa');
end
% //Colormapping *********************************************************
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

% //Generating Deformed Bead Image ****************************************
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

% //Saving outputs for comparison after TFM processing ********************
disp(strcat('Removed ',num2str(nanElements),'NaN elements.'))
save(strcat('outputs_',chr1,'.mat'),'bead_ux','bead_uy','force_x','force_y'); %chr1 defined initially by user, matches image file naming
