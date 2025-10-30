function [bead_x,bead_y,bead_ux,bead_uy,Av]=simHydrogelDeformationAbaqus(fx,fy,d,nPoints,E,v,dataPath,multiForce)
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
%       [bead_x,bead_y,bead_ux,bead_uy,Av]=simHydrogelDeformationAbaqus(0,500,40,8000,1000,0.49,'/home/sehaarma/matlab/abaqus bead image testing',false)
close all

%% Data path configuration
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

%% Generate reference bead image
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

%Save bead coords and displacements to movieData structure
displField(1).pos = [bead_x, bead_y];
displField(1).vec = [zeros(length(bead_x),1), zeros(length(bead_y),1)];

%Save bead coords and forces to movieData structure
forceField(1).pos = [bead_x, bead_y];
forceField(1).vec = [zeros(length(bead_x),1), zeros(length(bead_y),1)];
    
%noise addition
refimg = refimg+0.05*rand(numPix_y,numPix_x)*max(refimg(:));
%add 3D noise addition in future

%writing reference image
imwrite(uint16(refimg*2^16/max(max(refimg))),[refPath filesep refstring]);
%add saving 3D refimg
    
%% Generate force field
[x_mat_u, y_mat_u]=meshgrid(xmin:1:numPix_x,ymin:1:numPix_y);
if multiForce %multiple force applications
    [force_x,force_y] = generateMultiForce();
else %single force application
    forceType = 'groupForce';
    force_x = assumedForceAniso2D(1,x_mat_u,y_mat_u,numPix_x/2,numPix_y/2,fx,0,d,d,forceType);
    force_y = assumedForceAniso2D(2,x_mat_u,y_mat_u,numPix_x/2,numPix_y/2,0,fy,d,d,forceType);
end
force_z = zeros(size(force_x));

%% Fit ellipses to each traction island
% This is done to obtain the locations, dimensions, and orientations for 
% the mesh refinement and UTRACLOAD subroutine
% ** MIGHT NEED TO DILATE ELLIPSES TO ENSURE ALL FORCE AREA IS CAPTURED **
tractionMask = abs(force_x) + abs(force_y) + abs(force_z); %combine x,y,z tractions to create mask of all nonzero elements
tractionMask(tractionMask ~= 0) = 1; %convert nonzero elements to 1
cc = bwconncomp(tractionMask);
s = regionprops(cc,{'Centroid','Orientation','MajorAxisLength','MinorAxisLength'});
figure, imshow(tractionMask,[]), hold on
theta = linspace(0,2*pi);
orientationList = zeros(size(s,1),2);
for i = 1:length(s)
    col = (s(i).MajorAxisLength/2)*cos(theta);
    row = (s(i).MinorAxisLength/2)*sin(theta);
    M = makehgtform('translate',[s(i).Centroid, 0],'zrotate',deg2rad(-1*s(i).Orientation));
    [orientationList(i,1),orientationList(i,2)] = pol2cart(deg2rad(-1*s(i).Orientation),1);
    D = M*[col;row;zeros(1,numel(row));ones(1,numel(row))];
    plot(D(1,:),D(2,:),'r','LineWidth',2); hold on
    s(i).Orientation = -1 * s(i).Orientation; %convert orientation to xy coordinate frame
end
hold off

B = bwboundaries(tractionMask); %get boundaries of traction islands
tractionPolygons = cell(length(B)*2,1); %initialize new cell array
tractionPolygons(1:2:end,:) = B; %separate each traction island boundary point set
tractionPolygons(2:2:end,:) = {[NaN,NaN]}; %add NaN in between point sets
tractionPolygons = vertcat(tractionPolygons{:}); %collapse cell array into array of islands separated by NaNs

%% Define geometry
%use a matlab pde toolbox model to easily create and mesh geometry
abaqusModel=createpde('structural','static-solid'); 
if ~multiForce
    halfSide = 512/2;
    bound = [3; 4; -halfSide; halfSide; halfSide; -halfSide; halfSide; halfSide; -halfSide; -halfSide]; %rectangle of substrate size
    %create circles of decreasing radius in order 
    d = 64; %diameter of fine mesh section
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
    abaqusModel.Geometry = g; %assign 3D geometry to model

elseif multiForce
    %conduct first decomposition outside for loop to include bounding box
    halfSide = 512/2;
    bound = [3; 4; -halfSide; halfSide; halfSide; -halfSide; halfSide; halfSide; -halfSide; -halfSide]; %rectangle of substrate size
    elli = [4; s(1).Centroid(1)-halfSide; s(1).Centroid(2)-halfSide; s(1).MajorAxisLength/2; s(1).MinorAxisLength/2; deg2rad(s(1).Orientation); 0; 0; 0; 0];
    gd = [bound,elli];
    ns = char('bound','elli');
    ns = ns';
    sf = 'bound + elli';
    [dl, ~] = decsg(gd,sf,ns);

    %define ellipses and conduct decomposition
    for i = 2:length(s)
        elli = [4; s(i).Centroid(1)-halfSide; s(i).Centroid(2)-halfSide; s(i).MajorAxisLength/2; s(i).MinorAxisLength/2; deg2rad(s(i).Orientation); 0; 0; 0; 0];
        [dlt, ~] = decsg(elli);
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
    figure,pdegplot(pdem)
    %convert analytical model geometry to discrete geometry
    facets = facetAnalyticGeometry(pdem,g,0);
    gm = analyticToDiscrete(facets);
    pdem.Geometry = gm; %reassign discrete geometry to model

    %extrude 2D geometry into 3D of defined thickness
    g = extrude(gm,thickness);
    abaqusModel.Geometry = g; %assign 3D geometry to model
end

%% Mesh the 3D geometry
%Hmax = maximum element size
%Hmin = minimum element size
%Hgrad = mesh size gradient (how quickly size changes)
%generateMesh(abaqusModel,'Hmax',30, 'Hmin',1, 'Hgrad', 1.1); %very fine mesh
generateMesh(abaqusModel,'Hmax',50, 'Hmin',5, 'Hgrad', 2); %coarse mesh
pdemesh(abaqusModel)

%% Obtain node IDs and locations needed for BC applications
topFaceNodeIDs = findNodes(abaqusModel.Mesh,"box",[-halfSide,halfSide],[-halfSide,halfSide],[thickness-0.0001,thickness]); %get nodes on top face
topFaceNodes = abaqusModel.Mesh.Nodes(:,topFaceNodeIDs); %get top node locations
bottomFaceNodeIDs = findNodes(abaqusModel.Mesh,"box",[-halfSide,halfSide],[-halfSide,halfSide],[0, 0.0001]); %get nodes on bottom face
bottomFaceNodes = abaqusModel.Mesh.Nodes(:,bottomFaceNodeIDs); %get bottom node locations
nonBCNodeIDs = 1:length(abaqusModel.Mesh.Nodes); %define node set IDs for all other nodes
nonBCNodeIDs([topFaceNodeIDs,bottomFaceNodeIDs]) = []; %clear node IDs already contained in top and bottom face node sets
nonBCNodes = abaqusModel.Mesh.Nodes(:,nonBCNodeIDs); %associate non boudnary node locations with node set IDs

%% Use fitted ellipses to extract all nodes in the traction area
idx = inpolygon(topFaceNodes(1,:),topFaceNodes(2,:),tractionPolygons(:,2)-halfSide,tractionPolygons(:,1)-halfSide);
nodeIDsInFootprint = topFaceNodeIDs(idx);
%Plot footprint nodes alongside mesh and traction mask to verify
figure,plot(topFaceNodes(1,logical(idx)),topFaceNodes(2,logical(idx)),'.')
figure,imshow(tractionMask,[]); axis xy
figure,pdemesh(abaqusModel); axis xy

%% Obtain element IDs and locations for BC applications
topFaceElementIDs = findElements(abaqusModel.Mesh,'attached',topFaceNodeIDs); %get element IDs for all elements with a node on the top face
bottomFaceElementIDs = findElements(abaqusModel.Mesh,'attached',bottomFaceNodeIDs); %get element IDs for all elements with a node on the bottom face
allElementIDs = 1:length(abaqusModel.Mesh.Elements); %extract all element IDs
allElements = abaqusModel.Mesh.Elements; %get all element definitions using element IDs
topFaceElements = abaqusModel.Mesh.Elements(:,topFaceElementIDs); %extract elements using element IDs from the top face
bottomFaceElements = abaqusModel.Mesh.Elements(:,bottomFaceElementIDs); %extract elements using element IDs from the bottom face
nonBCElementIDs = 1:length(abaqusModel.Mesh.Elements); %define element set IDs for all other nodes
nonBCElementIDs([topFaceElementIDs,bottomFaceElementIDs]) = []; %exclude elements already included in top face and bottom face sets
nonBCElements = abaqusModel.Mesh.Elements(:,nonBCElementIDs); %associate non boundary elements with element set IDs

elemIDsInFootprint = findElements(abaqusModel.Mesh,"attached",nodeIDsInFootprint); %get IDs of all elements with at least 1 corner node on the top face
elemsInFootprint = abaqusModel.Mesh.Elements(:,elemIDsInFootprint); %get node IDs for all elements found on prev line
uniqueElemsInFootprint = unique(elemsInFootprint)'; %exclude nodes shared between elements and compress array
nodesInFootprint = abaqusModel.Mesh.Nodes(:,uniqueElemsInFootprint); %replace with line below if deep node search is not desired
%nodesInFootprint = modelMesh.Nodes(:,nodeIDsInFootprint);
nodesInFootprint2D = nodesInFootprint(1:2,:); %extract x and y dimensions

%% Interpolate force maps onto footprint node locations
tractionMag = sqrt(force_x.^2 + force_y.^2 + force_z.^2);
% forceNorm_x = force_x ./ tractionMag;
% forceNorm_y = force_y ./ tractionMag;
% forceNorm_z = force_z ./ tractionMag;

[Xgrid,Ygrid] = ndgrid(linspace(-meshPtsFwdSol/2,meshPtsFwdSol/2,512),linspace(-meshPtsFwdSol/2,meshPtsFwdSol/2,512));
tractionInterp = griddedInterpolant(Xgrid,Ygrid,tractionMag);
tractionAtNodesInFootprint = tractionInterp(nodesInFootprint2D(1,:),nodesInFootprint2D(2,:));

%% Set orientation information for each footprint node
% The footprint node inherits the orientation of the force area it is from
orientations = zeros(size(nodesInFootprint2D));
for i = 1:size(orientationList,1)
    curPoly = B{i} - halfSide;
    idx = inpolygon(nodesInFootprint2D(1,:),nodesInFootprint2D(2,:),curPoly(:,2),curPoly(:,1));
    for j = 1:length(idx)
        if idx(j)
            orientations(:,j) = orientationList(i,:);
        end
    end
end
orientations = [orientations; zeros(1,size(orientations,2))];
magOrients = sqrt(orientations(1,:).^2 + orientations(2,:).^2);
figure, plot3(nodesInFootprint2D(1,:),nodesInFootprint2D(2,:),magOrients,'.')

%Some footprint nodes are missed so we search nearest non-skipped nodes to
%apply the same orientation as its nearest neighbor has
missedPts = magOrients==0; %logical array of all pts missed during fitting
queryPts = nodesInFootprint2D(:,missedPts)'; %trim down to only those missed points
donePts = nodesInFootprint2D';
donePts(magOrients==0,:) = -512; %master list of all points with missed points sent outside substrate to be ignored during neighbor search
k = dsearchn(donePts,queryPts); %search for closest neighbor to missed points
%Plot search results
figure, plot(donePts(:,1),donePts(:,2),'ko'); hold on
plot(queryPts(:,1),queryPts(:,2),'*g'); hold on
plot(donePts(k,1),donePts(k,2),'*r'); hold off
xlim([-256 256]); ylim([-256 256]);
orientations(:,missedPts) = orientations(:,k); %set missed orientations to nearest neighbor's orientation
%Plot magnitude of orientations to verify all unit vectors add to one and
%none have been missed
magOrients = sqrt(orientations(1,:).^2 + orientations(2,:).^2);
figure, plot3(nodesInFootprint2D(1,:),nodesInFootprint2D(2,:),magOrients,'.')

%Plot the node locations used in traction definition for verification of
%deep node search
figure, plot3(nodesInFootprint(1,:),nodesInFootprint(2,:),nodesInFootprint(3,:),'b.'), grid on
title('Node Locations')

%% Create matrix of values for user subroutine creation
% This subroutineVals matrix contains the traction and orientation
% information for each element with nodes that fall within the traction areas 
% Format: [Element#, IntegrationPoint#, Traction Value, X-dir, Y-dir, Z-dir]
k = 1; %initialize loop variable
subroutineVals = zeros(6,size(topFaceElements,1)*size(topFaceElements,2)); %overinitialize to maximum possible size
for i = 1:length(uniqueElemsInFootprint) %loop over all nodes within force area
    %[row,col] = find(topFaceElements == nodeIDsInFootprint(i)); %search top face elements for current node
    [row,col] = find(topFaceElements == uniqueElemsInFootprint(i)); %search top face elements for current node
    for j = 1:length(row) %loop over all found elements
        % assemble subroutine information matrix as [elementID, intergrationPointID, tractionVal, orientation];
        subroutineVals(:,k) = [topFaceElementIDs(col(j));row(j);tractionAtNodesInFootprint(i);orientations(1,i);orientations(2,i);orientations(3,i)]; 
        k = k + 1; %iterate loop variable
    end
end
subroutineVals = subroutineVals(:,1:find(subroutineVals(1,:),1,'last')); %trim off any trailing zeros from overinitialization
subroutineVals = sortrows(subroutineVals',[1,2])'; %sort by element number and integration point number within element numbers

%% Write user-subroutine UTRACLOAD to define tractions and orientations
%if UTRACLOAD subroutine file exists already we must delete it 
if isfile([dataPath,'/abaqus_input.f'])
    delete([dataPath,'/abaqus_input.f'])
end
%create the subroutine file with write permissions
fileIDsub = fopen([dataPath,'/abaqus_input.f'],'w'); %name of subroutine must match name of .inp file

%write subroutine header
fprintf(fileIDsub,'      SUBROUTINE UTRACLOAD(ALPHA,TUSER,KSTEP,KINC,TIME,NOEL,NPT,COORDS,DIRCOS,JLTYP,SNAME)\n');
fprintf(fileIDsub,'      INCLUDE ''ABA_PARAM.INC''\n');
fprintf(fileIDsub,'      DIMENSION TUSER(3),TIME(2),COORDS(3), DIRCOS(3,3)\n');
fprintf(fileIDsub,'      CHARACTER*80 SNAME\n');

%writing subroutine logic
g = 1;
elemList = unique(subroutineVals(1,:));
fprintf(fileIDsub,'      if(noel .eq. %d) then\n',subroutineVals(1,1));
for i = 1:length(elemList)
    nextElem = find(subroutineVals(1,:) == elemList(i),1,'last');
    elemLength = nextElem - g + 1;
    if g ~= 1
        fprintf(fileIDsub,'      else if(noel .eq. %d) then\n',subroutineVals(1,g));
    end
    for j = 1:elemLength
        if j == 1
            fprintf(fileIDsub,'            if(npt .eq. %d) then\n',subroutineVals(2,g));
            fprintf(fileIDsub,'                  alpha=%G\n',subroutineVals(3,g));
            fprintf(fileIDsub,'                  TUSER(1)=%0.4f\n',subroutineVals(4,g));
            fprintf(fileIDsub,'                  TUSER(2)=%0.4f\n',subroutineVals(5,g));
            fprintf(fileIDsub,'                  TUSER(3)=%0.4f\n',subroutineVals(6,g));
            g = g + 1;
        else
            fprintf(fileIDsub,'            else if(npt .eq. %d) then\n',subroutineVals(2,g));
            fprintf(fileIDsub,'                  alpha=%G\n',subroutineVals(3,g));
            fprintf(fileIDsub,'                  TUSER(1)=%0.4f\n',subroutineVals(4,g));
            fprintf(fileIDsub,'                  TUSER(2)=%0.4f\n',subroutineVals(5,g));
            fprintf(fileIDsub,'                  TUSER(3)=%0.4f\n',subroutineVals(6,g));
            g = g + 1;
        end
        if j == elemLength
            fprintf(fileIDsub,'            endif\n');
        end
    end
end
fprintf(fileIDsub,'      endif\n');
fprintf(fileIDsub,'\n      RETURN\n      END');
fclose(fileIDsub); %close the subroutine now that we are done writing

%% Begin writing to input file
%if the input file exists delete it
if isfile([dataPath,'/abaqus_input.inp'])
    delete([dataPath,'/abaqus_input.inp'])
end
%create the input file with write permissions
fileID = fopen([dataPath,'/abaqus_input.inp'],'w');

%% Write *HEADING section to input file
fprintf(fileID,['*HEADING\n3-dimensional rectangular substrate tangential traction\n' ...
    'SI Units\n1-axis horizontal, 2-axis vertical, 3-axis OOP normal\n']);

%% Write *PREPRINT section to input file
fprintf(fileID,'*PREPRINT, ECHO=YES, MODEL=YES, HISTORY=YES\n');

%% Write *NODE sections to  input file
%this section defines the nodes and their locations
%top face node set
fprintf(fileID,'*NODE, NSET=TOPNODES\n');
for ii = 1:length(topFaceNodeIDs)
    fprintf(fileID,'%7d,%13d,%13d,%13d\n',topFaceNodeIDs(ii),topFaceNodes(:,ii));
end
%bottom face node set
fprintf(fileID,'*NODE, NSET=BOTTOMNODES\n');
for jj = 1:length(bottomFaceNodeIDs)
    fprintf(fileID,'%7d,%13d,%13d,%13d\n',bottomFaceNodeIDs(jj),bottomFaceNodes(:,jj));
end
%rest of nodes
fprintf(fileID,'*NODE, NSET=BULK\n');
for kk = 1:length(nonBCNodeIDs)
    fprintf(fileID,'%7d,%13d,%13d,%13d\n',nonBCNodeIDs(kk),nonBCNodes(:,kk));
end

%% Write *ELEMENT sections to input file
%this section defines the element type and the nodes that make up the elements
%10-node quadratic tet element type = C3D10
%10-node quadratic tet element with hybrid stress calculation = C3D10HS
%10-node modified tetrahedron, with hourglass control, hybrid with linear pressure = C3D10MH
elementType = 'C3D10HS';
%elementType = 'C3D10MH';
%top face elements
fprintf(fileID,['*ELEMENT, TYPE=',elementType,', ELSET=TOPELEMS\n']);
for te = 1:length(topFaceElementIDs)
    fprintf(fileID,'%7d,%7d,%7d,%7d,%7d,%7d,%7d,%7d,%7d,%7d,%7d\n',topFaceElementIDs(te),topFaceElements(:,te));
end
%bottom face element set
fprintf(fileID,['*ELEMENT, TYPE=',elementType,', ELSET=BOTTOMELEMS\n']);
for be = 1:length(bottomFaceElementIDs)
    fprintf(fileID,'%7d,%7d,%7d,%7d,%7d,%7d,%7d,%7d,%7d,%7d,%7d\n',bottomFaceElementIDs(be),bottomFaceElements(:,be));
end
%bulk elements
fprintf(fileID,['*ELEMENT, TYPE=',elementType,', ELSET=BULK\n']);
for ne = 1:length(nonBCElementIDs)
    fprintf(fileID,'%7d,%7d,%7d,%7d,%7d,%7d,%7d,%7d,%7d,%7d,%7d\n',nonBCElementIDs(ne),nonBCElements(:,ne));
end
%combine element sets into an element set containing all elements
fprintf(fileID,['*ELSET, ELSET=ALLELEMS\n' ...
    'TOPELEMS, BOTTOMELEMS, BULK\n']);

%% Write *SURFACE section to input file
%this section defines the surface names and elements it contains for
%application of traction
fprintf(fileID,['*SURFACE, NAME=TOPSURF, TYPE=ELEMENT\n' ...
    'TOPELEMS,\n']);

%% Write *SOLID SECTION to input file 
%this section defines the solid elements associated with a material name
%onto which, material properties are associated
fprintf(fileID,'*SOLID SECTION, ELSET=ALLELEMS, MATERIAL=SUBSTRATE\n');

%% Write *MATERIAL section to input file
%this section defines the material properties of the material name defined
%in the previous step
E = 1000;
v = 0.49;
matType = 'ogdenL3'; %linear, neohookean, ogdenPVA, ogdenL3, ogdenThetaFive, neohookeanThetaFive
switch matType
    case 'linear'
        %write material section to file
        fprintf(fileID,['*MATERIAL, NAME=SUBSTRATE\n' ...
           '*ELASTIC\n' ...
           '%7d, %7d\n'],E,v); %young's modulus and poisson's ratio respectively
    case 'neohookean'
        %given Young's Modulus of 1000 and Poisson's Ratio of 0.48, the bulk
        %modulus and shear modulus can be calculated. These two quantities can then
        %be used to calculate C1 and D1 for a neo hookean material.
        G = E / (2 * (1 + v)); %shear modulus from E and v
        K = E / (3 * (1 - 2*v)); %bulk modulus from E and V (aka 2nd Lame param)
        C10 = G / 2; %1st Lame param from bulk modulus
        D1 = ((K - (2/3 * G)) / 2); %
        %D1_equiv = ((E*v)/((1+v)*(1-2*v))) / 2;
        %D1 = 10;
        %C10 = 1.5*10^6;
        PR0 = -D1 / (2 * C10);
        disp(['Initial Poisson''s Ratio is ',num2str(PR0)])
        fprintf(fileID,['*MATERIAL, NAME=SUBSTRATE\n' ...
             '*HYPERELASTIC, NEO HOOKE\n' ...
             '%7d, %7d\n'],C10,D1);
        %important abaqus documentation regarding nonlinearity
        %https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/gss/default.htm?startat=ch07s03.html
    case 'ogdenPVA'
        mu1 = -2.434 * 10^6;
        al1 = 0.422;
        mu2 = 1.413 * 10^6;
        al2 = 1.237;
        mu3 = 1.042 * 10^6;
        al3 = -0.313;
        D1 = 0;
        D2 = 0;
        D3 = 0;
        fprintf(fileID,['*MATERIAL, NAME=SUBSTRATE\n' ...
             '*HYPERELASTIC, OGDEN\n' ...
             '%7d, %7d, %7d, %7d, %7d, %7d, %7d, %7d, %7d\n'],mu1,al1,mu2,al2,mu3,al3,D1,D2,D3);
    case 'ogdenL3'
        mu1 = 6.18 * 10^3;
        al1 = -4.44;
        mu2 = 12.15 * 10^3;
        al2 = 5.15;
        mu3 = -13.18 * 10^3;
        al3 = -9.81;
        D1 = 0;
        D2 = 0;
        D3 = 0;
        fprintf(fileID,['*MATERIAL, NAME=SUBSTRATE\n' ...
             '*HYPERELASTIC, OGDEN\n' ...
             '%7d, %7d, %7d, %7d, %7d, %7d, %7d, %7d, %7d\n'],mu1,al1,mu2,al2,mu3,al3,D1,D2,D3); 
    case 'ogdenThetaFive'
        mu1 = 1.6 * 10^(-3);
        al1 = 2.349;
        mu2 = 1.5 * 10^(-3);
        al2 = 2.69;
        mu3 = 1.5 * 10^(-3);
        al3 = 2.69;
        D1 = 11.52;
        D2 = 0;
        D3 = 0;
        fprintf(fileID,['*MATERIAL, NAME=SUBSTRATE\n' ...
             '*HYPERELASTIC, OGDEN\n' ...
             '%7d, %7d, %7d, %7d, %7d, %7d, %7d, %7d, %7d\n'],mu1,al1,mu2,al2,mu3,al3,D1,D2,D3);
    case 'neohookeanThetaFive'
        C10 = 2.23*10^-3;
        D1 = 11.96;
        fprintf(fileID,['*MATERIAL, NAME=SUBSTRATE\n' ...
             '*HYPERELASTIC, NEO HOOKE\n' ...
             '%7d, %7d\n'],C10,D1);
end

%% Write *STEP section to input file
%this section defines the incrementation scheme of the FEA algorithm as
%well as some options to improve convergence in nonlinear models
%*STEP,UNSYMM=YES can be used to improve convergence if problems arise
fprintf(fileID,'*STEP, AMPLITUDE=RAMP, INC=200, NLGEOM=YES\n');
fprintf(fileID,'*STATIC, STABILIZE, FACTOR=2\n'); %STABILIZE is used in conjunction with NLGEOM to reduce nonlinear instability (link below)
%http://130.149.89.49:2080/v6.13/books/key/default.htm?startat=ch18abk31.html#usb-kws-hstatic
fprintf(fileID,'0.005, 1.0, 1e-15, 0.1\n'); %Initial, Step Time, Min, Max (1e-07 might further improve convergence if needed) https://www.eng-tips.com/viewthread.cfm?qid=188499
fprintf(fileID,'*CONTROLS, PARAMETERS=LINE SEARCH\n');
fprintf(fileID,'10, 1.0, 0.0001, 0.25, 0.01\n');

%% Write *BOUNDARY section to input file
%this section defines the boundary conditions, in this case fixing all
%nodes belonging to the bottom face
fprintf(fileID,'*BOUNDARY\n');
fprintf(fileID,'BOTTOMNODES, ENCASTRE\n');

%% Write *DSLOAD section to input file
%this section utilizes nonuniform traction defined by the .f UTRACLOAD 
%subroutine via the DSLOAD and TRVECNU names
fprintf(fileID,['*DSLOAD\n' ...
    'TOPSURF, TRVECNU, 1., 1., 0., 0.\n']);

%% Write *NODE PRINT section to input file
%this section defines the parameters we want printed to the history file
%and which nodes we want to write those parameters from
%FREQUENCY value is chosen to be very big such that only the final results 
%are printed to the data file
fprintf(fileID,['*NODE PRINT, NSET=TOPNODES, FREQUENCY=10000\n' ...
    'U1,U2,U3\n']);

%% Finish writing and close the .inp input file
%end step
fprintf(fileID,'*END STEP');
%close file
fclose(fileID);

%% Datacheck (optional) to verify input file formatting
%run abaqus datacheck to verify input file
%system('abaqus job=test inp=abaqus_input.inp datacheck interactive');
%here we will check the test.dat file for ***ERROR, if no ***ERROR text is found then the model is ready to run

%% Run ABAQUS from command line
%this section specifies the locations of the input and subroutine files as
%well as the number of cpu cores to use and the verbosity of the .msg files
% ** ADD EXACT PATH TO INPUT AND SUBROUTINE FILES ** (using dataPath)
system(['/opt/abaqus/v2023/SIMULIA/Commands/abaqus job=abaqus_test inp=abaqus_input.inp ' ... 
    ,'user=abaqus_input.f verbose=3 cpus=3 interactive ask_delete=OFF'])

%% Read nodal displacement results from .dat file
%this section loads the .dat file and loads the nodal displacements from
%the node output section
nodeData = zeros(length(topFaceNodeIDs),4);
nodeID = zeros(length(topFaceNodeIDs),1);
dispOutput = zeros(length(topFaceNodeIDs),3);

file = [dataPath,'/abaqus_test.dat'];
fid = fopen(file);
i=0; j=0;
while ( ~feof(fid) )
    tline = fgetl(fid);
    i = i+1;
    if (regexpi(tline, 'N O D E   O U T P U T')>0)
        tline = fgetl(fid);
        i = i+1;
        while(isempty(str2num(tline))) %#ok<ST2NM>
            tline = fgetl(fid);
            i = i+1;
        end
        while(~isempty(str2num(tline))) %#ok<ST2NM>
            j=j+1;
            nodeData(j,:) = sscanf(tline, '%d %e %e %e', [1,4]);
            nodeID(j)=nodeData(j,1);
            tline = fgetl(fid);
            i = i+1;
            dispOutput(j,:) = nodeData(j,2:4);
        end
        break
    end
end
fclose(fid);
fclose all; %fclose all here to make sure the input, subroutine, and data files are closed

%% Plotting displacement results
[x,y] = meshgrid(-halfSide:1:halfSide,-halfSide:1:halfSide);
dispMatX = griddata(topFaceNodes(1,:),topFaceNodes(2,:),dispOutput(:,1),x,y);
dispMatY = griddata(topFaceNodes(1,:),topFaceNodes(2,:),dispOutput(:,2),x,y);
figure, imagesc(dispMatX);
figure, imagesc(dispMatY)

dispMatZ = griddata(topFaceNodes(1,:),topFaceNodes(2,:),dispOutput(:,3),x,y); %dispOutput(:,3) takes the z-direction of displacement
figure, title('Z-component Displacement')
surf(x,y,dispMatZ)

%% Shift bead coordinates to align with FEM substrate coordinate system
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
    bead_xshifted=bead_x-(numPix_x/2)-0.5;     %bead_x from simGaussianBeads.m
    bead_yshifted=bead_y-(numPix_x/2)-0.5;     %bead_y from simGaussianBeads.m
end
%3D
% bead_x3Dshifted=beadcenters3D(:,1)-(meshPtsFwdSol/2)+0.5; %beadcenters3D from simGaussianSpots3D.m
% bead_y3Dshifted=beadcenters3D(:,2)-(meshPtsFwdSol/2)+0.5; %beadcenters3D from simGaussianSpots3D.m
% bead_z3Dshifted=beadcenters3D(:,3)-0.5;   %beadcenters3D from simGaussianSpots3D.m

%% Interpolate nodal displacements to bead locations
xinterp = scatteredInterpolant(topFaceNodes(1,:)',topFaceNodes(2,:)',dispOutput(:,1));
bead_ux = xinterp(bead_x,bead_y);
yinterp = scatteredInterpolant(topFaceNodes(1,:)',topFaceNodes(2,:)',dispOutput(:,2));
bead_uy = yinterp(bead_x,bead_y);

%% Plotting interpolated bead displacments
if ~multiForce
    figure, grid on;
    q2=quiver(bead_xshifted,bead_yshifted,bead_ux,bead_uy);
    xlim([-(meshPtsFwdSol/2) (meshPtsFwdSol/2)]); ylim([-(meshPtsFwdSol/2) (meshPtsFwdSol/2)]);
    xlabel('X'); ylabel('Y'); title('Ground Truth Displacement Field - 1 kPa');
elseif multiForce
    figure, grid on;
    q2=quiver(bead_xshifted,bead_yshifted,bead_ux,bead_uy);
    xlim([-(meshPtsFwdSol/2) (meshPtsFwdSol/2)]); ylim([-(meshPtsFwdSol/2) (meshPtsFwdSol/2)]);
    xlabel('X'); ylabel('Y'); title('Ground Truth Displacement Field - 1 kPa');
end

%% Colormapping
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

%% Create deformed bead image and write to .tif file
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

%% Save results to .mat files for posterity
%Save bead coords and displacements to movieData structure
displField(2).pos = [bead_x, bead_y];
displField(2).vec = [bead_ux, bead_uy];

%Save bead coords and forces to movieData structure
forceField(2).pos = [bead_x, bead_y];
forceField(2).vec = [force_x, force_y];

disp(strcat('Removed ',num2str(nanElements),'NaN elements.'))
save(strcat('outputs_',chr1,'.mat'),'displField','forceField'); %chr1 defined initially by user, matches image file naming