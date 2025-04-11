function [ux,uy,x_grid,y_grid,meshPtsFwdSol]=fwdSolutionFEM(grooveHeight,x0,y0,E,xmin,xmax,ymin,ymax,force_x,force_y,method,opt,meshPtsFwdSol,h,v,refine,useSameSampling)
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

if nargin <14
    v=0.5;
    refine = true;
    useSameSampling = false;
elseif nargin <15
    refine = true;
    useSameSampling = false;
elseif nargin <16
    useSameSampling = false;
end

if strcmpi(method,'FEM')
    % -= GROOVE FEM =-
    tic
    %% Create traction values and substrate bounds
    %set force parameters
    forceOffsetX = 0; forceOffsetY = 0;
    %2 micron = 20 pixels, 0.1 um/pxl
    %set substrate dimensions
    %512 pixel substrate = 51.2 um
    numPix_x = abs(x0(1))+x0(2); numPix_y = abs(y0(1))+y0(2); %#ok<*NASGU>
    %thickness = 800; 
    thickness = x0(2)/2;
    halfSide = x0(2);
    %set groove dimensions
    grooveWidth = 5; %5 groove width = 0.5 um = 500 nm, 50 pxl groove width = 5 um
    %fprintf('Groove width: %1.1f micron \n',grooveWidth/10);
    %grooveHeight = 5;
    
    vx = linspace(x0(1),x0(2),meshPtsFwdSol);
    vy = linspace(y0(1),y0(2),meshPtsFwdSol);
    [x_grid,y_grid] = meshgrid(vx,vy);
    %z_grid = thickness*ones(size(x_grid));
    
    oneORtwo = force_x(0,0);
    if oneORtwo == 1
        %force_xy is defined using function handles for the central basis
        %function
        tractionLoad = @(location,state)[force_x(location.x,location.y); ... %force_x is force mesh of x forces
                                    force_y(location.x,location.y); ... %force_y is force mesh of y forces
                                    force_y(location.x,location.y);]; %forceInterpz is dummy variable for now         
    elseif oneORtwo == 0
        tractionLoad = @(location,state)[force_x(location.x,location.y); ... %force_x is force mesh of x forces
                                    force_y(location.x,location.y); ... %force_y is force mesh of y forces
                                    force_x(location.x,location.y);]; %forceInterpz is dummy variable for now  
    end
    
    [xLoc,yLoc]=meshgrid(linspace(x0(1),x0(2),meshPtsFwdSol+1),linspace(x0(1),x0(2),meshPtsFwdSol+1));
    xLin = xLoc(:); yLin = yLoc(:);

    r = xmax + 1;
    [xCirc,yCirc]=meshgrid(linspace(-r,r,r*4),linspace(-r,r,r*4));
    radiusFuzz = r / 50;
    xc = forceOffsetX; yc = forceOffsetY; n = 100; r = r - radiusFuzz;
    theta = (0:n-1)*(2*pi/n);
    x = xc + r*cos(theta);
    y = yc + r*sin(theta);
    circBoundary = polyshape(x,y);
    xv = circBoundary.Vertices(:,1);
    yv = circBoundary.Vertices(:,2);
    xLinCirc = xCirc(:); yLinCirc = yCirc(:);
    ptsInForceArea = inpolygon(xLinCirc,yLinCirc,xv,yv); 
    coordsInForceArea = [xLinCirc(ptsInForceArea), yLinCirc(ptsInForceArea), thickness*ones(size(xLinCirc(ptsInForceArea)))];
    
    %% Create geometry
    bound = [3; 4; -halfSide; halfSide; halfSide; -halfSide; halfSide; halfSide; -halfSide; -halfSide];
    circ = [1; forceOffsetX; forceOffsetY; r; 0; 0; 0; 0; 0; 0];
    ns = char('bound','circ');
    sf = '(bound+circ)';
    
    %generate groove limits based on geometry size and groove dimensions
    numGrooves = floor(ceil(numPix_x / grooveWidth)/2); %determine number of grooves
    edgePad = mod(numPix_x,grooveWidth)/2; %calculate dead space at edges of substrate where no full groove can fit
    grooveEdges = -halfSide+edgePad-grooveWidth/2:grooveWidth:halfSide-edgePad+grooveWidth/2; %calculate groove edges in x-direction
    if ~mod(numel(grooveEdges), 2) == 0 && grooveWidth ~= 8 %if number of groove edges is odd we have to refine
        grooveEdges = -halfSide+edgePad:grooveWidth:halfSide-edgePad;
        numGrooves = numGrooves - 1;
    end
    if grooveWidth == 20
        grooveEdges = grooveEdges(3:end-2);
        numGrooves = numGrooves - 1;
    elseif grooveWidth == 8
        grooveEdges = grooveEdges(3:end-2);
        numGrooves = numGrooves - 1;
    else
        grooveEdges = grooveEdges(2:end-1); %ensure a groove is placed directly in the center of the substrate
    end
    
    %calculate number of grooves touched by force area
    numPreservedGrooves = 1; %desired number of extra grooves adjacent to outermost force-touched groove
    if grooveWidth / 2 <= xmax || grooveWidth / 2 <= ymax
        overlap = ceil((xmax / grooveWidth)) / 2 - 1;
        numPreservedGrooves = overlap + numPreservedGrooves;
    end
    %numPreservedGrooves must always be 2 or a multiple of 2
    numPreservedGrooves = numPreservedGrooves * 2; 
    if numPreservedGrooves == 1
        numPreservedGrooves = numPreservedGrooves + 1;
    end
    %trim outer grooves beyond force area
    grooveEdges = grooveEdges(length(grooveEdges)/2 - numPreservedGrooves:length(grooveEdges)/2 + 1 + numPreservedGrooves);
    numGrooves = length(grooveEdges) / 2;

    %get point cloud within the boundaries of groove peaks so that the grooves can be extruded in the following steps
    groovePolys{numGrooves} = []; groovePolysUnfuzzed{numGrooves} = []; g = 1;
    edgeFuzz = grooveWidth / 100; %distance to shrink the point cloud to avoid selecting excess faces
    % EDGE FUZZING OF GROOVES MAY CAUSE ISSUES IF THE ADHESION SLIPS JUST BARELY INTO THE NEXT GROOVE AND IS MISSED BY FUZZING AWAY EDGE POINTS
    ptsOnGrooves = false(size(xLin));
    ptsOnGroovesCirc = false(size(xLinCirc));
    for i = 1:2:length(grooveEdges)
        groovePolys{g} = polyshape([grooveEdges(i)+edgeFuzz, grooveEdges(i)+edgeFuzz, grooveEdges(i+1)-edgeFuzz, grooveEdges(i+1)-edgeFuzz],[halfSide-edgeFuzz, -halfSide+edgeFuzz, -halfSide+edgeFuzz, halfSide-edgeFuzz]);
        groovePolysUnfuzzed{g} = polyshape([grooveEdges(i), grooveEdges(i), grooveEdges(i+1), grooveEdges(i+1)],[halfSide, -halfSide, -halfSide, halfSide]);
        xv = groovePolys{g}.Vertices(:,1);
        yv = groovePolys{g}.Vertices(:,2);
        xv2 = groovePolysUnfuzzed{g}.Vertices(:,1);
        yv2 = groovePolysUnfuzzed{g}.Vertices(:,2);
        ptsOnGrooves = ptsOnGrooves | inpolygon(xLin,yLin,xv,yv);
        ptsOnGroovesCirc = ptsOnGroovesCirc | inpolygon(xLinCirc,yLinCirc,xv2,yv2);
        g = g + 1;
    end
    ptsInGrooves = ~ptsOnGrooves; %invert logical array of points at groove peaks to obtain points in groove valleys
    ptsInGroovesCirc = ~ptsOnGroovesCirc;
    
    %obtain coordinates for all points on groove peaks
    coordsOnGrooves = [xLin(ptsOnGrooves), yLin(ptsOnGrooves), thickness*ones(size(yLin(ptsOnGrooves)))-(0.3 * grooveHeight)]; 
    %0.6 is chosen to place the pts ~2/3rds of the groove width above the initial extrusion to ensure inside faces are not accidentally selected
    
    %obtain coordinates for all points in groove valleys and inside the adhesion
    ptsInGroovesAndForceArea = ptsInGroovesCirc & ptsInForceArea;
    coordsInGroovesAndForceArea = [xLinCirc(ptsInGroovesAndForceArea), yLinCirc(ptsInGroovesAndForceArea), thickness*ones(size(yLinCirc(ptsInGroovesAndForceArea)))-(0.5 * grooveHeight)];
    %0.5 is chosen to place the pts exactly halfway between the groove valley floor and groove peak ceiling
    
    %obtain coordinates for all points on groove peaks and inside the adhesion
    ptsOnGroovesAndForceArea = ptsOnGroovesCirc & ptsInForceArea;
    coordsOnGroovesAndForceArea = [xLinCirc(ptsOnGroovesAndForceArea), yLinCirc(ptsOnGroovesAndForceArea), thickness*ones(size(yLinCirc(ptsOnGroovesAndForceArea)))];
    
    grooveEdge(1,:) = grooveEdges(1:2:end); %order groove edges into pairs which define each peak
    grooveEdge(2,:) = grooveEdges(2:2:end);
    grooveEdge = repmat(grooveEdge,2,1); %duplicate both x-coordinates to allow for 4 coordinate pairs once y-coordinates are added
    grooveEdge([3 4],:)=grooveEdge([4 3],:); %swap 3rd and 4th rows to ensure proper ccw point order
    %write groove coordinates into decsg format
    grooveLims = [grooveEdge; repmat([-halfSide; -halfSide; halfSide; halfSide],1,numGrooves)]; %add y-coordinates of groove corner points
    grooves = [repmat(3,1,numGrooves); repmat(4,1,numGrooves); grooveLims]; %add rectangle identifier for decsg
    
    %fill the shape logic and name-space variables with numbered grooves
    if grooveHeight ~= 0
        for i = 1:numGrooves
            sf = [sf,'+groo',num2str(i)]; %#ok<AGROW>
            ns = char(ns,['groo',num2str(i)]);
        end
        gd = [bound,circ,grooves]; %combine geometry descriptions
    else
        gd = [bound,circ];
    end
    ns = ns'; %flip name-space to column orientation
    [dl,~] = decsg(gd,sf,ns); %decompose geometry
    figure,pdegplot(dl,'EdgeLabels','on','FaceLabels','on')
    
    temp = createpde; %create pde to contain geometry
    gtemp = geometryFromEdges(temp,dl); %convert edges into pdetool geometry description
    figure,pdegplot(temp)
    facets = facetAnalyticGeometry(temp,gtemp,0); %grab facets from the geometry
    gm = analyticToDiscrete(facets); %discretize facets to prep for extrusion
    temp.Geometry = gm; %reassociate discretized geometry
    
    pdem = createpde('structural','static-solid'); %create main pde
    g = extrude(gm,thickness-grooveHeight); %initial extrusion of bulk substrate
    pdem.Geometry = g; %reassociate extruded geometry
    figure,pdegplot(pdem,'FaceLabels','on','FaceAlpha',0.5)
    if grooveHeight ~= 0
        if grooveWidth == 50
            grooveFaceIDs = [11,12,14,16];
            g = extrude(g,grooveFaceIDs,grooveHeight); %extrude all grooves to desired height
        else
            grooveFaceIDs = nearestFace(g,coordsOnGrooves); %use centroids to find groove faces IDs
            grooveFaceIDs = unique(grooveFaceIDs);
            g = extrude(g,grooveFaceIDs,grooveHeight); %extrude all grooves to desired height
        end
    end

    pdem.Geometry = g; %reassociate grooved geometry
    figure,pdegplot(pdem,'FaceLabels','on')

    % FLAT SUBSTRATE GEOMETRY
    % bound = [3; 4; -halfSide; halfSide; halfSide; -halfSide; halfSide; halfSide; -halfSide; -halfSide];
    % circ = [1; forceOffsetX; forceOffsetY; d/2; 0; 0; 0; 0; 0; 0];
    % ns = char('bound','circ');
    % sf = '(bound+circ)';
    % 
    % gd = [bound,circ]; %combine geometry descriptions
    % ns = ns'; %flip name-space to column orientation
    % [dl,~] = decsg(gd,sf,ns); %decompose geometry
    % figure,pdegplot(dl,'EdgeLabels','on','FaceLabels','on')
    % 
    % temp = createpde; %create pde to contain geometry
    % gtemp = geometryFromEdges(temp,dl); %convert edges into pdetool geometry description
    % figure,pdegplot(temp)
    % facets = facetAnalyticGeometry(temp,gtemp,0); %grab facets from the geometry
    % gm = analyticToDiscrete(facets); %discretize facets to prep for extrusion
    % temp.Geometry = gm; %reassociate discretized geometry
    % 
    % pdem = createpde('structural','static-solid'); %create main pde
    % g = extrude(gm,thickness); %initial extrusion of bulk substrate
    % pdem.Geometry = g; %reassociate grooved geometry
    % figure, pdegplot(pdem)  

    %% Generate mesh
    try
        generateMesh(pdem,'Hmax',50, 'Hmin',4, 'Hgrad', 2);
    catch
        try
            generateMesh(pdem,'Hmax',35, 'Hmin',3, 'Hgrad', 2);
        catch
            generateMesh(pdem,'Hmax',30, 'Hmin',2, 'Hgrad', 2);
        end
    end
    
    %% Surface traction
    fset2 = nearestFace(g,coordsOnGroovesAndForceArea);
    loadAreaFaceIDs = unique(fset2);
    structuralBoundaryLoad(pdem,'Face',loadAreaFaceIDs,'SurfaceTraction',tractionLoad,'Vectorize','on');
    
    %% Boundary conditions
    coordsOnBase = [xLin, yLin, zeros(size(xLin))];
    constraintFaceIDs = nearestFace(g,coordsOnBase); %use centroids to find groove faces IDs
    constraintFaceIDs = unique(constraintFaceIDs);
    structuralBC(pdem,'Face',constraintFaceIDs,'Constraint','fixed');
    
    %% Substrate material properties
    v = 0.49;
    structuralProperties(pdem,'YoungsModulus',E,'PoissonsRatio',v);
    
    %% Solve model
    pdemResults=solve(pdem);
    disp('Grooved Basis Solution Calculated Successfully!')

    %% Visualize results
    mapPad = 50;
    figure,
    pdeplot3D(pdem,'ColorMapData',pdemResults.Displacement.Magnitude, ...
        'Deformation',pdemResults.Displacement,'DeformationScaleFactor',1);
    view(0,90)
    xlim([128-mapPad 128+mapPad]), ylim([128-mapPad 128+mapPad])
    figure,
    pdeplot3D(pdem,'ColorMapData',pdemResults.Displacement.Magnitude,'Mesh','on')

    %% Interpolate results
    z_mat = thickness*ones(size(x_grid));
    interpdresults = interpolateDisplacement(pdemResults,x_grid,y_grid,z_mat);

    ux = interpdresults.ux;
    ux = reshape(ux,size(x_grid));
    ux(isnan(ux)) = 0;
    figure, imshow(ux,[])

    uy = interpdresults.uy;
    uy = reshape(uy,size(y_grid));
    uy(isnan(uy)) = 0;
    figure, imshow(uy,[])
    
    %% Save fwdSol ux and uy maps
    if grooveHeight ~= 0
        if oneORtwo == 1
            save('FEMOutputsXforce.mat','ux','uy')
        elseif oneORtwo == 0
            save('FEMOutputsYforce.mat','ux','uy')
        end
    else
        if oneORtwo == 1
            save('FEMflatOutputsXforce.mat','ux','uy')
        elseif oneORtwo == 0
            save('FEMflatOutputsYforce.mat','ux','uy')
        end
    end

    toc
end

return;















