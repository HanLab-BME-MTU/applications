function [ux,uy,x_grid,y_grid,meshPtsFwdSol]=fwdSolution(x0,y0,E,xmin,xmax,ymin,ymax,force_x,force_y,method,opt,meshPtsFwdSol,h,v,refine,useSameSampling)
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
if strcmpi(method,'conv_free')
    tic;
    disp('Calulate the convolution explicitely in free triangulated mesh')
    [nRow,~]=size(x0);

    ux = zeros(nRow,1);
    uy = zeros(nRow,1);
    for i=1:nRow
        integrandx = @(x,y) boussinesqGreens(1,1,x0(i)-x,y0(i)-y,E,v).*force_x(i) + boussinesqGreens(1,2,x0(i)-x,y0(i)-y,E,v).*force_y(i);
        integrandy = @(x,y) boussinesqGreens(2,1,x0(i)-x,y0(i)-y,E,v).*force_x(i) + boussinesqGreens(2,2,x0(i)-x,y0(i)-y,E,v).*force_y(i);

        ux(i) = quad2d(integrandx,xmin,xmax,ymin,ymax,'MaxFunEvals',10^5,'AbsTol',5e-6);% RelTol sucks! 'RelTol',5e-13);
        uy(i) = quad2d(integrandy,xmin,xmax,ymin,ymax,'MaxFunEvals',10^5,'AbsTol',5e-6);% RelTol sucks! 'RelTol',5e-13);
    end
    toc;
    
elseif nargin<10 || strcmpi(method,'conv')
    tic;
    disp('Calulate the convolution explicitely')
    [nRow,nCol]=size(x0);

    for i=1:nRow
        for jj=1:nCol  
            integrandx = @(x,y) boussinesqGreens(1,1,x0(i,jj)-x,y0(i,jj)-y,E,v).*force_x(i,jj) + boussinesqGreens(1,2,x0(i,jj)-x,y0(i,jj)-y,E,v).*force_y(i,jj);
            integrandy = @(x,y) boussinesqGreens(2,1,x0(i,jj)-x,y0(i,jj)-y,E,v).*force_x(i,jj) + boussinesqGreens(2,2,x0(i,jj)-x,y0(i,jj)-y,E,v).*force_y(i,jj);

            ux(i,jj) = quad2d(integrandx,xmin,xmax,ymin,ymax,'MaxFunEvals',10^10,'AbsTol',5e-10);% RelTol sucks! 'RelTol',5e-13);
            uy(i,jj) = quad2d(integrandy,xmin,xmax,ymin,ymax,'MaxFunEvals',10^10,'AbsTol',5e-10);% RelTol sucks! 'RelTol',5e-13);
        end
    end
    toc;
    
elseif strcmpi(method,'FEM')
    nonlin = 1; %nonlinear solver switch, 1 uses nonlin neohookean FEBio, 0 uses linear matlab PDEtool
    % //Work in progress FEM fwd solution method
    if nonlin == 0
    tic;
    disp('Utilizing FEM for forward solution')

    fprintf(1,'Building substrate...')
    vx = linspace(x0(1),x0(2),meshPtsFwdSol);
    vy = linspace(y0(1),y0(2),meshPtsFwdSol);
    [x_grid,y_grid] = meshgrid(vx,vy);

    %Creating PDE Container
    femFwdModel = createpde('structural','static-solid');
    
    %Define geometry
    h = 256; %thickness of substrate
    r = 16; %radius of fine mesh section
    
    %extrude method
    halfSide = x0(2);
    bound = [3; 4; -halfSide; halfSide; halfSide; -halfSide; halfSide; halfSide; -halfSide; -halfSide]; %rectangle of substrate size
    %create circles of decreasing radius in order
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
    g = extrude(gm,h);
    %reassign extruded geometry to model for fwdSolution
    femFwdModel.Geometry = g;
    fprintf(' done.\n') %Done generating substrate
    
    %Mesh the model
    fprintf(1,'Meshing model...')
    generateMesh(femFwdModel,'Hmax',40, 'Hmin',1, 'Hgrad', 1.2);
    fprintf(1,' done.\n') %Done meshing
    
    %Material properties
    fprintf(1,'Assigning BCs and material properties...')
    %FEM will fail when given a Poisson's Ratio of 0.5, therefore we detect
    %that case and slightly reduce the ratio to 0.49
    if v < 0.5
    structuralProperties(femFwdModel,'YoungsModulus',E,'PoissonsRatio',v);
    elseif v == 0.5
    structuralProperties(femFwdModel,'YoungsModulus',E,'PoissonsRatio',v-0.001);
    end
    
    %Constrain bottom face
    structuralBC(femFwdModel,'Face',1:5,'Constraint','fixed');
    
    %Setting up forces
    %force_z = zeros(meshPtsFwdSol); %can be replaced by real data if a normal force is desired
    %To utilize griddedInterpolant, the input matrices must be NGRID format,
    %therefore the dummy variables x_zgrid and y_zgrid are created.
    %[x_zgrid,y_zgrid] = ndgrid(vx,vy);
    %forceInterpZ = griddedInterpolant(x_zgrid,y_zgrid,force_z);
    
    oneORtwo = force_x(0,0);
    if oneORtwo == 1
        %force_xy is defined using function handles for the central basis
        %function
        force_xy = @(location,state)[force_x(location.x,location.y); ... %force_x is force mesh of x forces
                                    force_y(location.x,location.y); ... %force_y is force mesh of y forces
                                    force_y(location.x,location.y);]; %forceInterpz is dummy variable for now         
    elseif oneORtwo == 0
        force_xy = @(location,state)[force_x(location.x,location.y); ... %force_x is force mesh of x forces
                                    force_y(location.x,location.y); ... %force_y is force mesh of y forces
                                    force_x(location.x,location.y);]; %forceInterpz is dummy variable for now  
    end
    
    %Apply loading condition
    structuralBoundaryLoad(femFwdModel,'Face',10,'SurfaceTraction',force_xy,'Vectorize','on');
    fprintf(1,' done.\n') %Done setting BCs and material properties
    
    %Solve model
    fprintf(1,'Solving model...')
    femFwdResults = solve(femFwdModel);
    fprintf(1,' done.\n') %Done solving

    %Interpolate displacements
    fprintf(1,'Interpolating results...')
    z_grid = h*ones(meshPtsFwdSol);
    interpDisp=interpolateDisplacement(femFwdResults,x_grid,y_grid,z_grid);

    %Fill in output displacement variables
    ux = reshape(interpDisp.ux,meshPtsFwdSol,[]);
    uy = reshape(interpDisp.uy,meshPtsFwdSol,[]);
    fprintf(1,' done.\n') %Done interpolating results
    toc;    
    disp('Completed FEM fwd solution.')
    
    save('FEMOutputs.mat','ux','uy')
    
%%%%%NONLINEAR FEM SOLUTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif nonlin == 1
        tic;
        disp('Utilizing Non-linear FEM for forward soluton')

        fprintf(1,'Building substrate...')
        % //Initializing paths and filenames
        %Path names
        folderName = 'FEBio_Data';
        currentFolder = pwd;
        %defaultFolder = fileparts(fileparts(mfilename('fullpath')));
        savePath=fullfile(currentFolder,folderName);

        %Defining file names
        febioFebFileNamePart='nonLinearModel';
        febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
        febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
        febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
        %febioLogFileName_stress=[febioFebFileNamePart,'_stress_out.txt']; %Log file name for exporting stress

        % //Create substrate geometry and tetrahedral mesh
        %Specifying dimensions of substrate
        sampleWidth = x0(2)*2; %Width
        sampleThickness = y0(2)*2; %Thickness
        sampleHeight = x0(2); %Height
        
        %Generating preliminary mesh to use as the basis for the sizing function
        pointSpacing = sampleWidth * 40 / 512; %determine point spacing from experimentally derived constant
        boxDim = [sampleWidth sampleThickness sampleHeight]; %box dimensions
        [F,V,C] = triBox(boxDim,pointSpacing); %preliminary mesh creation
        fprintf(1,' done.\n') %Done generating substrate and surface mesh
        %C = ones(size(F,1),1); %unify region IDs
        
        %Setup tetGen input structure
        V_regions = getInnerPoint(F,V);
        V_holes = [];
        regionTetVolumes = tetVolMeanEst(F,V);
        
        %Generate size field
        %Initialize mask of all nodes on top surface of substrate
        topFaceHeight = max(V(:,3)); %Get height in terms of mesh coords
        edgeSizeField = ones(length(V(:,1)),1); %initialize size field
        topFaceMask = V(:,3) == topFaceHeight; %locate all nodes on top face
        %Assign decreasing edge length to inner radius on top
        distanceFromCenter = sum(bsxfun(@minus,V(:,1:2),[0,0]).^2,2); %calculates the distance from the center to all nodes
        distanceFromCenter = 1.5 * ((distanceFromCenter * ((pointSpacing) / max(distanceFromCenter)))) + 0.5; %scale distance from 2 to max edge size
        edgeSizeField(topFaceMask) = distanceFromCenter(topFaceMask); %assign sizing function to top face
        edgeSizeField(~topFaceMask) = pointSpacing; %set other nodes to original mesh size
        
        %tetGen string modifiers
        stringOpt = '-pq1.2mi';
                    
        %Setup tetGen options
        inputStruct.stringOpt = stringOpt; %tetGen meshing options
        inputStruct.Faces = F; %combined surface faces
        inputStruct.Nodes = V; %combined boundary nodes
        inputStruct.faceBoundaryMarker = C;
        %inputStruct.regionPoints = V_regions; %interior separation points for regions
        inputStruct.holePoints = V_holes; %interior points for hole boundaries
        %inputStruct.regionA = regionTetVolumes; %desired tetrahedral volumes
        inputStruct.sizeData=edgeSizeField; %The size data
        
        %Volume mesh the substrate using tetGen
        meshOutput = runTetGen(inputStruct);
               
        %Write the output data to usable variables for FEBio/GIBBON
        El = meshOutput.elements; %volume mesh elements
        V = meshOutput.nodes; %mesh nodal locations
        CE = meshOutput.elementMaterialID; %region ID for different materials
        Fb = meshOutput.facesBoundary; %boundary faces
        Cb = meshOutput.boundaryMarker; %boundary markers
        
        if isempty(V)
            error('Meshing failed, node list is empty! Try altering tetGen options.')
        end

        %Specifying load type
        %loadType = 'traction';

        %Material Properties
        fprintf('Assigning BCs and material properties...')
        E_youngs1=E; %Material Young's modulus
        if v == 0.5
            nu1=v - 0.001; %Material Poisson's ratio
        elseif v ~= 0.5
            nul=v;
        end

        %FEA control settings
        numTimeSteps=15; %Number of time steps desired
        max_refs=25; %Max reforms
        max_ups=0; %Set to zero to use full-Newton iterations
        opt_iter=6; %Optimum number of iterations
        max_retries=5; %Maximum number of retires
        dtmin=(1/numTimeSteps)/250; %Minimum time step size
        dtmax=1/numTimeSteps; %Maximum time step size

        % //Prepare face nodes for BC application *************************
        %Define supported node sets
        %logicFace=Cb==1; %Logic for current face set
        %Fr=Fb(logicFace,:); %The current face set
        %bcSupportList_X=unique(Fr(:)); %Node set part of selected face

        %logicFace=Cb==3; %Logic for current face set
        %Fr=Fb(logicFace,:); %The current face set
        %bcSupportList_Y=unique(Fr(:)); %Node set part of selected face

        %logicFace=Cb==5; %Logic for current face set
        %Fr=Fb(logicFace,:); %The current face set
        bcSupportList = unique(Fb(Cb==5,:)); %Extract unique nodes from each face in set
        %Since applying BCs only takes the nodal number not the facet of
        %three nodes we select the unique nodes to get a list of support nodes
        
        %bcSupportList = unique(Fb(Cb==1,:));
        %bcSupportList = [bcSupportList; unique(Fb(Cb==2,:))];
        %bcSupportList = [bcSupportList; unique(Fb(Cb==3,:))];
        %bcSupportList = [bcSupportList; unique(Fb(Cb==4,:))];
        %bcSupportList = [bcSupportList; unique(Fb(Cb==5,:))];

        %Prescribed force nodes
        %logicPrescribe=Cb==6; %Logic for current face set
        %Fr=fliplr(Fb(Cb==6,:)); %The current face set
        bcPrescribeList=Fb(Cb==6,:);
        %bcPrescribeList=fliplr(bcPrescribeList);
        bcPrescribeList=unique(bcPrescribeList);
        
        % //Create traction matrix ****************************************
        %Remake topFaceMask
        %topFaceMaskFine = V(:,3) == topFaceHeight; %locate all nodes on top face
        
        %Extract traction value at each surface node
        %topFaceNodeLocations = V((Fb(Cb==6,:)),1:2); %locations of all top nodes
        topFaceNodeLocations = V(bcPrescribeList(:),1:2);
        topFaceNodeLocationX = topFaceNodeLocations(:,1);
        topFaceNodeLocationY = topFaceNodeLocations(:,2);
        C_traction_x = force_x(topFaceNodeLocationX(:),topFaceNodeLocationY(:));
        C_traction_y = force_y(topFaceNodeLocationX(:),topFaceNodeLocationY(:));
        C_traction_z = zeros(length(bcPrescribeList(:)),1);
        
        nodalTractionForce = [C_traction_x(:) C_traction_y(:) C_traction_z(:)];
        
        %**********************
        %bcPrescribeNodalForce=unique(Fb(Cb==6,:));
        %nodal_x = V(bcPrescribeNodalForce,1);
        %nodal_y = V(bcPrescribeNodalForce,2);
        %C_traction_x_nodal = force_x(nodal_x,nodal_y);
        %C_traction_y_nodal = force_y(nodal_x,nodal_y);
        %**********************
        
        debugFig = 1;
        if debugFig == 1
            for i = 1
                %Debugging figures
                %figure,histogram(C_traction_x); ylim([0 5])
                %figure,histogram(C_traction_y);
                
                %Interpolated colormap on top surface
                figure
                x = topFaceNodeLocationX; y = topFaceNodeLocationY; z = C_traction_x;
                dt = delaunayTriangulation(x,y) ;
                tri = dt.ConnectivityList ;
                xi = dt.Points(:,1) ; 
                yi = dt.Points(:,2) ; 
                F = scatteredInterpolant(x,y,z);
                zi = F(xi,yi) ;
                trisurf(tri,xi,yi,zi) 
                hcb=colorbar; colormap(warmcold(250)); caxis([min(C_traction_x) max(C_traction_x)]);
                view(2)
                shading interp
                
                %2D Nodal location colormap on top surface
                figure
                x = topFaceNodeLocationX; y = topFaceNodeLocationY; c = C_traction_x;
                sz = 40;
                scatter(x,y,sz,c,'filled')
                colormap(warmcold(250));

                %Plotting BCs
                fontSize=18;
                %faceAlpha1=0.8;
                markerSize=20;
                markerSize2=15;
                %lineWidth=3;

                hf=cFigure;
                title('Boundary conditions','FontSize',fontSize);
                xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
                hold on;

                gpatch(Fb,V,'kw','k',0.5);

                %hl(1)=plotV(V(bcSupportList,:),'k.','MarkerSize',markerSize);
                
                %if force_x(0,0) == 1
                %    hl(2)=scatterV(V(bcPrescribeList,:),75,C_traction_x,'filled');
                %else
                %    hl(2)=scatterV(V(bcPrescribeList,:),75,C_traction_y,'filled');
                %end

                %legend(hl,{'BC fix support','BC prescribed force'});

                axisGeom(gca,fontSize);
                %hcb=colorbar; colormap(warmcold(250)); caxis([min(C_traction_x) max(C_traction_x)]);
                %yl = ylabel(hcb,'Force Magnitude','FontSize',fontSize);
                %ylp = get(yl,'Position');
                %ylp(1) = 1.4 * ylp(1);
                %set(get(hcb,'ylabel'),'rotation',270,'position',ylp)
                camlight headlight;
                drawnow;
                lol = 1;
            end
        end
        
        % //Define FEBio input structure **********************************
        %Get a template with default settings
        [febio_spec]=febioStructTemplate;

        %febio_spec version
        febio_spec.ATTR.version='3.0';

        %Module section
        febio_spec.Module.ATTR.type='solid';

        %Control section
        febio_spec.Control.analysis='STATIC';
        febio_spec.Control.time_steps=numTimeSteps;
        febio_spec.Control.step_size=1/numTimeSteps;
        febio_spec.Control.solver.max_refs=max_refs;
        febio_spec.Control.solver.max_ups=max_ups;
        febio_spec.Control.time_stepper.dtmin=dtmin;
        febio_spec.Control.time_stepper.dtmax=dtmax;
        febio_spec.Control.time_stepper.max_retries=max_retries;
        febio_spec.Control.time_stepper.opt_iter=opt_iter;

        %Material section
        %Ogden material parameters
        c1 = 1e-3;
        m1 = 8;
        k_factor = 1e2;
        k = c1*k_factor;
        %Assign material parameters
        materialName1='Material1';
        febio_spec.Material.material{1}.ATTR.name=materialName1;
        
        materialType = 0;
        switch materialType
            case 0
                %Neo-Hookean Material
                febio_spec.Material.material{1}.ATTR.type='neo-Hookean';
                febio_spec.Material.material{1}.ATTR.id=1;
                febio_spec.Material.material{1}.E=E_youngs1;
                febio_spec.Material.material{1}.v=nu1;
            case 1
                %Ogden Material
                febio_spec.Material.material{1}.ATTR.type='Ogden';
                febio_spec.Material.material{1}.ATTR.id=1;
                febio_spec.Material.material{1}.c1=c1;
                febio_spec.Material.material{1}.m1=m1;
                febio_spec.Material.material{1}.c2=c1;
                febio_spec.Material.material{1}.m2=-m1;
                febio_spec.Material.material{1}.k=k;
        end

        % Mesh section
        % -> Nodes
        febio_spec.Mesh.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
        febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
        febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates

        % -> Elements
        partName1='Part1';
        febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
        febio_spec.Mesh.Elements{1}.ATTR.type='tet4'; %Element type
        febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(El,1))'; %Element id's
        febio_spec.Mesh.Elements{1}.elem.VAL=El; %The element matrix

        % -> NodeSets
        %nodeSetName2='bcSupportList_X';
        %nodeSetName3='bcSupportList_Y';
        nodeSetName1='bcSupportList';
        %nodeSetName4='bcPrescribeList';
        
        febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
        febio_spec.Mesh.NodeSet{1}.node.ATTR.id=bcSupportList(:);

        %febio_spec.Mesh.NodeSet{2}.ATTR.name=nodeSetName2;
        %febio_spec.Mesh.NodeSet{2}.node.ATTR.id=bcSupportList_X(:);

        %febio_spec.Mesh.NodeSet{3}.ATTR.name=nodeSetName3;
        %febio_spec.Mesh.NodeSet{3}.node.ATTR.id=bcSupportList_Y(:);
        
        %febio_spec.Mesh.NodeSet{2}.ATTR.name=nodeSetName2;
        %febio_spec.Mesh.NodeSet{2}.node.ATTR.id=bcPrescribeList;
        
        %**********************
        nodeSetName2='bcPrescribeList';
        febio_spec.Mesh.NodeSet{2}.ATTR.name=nodeSetName2;
        febio_spec.Mesh.NodeSet{2}.node.ATTR.id=bcPrescribeList(:);
        %**********************

        % -> Surfaces
        %surfaceName1='LoadedSurface';
        %febio_spec.Mesh.Surface{1}.ATTR.name=surfaceName1;
        %febio_spec.Mesh.Surface{1}.tri3.ATTR.id=(1:1:size(bcPrescribeList,1))';
        %There is some debate regarding using Fr over bcPrescribe here
        %febio_spec.Mesh.Surface{1}.tri3.VAL=bcPrescribeList;%Fr;

        %MeshDomains section
        febio_spec.MeshDomains.SolidDomain.ATTR.name=partName1;
        febio_spec.MeshDomains.SolidDomain.ATTR.mat=materialName1;

        %Boundary condition section
        % -> Fix boundary conditions
        febio_spec.Boundary.bc{1}.ATTR.type='fix';
        febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
        febio_spec.Boundary.bc{1}.dofs='x,y,z';
        
        %febio_spec.Boundary.bc{2}.ATTR.type='fix';
        %febio_spec.Boundary.bc{2}.ATTR.node_set=nodeSetName2;
        %febio_spec.Boundary.bc{2}.dofs='x';

        %febio_spec.Boundary.bc{3}.ATTR.type='fix';
        %febio_spec.Boundary.bc{3}.ATTR.node_set=nodeSetName3;
        %febio_spec.Boundary.bc{3}.dofs='y';

        %MeshData section
        %************************
        % -> Node data
        loadDataName1='traction_force';
        febio_spec.MeshData.NodeData{1}.ATTR.name=loadDataName1;
        febio_spec.MeshData.NodeData{1}.ATTR.node_set=nodeSetName2;
        febio_spec.MeshData.NodeData{1}.ATTR.datatype='vec3';
        febio_spec.MeshData.NodeData{1}.node.ATTR.lid=(1:1:numel(bcPrescribeList))';
        febio_spec.MeshData.NodeData{1}.node.VAL=nodalTractionForce;
        
        % -> Prescribed nodal forces
        febio_spec.Loads.nodal_load{1}.ATTR.name='topface_nodal_force';
        febio_spec.Loads.nodal_load{1}.ATTR.type='nodal_force';
        febio_spec.Loads.nodal_load{1}.ATTR.node_set=nodeSetName2;
        febio_spec.Loads.nodal_load{1}.value.ATTR.lc=1;
        febio_spec.Loads.nodal_load{1}.value.ATTR.type='map';
        febio_spec.Loads.nodal_load{1}.value.VAL=loadDataName1;
        %************************
        
        % -> Surface data
        %loadDataName1='LoadData1';
        %febio_spec.MeshData.SurfaceData{1}.ATTR.name=loadDataName1;
        %febio_spec.MeshData.SurfaceData{1}.ATTR.surface=surfaceName1;
        %febio_spec.MeshData.SurfaceData{1}.ATTR.datatype='vec3';
        %febio_spec.MeshData.SurfaceData{1}.face.ATTR.lid=(1:1:numel(C_traction_x))';
        %febio_spec.MeshData.SurfaceData{1}.face.VAL=[C_traction_x(:) C_traction_y(:) zeros(size(C_traction_x(:)))];

        %Loads section
        % -> Surface load
        %febio_spec.Loads.surface_load{1}.ATTR.type='traction';
        %febio_spec.Loads.surface_load{1}.ATTR.surface=surfaceName1;
        %febio_spec.Loads.surface_load{1}.traction.ATTR.lc=1;
        %febio_spec.Loads.surface_load{1}.traction.ATTR.type='map';
        %febio_spec.Loads.surface_load{1}.traction.VAL=loadDataName1;
        
        %LoadData section
        % -> load_controller
        febio_spec.LoadData.load_controller{1}.ATTR.id=1;
        febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
        febio_spec.LoadData.load_controller{1}.interpolate='LINEAR';
        %Conduct a run with all points set to 1
        febio_spec.LoadData.load_controller{1}.points.point.VAL=[0 0; 1 1];

        %Output section
        % -> log file
        febio_spec.Output.logfile.ATTR.file=febioLogFileName;
        febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
        febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
        febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
        febio_spec.Output.logfile.node_data{1}.VAL=1:size(V,1);

        %febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress;
        %febio_spec.Output.logfile.element_data{1}.ATTR.data='s1';
        %febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';
        %febio_spec.Output.logfile.element_data{1}.VAL=1:size(V,1);

        % //FEBio running *************************************************
        %Run Settings
        febioAnalysis.run_filename=febioFebFileName; %The input file name
        febioAnalysis.run_logname=febioLogFileName; %The name for the log file
        febioAnalysis.disp_on=1; %Display information on the command window
        febioAnalysis.runMode='internal';%'internal';preferred
        febioAnalysis.maxLogCheckTime=120; %Max log file checking time
        optionStruct.arrayParseMethod = 1; %1,2, or 3 (1 preferred)
        
        %Check and create output folders
        if ~exist(folderName,'dir')
            mkdir(folderName)
        end
        fprintf(1,' done.\n') %Done defining BCs and material conditions

        %Export FEBio Structure
        [domNode] = febioStruct2xml(febio_spec,febioFebFileName,optionStruct); %Exporting to file
        %system(['gedit ',febioFebFileName,' &']);

        %Run Model
        skipFEBio = 0;
        if skipFEBio ~= 1
        [runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!
        end

        % //Import Displacement from Logfile ******************************
        if runFlag == 1
        dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp),1,1);

        N_disp_mat=dataStruct.data; % nNodes-by-(u,v,w)-by-nTimeSteps
        timeVec=dataStruct.time; %Time
        
        % //Plotting results **********************************************
        drawFEBioFig = 0;
        if drawFEBioFig == 1
            for i = 1
        %Create deformed coordinate set
        V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);
        DN_magnitude=sqrt(sum(N_disp_mat(:,:,end).^2,2)); %Current displacement magnitude
        % Create basic view and store graphics handle to initiate animation
        hf=cFigure; %Open figure
        gtitle([febioFebFileNamePart,': Press play to animate']);
        title('Displacement magnitude [mm]','Interpreter','Latex')
        hp=gpatch(Fb,V_DEF(:,:,end),DN_magnitude,'k',1,2); %Add graphics object to animate
        hp.Marker='.';
        hp.MarkerSize=markerSize2;
        hp.FaceColor='interp';
        gpatch(Fb,V,0.5*ones(1,3),'none',0.25); %A static graphics object

        axisGeom(gca,fontSize);
        colormap(cMap); colorbar;
        caxis([0 max(DN_magnitude)]); caxis manual;
        axis(axisLim(V_DEF)); %Set axis limits statically
        view(140,30);
        camlight headlight;

        % Set up animation features
        animStruct.Time=timeVec; %The time vector
        for qt=1:1:size(N_disp_mat,3) %Loop over time increments
            DN_magnitude=sqrt(sum(N_disp_mat(:,:,qt).^2,2)); %Current displacement magnitude

            %Set entries in animation structure
            animStruct.Handles{qt}=[hp hp]; %Handles of objects to animate
            animStruct.Props{qt}={'Vertices','CData'}; %Properties of objects to animate
            animStruct.Set{qt}={V_DEF(:,:,qt),DN_magnitude}; %Property values for to set in order to animate
        end
        anim8(hf,animStruct); %Initiate animation feature
        drawnow;
            end
        end

        % //Format FEBio output into fwdSolution compatible form **********
        %Component displacement at final time step
        dispMat = N_disp_mat(:,:,end);
        node_ux = dispMat(:,1);
        node_uy = dispMat(:,2);

        %Extract top face nodal displacements
        topNode_ux = node_ux(bcPrescribeList(:));
        topNode_uy = node_uy(bcPrescribeList(:));
        %topNode_ux = node_ux(V(:,3) == topFaceHeight);
        %topNode_uy = node_uy(V(:,3) == topFaceHeight);
        
        %Define interpolation locations
        vx = linspace(x0(1),x0(2),meshPtsFwdSol);
        vy = linspace(y0(1),y0(2),meshPtsFwdSol);
        [x_grid,y_grid] = meshgrid(vx,vy);

        reformatOut = 0;
        %Reformatting output matrices to ensure correct orientation
        if reformatOut == 1
            %Initialize loop variables
            nodeDispMap_ux = zeros(130);
            nodeDispMap_uy = zeros(130);
            outMatWidth = numElementsWidth + 1;
            startID = 1;
            endID = outMatWidth;
            idx = 1;
            %Loop through every column of the output matrices to place the
            %displacements in the correct order
            for k = 1:outMatWidth
                %Get current column of UX and UY displacement
                %Rotate 180 to orient Y coordinates in proper order
                curXcol = rot90(topNode_ux(startID:endID),2);
                curYcol = rot90(topNode_uy(startID:endID),2);
                %Place current columns into output matrices 1 column at a time
                %from left to right
                nodeDispMap_ux(1:outMatWidth,idx) = curXcol;
                nodeDispMap_uy(1:outMatWidth,idx) = curYcol;
                %Increment counters
                startID = startID + outMatWidth;
                endID = endID + outMatWidth;
                idx = idx + 1;
            end
        end
        
        %x = topFaceNodeLocationX ; y = topFaceNodeLocationY ; z = topNode_ux ;
        %dt = delaunayTriangulation(x,y) ;
        %tri = dt.ConnectivityList ;
        %xi = dt.Points(:,1) ; 
        %yi = dt.Points(:,2) ; 
        %F = scatteredInterpolant(x,y,z);
        %zi = F(xi,yi) ;
        %trisurf(tri,xi,yi,zi) 
        %view(2)
        %shading interp

        %Interpolate disp map at TFM nodal locations into output variables
        [~,~,ux] = griddata(topFaceNodeLocationX,topFaceNodeLocationY,topNode_ux,x_grid,y_grid);
        [~,~,uy] = griddata(topFaceNodeLocationX,topFaceNodeLocationY,topNode_uy,x_grid,y_grid);
        end
        toc;    
        disp('Completed Non-linear FEM fwd solution.')
        switch materialType
            case 0
                save('FEMNonLinNeoHookeanOutputs.mat','ux','uy')
            case 1
                save('FEMNonLinOgdenOutputs.mat','ux','uy')
        end     
    end

elseif strcmpi(method,'fft')
    %display('Use fast convolution')    

    %***************************************************************
    % Here starts the calculation using the fast fourier transform *
    %***************************************************************
    
    % Number of points to calculate the force field and the Greensfunction    % Since below Nx_G and Ny_G are odd, the Greensfunctions will be
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
    if (nargin < 12 || isempty(meshPtsFwdSol)) && ~useSameSampling
        disp('Use meshPtsFwdSol=2^10. This value should be given with the function call!!!');
        meshPtsFwdSol=2^10;
    end
    
    if useSameSampling
        Nx_F=size(x0,2); % 2^10 is the densest sampling possible.
        Ny_F=size(y0,1);
    else
        Nx_F=meshPtsFwdSol; % 2^10 is the densest sampling possible.
        Ny_F=Nx_F;
    end
    
    % To account for dx*dy in the convolution integral one has to finally
    % rescale the result by the following scaling factor:
    xRange=(max(max(x0))-min(min(x0)));
    yRange=(max(max(y0))-min(min(y0)));
    scalingFactor=(xRange*yRange)/(Nx_F*Ny_F);
%     scalingFactor=1;
    
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
    [xgrid_F,ygrid_F]=meshgrid(xvec_F,yvec_F);
    
    % create a mesh centered at zero with Nx_G*Ny_G meshpoints, where the
    % Greensfunctions are calculated.
    xvec_G=linspace(-xRange,xRange,Nx_G);
    yvec_G=linspace(-yRange,yRange,Ny_G);
    
%     [xgrid_G,ygrid_G]=meshgrid(yvec_G,xvec_G);
    [xgrid_G,ygrid_G]=meshgrid(xvec_G,yvec_G);
      
    %calculate the force values at the grid_F positions:
    if useSameSampling
        discrete_Force_x_unPadded=force_x; %this has only to be calculated over the support xmin,xmax,ymin,ymax rest is zero
        discrete_Force_y_unPadded=force_y;
    else
    % Make force_x and force_y a function handle if it is a matrix
        if ismatrix(force_x) && ismatrix(force_y) && ~isa(force_x,'TriScatteredInterp') && ~isa(force_x,'function_handle')
        %     [xmat,ymat]=meshgrid(xmin:xmax,ymin:ymax);
        %     xvec=xmat(:);
        %     yvec=ymat(:);
            xvec=x0(:);
            yvec=y0(:);
            force_x_vec=force_x(:);
            force_y_vec=force_y(:);
            force_x = scatteredInterpolant(xvec,yvec,force_x_vec);
            force_y = scatteredInterpolant(xvec,yvec,force_y_vec);
        end
        discrete_Force_x_unPadded=force_x(xgrid_F,ygrid_F); %this has only to be calculated over the support xmin,xmax,ymin,ymax rest is zero
        discrete_Force_y_unPadded=force_y(xgrid_F,ygrid_F);
        % taking care of nans
        checkVec=isnan(discrete_Force_x_unPadded);
        discrete_Force_x_unPadded(checkVec)=0;
        checkVec=isnan(discrete_Force_y_unPadded);
        discrete_Force_y_unPadded(checkVec)=0;
    end
    % Calculate the Greens-function values at the grid_G positions. This can
    % be improved since the Greensfunction never change for a given grid
    % size. When the Basis functions are calculated this has to be done
    % only once (as well as the FFT for these fields!!!):
    discrete_boussinesqGreens11=boussinesqGreens(1,1,xgrid_G,ygrid_G,E,v);
    discrete_boussinesqGreens12=boussinesqGreens(1,2,xgrid_G,ygrid_G,E,v);
   %discrete_boussinesqGreens21=discrete_boussinesqGreens12;
    discrete_boussinesqGreens22=boussinesqGreens(2,2,xgrid_G,ygrid_G,E,v);
    
    % Pad the calculated fields with zero to the next power larger than 
    % (2*N-1), see above. For this setup, the FFT is fastest.
%     discrete_Force_x=padarray(discrete_Force_x_unPadded,[Nx_pad-Nx_F Ny_pad-Ny_F],0,'post');%'symmetric','post');
%     discrete_Force_y=padarray(discrete_Force_y_unPadded,[Nx_pad-Nx_F Ny_pad-Ny_F],0,'post');
    discrete_Force_x=padarray(discrete_Force_x_unPadded,[Ny_pad-Ny_F Nx_pad-Nx_F],0,'post');%'symmetric','post');
    discrete_Force_y=padarray(discrete_Force_y_unPadded,[Ny_pad-Ny_F Nx_pad-Nx_F],0,'post');
    
%     discrete_boussinesqGreens11=padarray(discrete_boussinesqGreens11,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
%     discrete_boussinesqGreens12=padarray(discrete_boussinesqGreens12,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
%    %discrete_boussinesqGreens22=padarray(discrete_boussinesqGreens22,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post'); 
%     discrete_boussinesqGreens22=padarray(discrete_boussinesqGreens22,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
    
    discrete_boussinesqGreens11=padarray(discrete_boussinesqGreens11,[Ny_pad-Ny_G Nx_pad-Nx_G],0,'post');
    discrete_boussinesqGreens12=padarray(discrete_boussinesqGreens12,[Ny_pad-Ny_G Nx_pad-Nx_G],0,'post');
    discrete_boussinesqGreens22=padarray(discrete_boussinesqGreens22,[Ny_pad-Ny_G Nx_pad-Nx_G],0,'post');
    
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
    
    % Now calculate the solution:                
    ux_grid=ifft2(dFT_boussinesqGreens11.*dFT_Force_x+dFT_boussinesqGreens12.*dFT_Force_y);
    clear dFT_boussinesqGreens11 dFT_boussinesqGreens12;
    uy_grid=ifft2(dFT_boussinesqGreens21.*dFT_Force_x+dFT_boussinesqGreens22.*dFT_Force_y);
    clear dFT_boussinesqGreens21 dFT_Force_x dFT_boussinesqGreens22 dFT_Force_y;
    
    % Plot the solution:
%     figure(10)
%     imshow(ux_grid,[])
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
%     ux_grid=real(ux_grid(startIndex_x:endIndex_x,startIndex_y:endIndex_y));
%     uy_grid=real(uy_grid(startIndex_x:endIndex_x,startIndex_y:endIndex_y));
    ux_grid=real(ux_grid(startIndex_y:endIndex_y,startIndex_x:endIndex_x));
    uy_grid=real(uy_grid(startIndex_y:endIndex_y,startIndex_x:endIndex_x));

%!!! This could be improved by using the analytical solution for the Fourie
%!!! Transform of the Greensfunction!
    % Add the solution for G(0,0). This is a correction term which becomes
    % irrelevant for very dense sampling. But for small Nx_F it is REALLY
    % essential!
    % Set the Poisson's ratio to 0.5:
    v=0.5;
    dx=abs(xvec_G(2)-xvec_G(1));
    dy=abs(yvec_G(2)-yvec_G(1));
    
    int_x2_over_r3=2*dy*log((dy^2+2*dx*(dx+sqrt(dx^2+dy^2)))/(dy^2));    
    int_y2_over_r3=2*dx*log((dx^2+2*dy*(dy+sqrt(dx^2+dy^2)))/(dx^2));    
    int_1_over_r  =int_x2_over_r3 + int_y2_over_r3;
        
    corrTerm_11=(1+v)/(pi*E)*((1-v)*int_1_over_r+v*int_x2_over_r3);
    corrTerm_22=(1+v)/(pi*E)*((1-v)*int_1_over_r+v*int_y2_over_r3);
    
    ux_grid=ux_grid+discrete_Force_x_unPadded*corrTerm_11;
    clear discrete_Force_x_unPadded;
    uy_grid=uy_grid+discrete_Force_y_unPadded*corrTerm_22;
    clear discrete_Force_y_unPadded;
    
    % scale the solution appropriately!
    ux_grid=scalingFactor*ux_grid;
    uy_grid=scalingFactor*uy_grid;
    
    
    
    % Recursive call to fwdSolution:   
%     refine =true;
    if ~isempty(xmin) && ~isempty(xmax) && ~isempty(ymin) && ~isempty(ymax) && refine
        % and the support of the force is much small than ROI for the
        % displacement, this check should be included otherwise one calculate
        % the same thing twice
        % extend the size of the region a little bit, here by a factor of 1.2,
        % but this is arbitrary.
        % The range of the support is:
        xsupp=xmax-xmin;
        ysupp=ymax-ymin;
        
        % Now expand the x- and y-range:
        xminExp=xmin-xsupp/10;
        xmaxExp=xmax+xsupp/10;
        yminExp=ymin-ysupp/10;
        ymaxExp=ymax+ysupp/10;
        
        [ux_fine uy_fine x_grid_fine y_grid_fine]=fwdSolution([xminExp xmaxExp],[yminExp ymaxExp],E,[],[],[],[],force_x,force_y,method,'noIntp',meshPtsFwdSol);
        
        % Later on we want to have ux and uy defined on a regular
        % grid. for this reason we now interpolate the fine solution onto
        % the regular xgrid_F and ygrid_F. This grid is usually so fine
        % that it already corresponds to subpixel sampling: (e.g. 2^10=1024
        % positions along each image dimension)
        
        % find the positions that are within
        
        iux_fine = griddata(x_grid_fine,y_grid_fine,ux_fine,xgrid_F,ygrid_F,'linear');%'*cubic'
        iuy_fine = griddata(x_grid_fine,y_grid_fine,uy_fine,xgrid_F,ygrid_F,'linear');%'*linear'
        
        % those contain a lot of NaNs since most of the points are outside
        % of [xminExp xmaxExp] and [yminExp ymaxExp]. This will give us two
        % masks that we can now use to update ux_grid and uy_grid:        
        valMat=~isnan(iux_fine);
        
        diff=2*abs(uy_grid(valMat)-iuy_fine(valMat))./abs(uy_grid(valMat)+iuy_fine(valMat));
        display(['max. correction: ',num2str(max(diff(:)))]);
        
        % Update the values in the coarse solution:
        ux_grid(valMat)=iux_fine(valMat);
        uy_grid(valMat)=iuy_fine(valMat);
    end
    
        
    % Now interpolate the displacement field from the regular grid to the irregular
    % measured grid. Since we have the force field defined on a regular grid
    % we can use the fast *option. 'linear' is about a factor of two faster than
    % 'cubic'. Hard to tell if cubic performs better than linear.    
    if nargin>10 && strcmp(opt,'noIntp')
        ux = ux_grid;%'*cubic'
        uy = uy_grid;%'*linear'
        x_grid=xgrid_F;
        y_grid=ygrid_F;
    else
        ux = interp2(xgrid_F,ygrid_F,ux_grid,x0,y0,'*cubic');%'*cubic'
        uy = interp2(xgrid_F,ygrid_F,uy_grid,x0,y0,'*cubic');%'*linear'
        x_grid=x0;
        y_grid=y0;
    end
    
    save('fastBEMoutputs.mat','ux','uy')
    
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
elseif strcmpi(method,'fft_finite')
    %display('Use fast convolution')    

    % This determines the sampling of the force field:
    if nargin < 12 || isempty(meshPtsFwdSol)
        display('Use meshPtsFwdSol=2^10. This value should be given with the function call!!!');
        meshPtsFwdSol=2^10;
    end
        
    Nx_F=meshPtsFwdSol; % 2^10 is the densest sampling possible.
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
    [xgrid_F,ygrid_F]=meshgrid(xvec_F,yvec_F);
    
    % create a mesh centered at zero with Nx_G*Ny_G meshpoints, where the
    % Greensfunctions are calculated.
    xvec_G=linspace(-xRange,xRange,Nx_G);
    yvec_G=linspace(-yRange,yRange,Ny_G);
    
    [xgrid_G,ygrid_G]=meshgrid(xvec_G,yvec_G);
      
    %calculate the force values at the grid_F positions:
    discrete_Force_x_unPadded=force_x(xgrid_F,ygrid_F); %this has only to be calculated over the support xmin,xmax,ymin,ymax rest is zero
    discrete_Force_y_unPadded=force_y(xgrid_F,ygrid_F);

    % This part if what's different from 'fft' method that uses
    % boussinesque greens function
    discrete_boussinesqGreens11=finiteThicknessGreens(1,1,xgrid_G,ygrid_G,E,h);
    discrete_boussinesqGreens12=finiteThicknessGreens(1,2,xgrid_G,ygrid_G,E,h);
   %discrete_boussinesqGreens21=discrete_boussinesqGreens12;
    discrete_boussinesqGreens22=finiteThicknessGreens(2,2,xgrid_G,ygrid_G,E,h);
    
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
    
    % Now calculate the solution:                
    ux_grid=ifft2(dFT_boussinesqGreens11.*dFT_Force_x+dFT_boussinesqGreens12.*dFT_Force_y);
    clear dFT_boussinesqGreens11 dFT_boussinesqGreens12;
    uy_grid=ifft2(dFT_boussinesqGreens21.*dFT_Force_x+dFT_boussinesqGreens22.*dFT_Force_y);
    clear dFT_boussinesqGreens21 dFT_Force_x dFT_boussinesqGreens22 dFT_Force_y;
    
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

%     figure
%     imshow(ux_grid,[])
    % scale the solution appropriately!
    ux_grid=scalingFactor*ux_grid;
    uy_grid=scalingFactor*uy_grid;
    
    % Now interpolate the displacement field from the regular grid to the irregular
    % measured grid. Since we have the force field defined on a regular grid
    % we can use the fast *option. 'linear' is about a factor of two faster than
    % 'cubic'. Hard to tell if cubic performs better than linear.    
    if nargin>10 && strcmp(opt,'noIntp')
        ux = ux_grid;%'*cubic'
        uy = uy_grid;%'*linear'
        x_grid=xgrid_F;
        y_grid=ygrid_F;
    else
        ux = interp2(xgrid_F,ygrid_F,ux_grid,x0,y0,'*cubic');%'*cubic'
        uy = interp2(xgrid_F,ygrid_F,uy_grid,x0,y0,'*cubic');%'*linear'
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

% strange!!!:
%x0_vec=linspace(-5,20,30);
%y0_vec=linspace(1,6,15);

x0_vec=linspace(-50,200,51);
y0_vec=linspace(0,600,61);

[x0 y0]=meshgrid(x0_vec,y0_vec);

E=10000;
meshPtsFwdSol=2^10;

xmin=min(x0_vec);
xmax=max(x0_vec);
ymin=min(y0_vec);
ymax=max(y0_vec);

xmin=0;
xmax=2;
ymin=1;
ymax=3;

force_x=@(x,y) assumedForce(1,x,y);
force_y=@(x,y) assumedForce(2,x,y);
force_x =@(x,y)  assumedForceAniso2D(1,x,y,70,28,150,620,400/72,500/72,'groupForce');
force_y =@(x,y)  assumedForceAniso2D(2,x,y,70,28,150,620,400/72,500/72,'groupForce');


[ux_conv uy_conv]=fwdSolution(x0,y0,E,xmin,xmax,ymin,ymax,force_x,force_y,'conv');

figure(1)
quiver(x0,y0,ux_conv,uy_conv);


[ux_fft uy_fft]=fwdSolution(x0,y0,E,xmin,xmax,ymin,ymax,force_x,force_y,'fft',[],meshPtsFwdSol);

% figure(11)
% quiver(x0,y0,ux_fft,uy_fft,'r');

%compare the two results:
scalePlot=0.3*sqrt(max(max(ux_conv.^2+uy_conv.^2)));
figure(40)
quiver(x0,y0,ux_fft/scalePlot,uy_fft/scalePlot,0,'r');
hold on
quiver(x0,y0,ux_conv/scalePlot,uy_conv/scalePlot,0,'g');
hold off

%figure(5)
%quiver(x0,y0,assumedForce(1,x0,y0),assumedForce(2,x0,y0),0,'g');

corr_x=2*abs(ux_fft-ux_conv)./abs(ux_fft+ux_conv);
corr_y=2*abs(uy_fft-uy_conv)./abs(uy_fft+uy_conv);

figure(3)
surf(x0,y0,corr_x)

figure(4)
surf(x0,y0,corr_y)

display(['mean rel. deviation in %: ',num2str(mean(corr_y(:)))]);
display([' max rel. deviation in %: ',num2str(max(corr_y(:)))]);
%uncorr_x=2*abs(ux_fft-ux_conv)./abs(ux_fft+ux_conv)
%uncorr_y=2*abs(uy_fft-uy_conv)./abs(uy_fft+uy_conv)

% uncorr_y-corr_y

display('This should be 0')



% in case of:
% x0_vec=linspace(-10,10,25);
% y0_vec=linspace(-2,2,15);
% Important note: for a sampling of 2^8 in the fast solution,
% mean(corr_y(:)) decreases for increasing numerical precision in the
% direct integration of the fwd solution (up to 'AbsTol' 10^-(8)). Thus, it seems 
% like as if 2^8 sampling in the Fourier solution is almost as precise as what
% can be achieved by direct numerical intergeation (up to 'AbsTol' 10^-(10))... which is amazing!!!
% Fourier sampling: 2^(8)
% AbsTol:       mean(corr_y(:))
% 10^(-10)      0.0067
% 10^(- 9)      0.0070
% 10^(- 8)      0.0074  (here precision: fast sol ~= direct conv)
% 10^(- 7)      0.0210
% 10^(- 6)      0.0442
% 10^(- 5)      0.1079
% 10^(- 4)      0.1179

% Fourier sampling: 2^(10)
% AbsTol:       mean(corr_y(:))
% 10^(-10)      0.0011
% 10^(- 9)      0.0013  (here precision: fast sol ~= direct conv)
% 10^(- 8)      0.0023
% 10^(- 7)      0.0162
% 10^(- 6)      0.0391
% 10^(- 5)      0.1024
% 10^(- 4)      0.1123

% Note also that the best precision is where the displacement is predicted
% to be highest, that is at the force center! The mean rel. difference
% between sampling 2^(8) and 2^(10) is 0.0056 which roughly indicates the
% accuracy of the 2^(8) sampling, which is inline with the results above.

















