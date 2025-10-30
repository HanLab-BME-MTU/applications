tic;
disp('Utilizing FEM for forward soluton')

%Hardcoded gel thickness value is temporary
thickness = 15;

vx = linspace(x0(1),x0(2),meshPtsFwdSol);
vy = linspace(y0(1),y0(2),meshPtsFwdSol);
[x_grid,y_grid] = meshgrid(vx,vy);

%Creating PDE Container
femFwdModel = createpde('structural','static-solid');

%Define geometry
gm = multicuboid(x0(2)*2,y0(2)*2,thickness); %multicuboid(x,y,z)
femFwdModel.Geometry = gm;

%Mesh the model
generateMesh(femFwdModel,'Hmax',10);
if length(femFwdModel.Mesh.Nodes(1,:)) < meshPtsFwdSol*h
generateMesh(femFwdModel,'Hmax',8);
end

%Material properties
%FEM will fail when given a Poisson's Ratio of 0.5, therefore we detect
%that case and slightly reduce the ratio to 0.4999
if v < 0.5
structuralProperties(femFwdModel,'YoungsModulus',E,'PoissonsRatio',v);
elseif v == 0.5
structuralProperties(femFwdModel,'YoungsModulus',E,'PoissonsRatio',v-0.001);
end

%Constrain bottom face
structuralBC(femFwdModel,'Face',1,'Constraint','fixed');

%Setting up forces
force_z = zeros(meshPtsFwdSol); %zeros for z forces as we currently are only operating in 2D
%To utilize griddedInterpolant, the input matrices must be NGRID format,
%therefore the dummy variables x_zgrid and y_zgrid are created.
[x_zgrid,y_zgrid] = ndgrid(vx,vy);
forceInterpZ = griddedInterpolant(x_zgrid,y_zgrid,force_z);

%Define variables used to generate the intrpolants that are to be the new
%force_x and force_y variables. This is done because force_x and force_y as
%defined previously do not interact properly with the boundary condition
%function.
allPts   = [0,0;xmin,ymin;xmin,0;0,ymin;0,ymax;xmax,0;xmax,ymax];
dtBaseSup= delaunayTriangulation(allPts);
f_disc_0 = zeros(length(allPts),1);
f_disc_1  = f_disc_0;
f_disc_1(1) = 1;  

%Detect if the fwdSolution function is being run in the X direction or Y
%direction and define the new force_x and force_y variabls accordingly.
if force_x(0,0)==1 && force_y(0,0)==0%this detects X direction
%     f_intp_x= @(location,state) nan2zeroTriScatteredInterp(location.x,location.y,dtBaseSup,f_disc_1,'linear');
%     f_intp_y= @(location,state) nan2zeroTriScatteredInterp(location.x,location.y,dtBaseSup,f_disc_0,'linear'); % only zeros
    
    force_xy = @(location,state)[femFwdInterp(location.x,location.y,allPts,dtBaseSup,f_disc_1,'linear'); ... %force_x is force mesh of x forces
                                 femFwdInterp(location.x,location.y,allPts,dtBaseSup,f_disc_0,'linear'); ... 
                                 forceInterpZ(location.x,location.y);]; %forceInterpz is dummy variable for now

elseif force_x(0,0)==0 && force_y(0,0)==1%this detects Y direction
%     f_intp_x= @(location,state) nan2zeroTriScatteredInterp(location.x,location.y,dtBaseSup,f_disc_0,'linear'); % only zeros
%     f_intp_y= @(location,state) nan2zeroTriScatteredInterp(location.x,location.y,dtBaseSup,f_disc_1,'linear'); 
    
    force_xy = @(location,state)[femFwdInterp(location.x,location.y,allPts,dtBaseSup,f_disc_0,'linear'); ... %only zeros
                                 femFwdInterp(location.x,location.y,allPts,dtBaseSup,f_disc_1,'linear'); ... 
                                 forceInterpZ(location.x,location.y);]; %forceInterpz is dummy variable for now

end
    
%Apply boundary condition to the top face
structuralBoundaryLoad(femFwdModel,'Face',2,'SurfaceTraction',force_xy,'Vectorize','on');

%Solve the structural model
femFwdResults = solve(femFwdModel);

%Interpolate displacements
z_grid = thickness*ones(meshPtsFwdSol);
vz = thickness*ones(1,meshPtsFwdSol);
interpDisp=interpolateDisplacement(femFwdResults,x_grid,y_grid,z_grid);

%Fill in output displacement variables
ux = interpDisp.ux;
ux = reshape(ux,meshPtsFwdSol,[]);
uy = interpDisp.uy;
uy = reshape(uy,meshPtsFwdSol,[]);

%Testing code used for prototyping, unused now
sums = [sum(any(ux~=0)),sum(any(uy~=0)),sum(femFwdResults.Displacement.Magnitude)]
pdeplot3D(femFwdModel,'Deformation',femFwdResults.Displacement,'ColorMap',femFwdResults.VonMisesStress)
                      
toc;    

function vOut=femFwdInterp(x,y,allPts,dtIn,vIn,method)
     F=TriScatteredInterp(allPts,vIn,method);
    vOut = F(x,y);
    checkVec = isnan(vOut);
    vOut(checkVec) = 0;
end