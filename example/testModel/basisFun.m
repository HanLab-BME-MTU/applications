fem.dim = {'v'};
fem.sdim = {'x','y'};

%Load the geometry.
load testGeom.mat; %You get geometry object, 'c'.
fem.geom = c;

%Mesh the domain
fem.mesh = meshinit(fem);

