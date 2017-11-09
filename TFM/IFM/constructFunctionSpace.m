function  fs = constructFunctionSpace(msh,numSubDoms)
%  fs = constructFunctionSpace(msh) constructs fucntion space from given
%  mesh (msh). 
        fs.fem.dim = {'v'};

        fs.fem.shape = 1;
        fs.fem.bnd.h = 1;
        fs.fem.bnd.r = 0;

        fs.fem.mesh  = msh;
        fs.fem.geom  = fs.fem.mesh;

        fs.fem.equ.ind = ones(1,numSubDoms);
         % I will come back to this point when it's needed SH
%         fs.fem.xmesh = meshextend(fs.fem); % not supported in pde toolbox

        %Find the indices of the DOFs whose shape functions have support disconnected 
        % from the boundary.
%         N = assemble(fs.fem,'out','N');
%         indDomDOF = find(full(sum(N,1))==0); % I will come back to this
%         point when it's needed

        %The Degree of Freedom vector for the basis or shape function in the finite
        % element space.
        %Dimension of the function space where the basis can be nonzero on the 
        % boundary.
%         dimFS  = flngdof(fs.fem);

        %Dimention of the function space where the basis function is kept zero on the
        % boundary.
%         dimBF = length(indDomDOF);

%         fs.dimFS     = dimFS;
%         fs.dimBF     = dimBF;
%         fs.indDomDOF = indDomDOF;
end
