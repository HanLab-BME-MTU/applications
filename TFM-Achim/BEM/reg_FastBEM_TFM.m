function [pos_f force forceMesh M pos_u u sol_coef]=reg_FastBEM_TFM(grid_mat, displField, frame, yModu_Pa, pRatio, regParam, meshPtsFwdSol)
if nargin < 7
    meshPtsFwdSol=[];
end

if isempty(pRatio) || pRatio~=0.5
    display('The poisson ratio has been set to 0.5, the script has not been generalized yet!');
end

if isempty(grid_mat)
    % If no mesh is specified for the forces, we create a hexagonal mesh
    % that will be centered in the field of view. Here, only the first
    % frame will be used:
    [xvec,yvec]=createHexGridFromDisplField(displField,1);

    % plot the grid points:
    % figure(1)
    % plot(xvec,yvec,'o');
else
    xvec=grid_mat(:,:,1);
    yvec=grid_mat(:,:,2);
    
    xvec=xvec(:);
    yvec=yvec(:);
end

display('1.) Creating mesh & basis [~30sec]:...');
tic;
doPlot=0;
forceMesh=createMeshAndBasisFastBEM(xvec,yvec,0,[],doPlot);
toc;
display('Done: mesh & basis!');

[fx fy x_out y_out M pos_u u sol_coef] = BEM_force_reconstruction(displField(frame).pos(:,1),displField(frame).pos(:,2),displField(frame).vec(:,1),displField(frame).vec(:,2),forceMesh,yModu_Pa,regParam,[],[],'fast',meshPtsFwdSol);
% The units of fx and fy are the same as the input E, that is ususally Pa!

pos_f  =horzcat(x_out,y_out);
force=horzcat(   fx,   fy);


figure(100)
quiver(x_out,y_out,fx,fy,'b')
hold on
quiver(displField(1).pos(:,1),displField(1).pos(:,2),displField(1).vec(:,1),displField(1).vec(:,2),'r')
hold off