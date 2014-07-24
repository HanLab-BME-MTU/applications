function [M]=calcFwdMapFastBEM3D(x_vec_u, y_vec_u, z, forceMesh, E, v)

numBasis=length(forceMesh.basis);

% get the boundary for the force mesh:
xmin_b=forceMesh.basis(1).node(1);
xmax_b=forceMesh.basis(1).node(1);
ymin_b=forceMesh.basis(1).node(2);
ymax_b=forceMesh.basis(1).node(2);
for j=1:numBasis
  if forceMesh.basis(j).node(1)<xmin_b
    xmin_b=forceMesh.basis(j).node(1);
  end
  if forceMesh.basis(j).node(1)>xmax_b
    xmax_b=forceMesh.basis(j).node(1);
  end
  if forceMesh.basis(j).node(2)<ymin_b
    ymin_b=forceMesh.basis(j).node(2);
  end
  if forceMesh.basis(j).node(2)>ymax_b
    ymax_b=forceMesh.basis(j).node(2);
  end
end


% here actually each basis function hast to be evaluated twice, which means
% it actually give two values:
ux=zeros(length(x_vec_u),2*numBasis);
uy=zeros(length(y_vec_u),2*numBasis);

% calculate the basis solution for all basis classes:
numClass=length(forceMesh.basisClass);
for class=1:numClass
    % To make sure that the range over which the solution is calculated,
    % take the double of the initial x and y ranges:
    xmin=min(min(x_vec_u),xmin_b); xmax=max(max(x_vec_u),xmax_b);
    ymin=min(min(y_vec_u),ymin_b); ymax=max(max(y_vec_u),ymax_b);
    dx=xmax-xmin;
    dy=ymax-ymin;
    
    xrange=[-dx dx]';
    yrange=[-dy dy]';
    
    % first basis solution:
    [ux1, uy1, ~, x_grid1, y_grid1]=fwdSolution3D(xrange ,yrange, z , E, v, [] , [] , [] , [] , forceMesh.basisClass(class).basisFunc(1).f_intp_x, forceMesh.basisClass(class).basisFunc(1).f_intp_y, [], 'fft', 'noIntp');

    % These field are good for intp2('*cubic'):
    forceMesh.basisClass(class).basisFunc(1).sol.x  = x_grid1;
    forceMesh.basisClass(class).basisFunc(1).sol.y  = y_grid1;
    forceMesh.basisClass(class).basisFunc(1).sol.ux = ux1;
    forceMesh.basisClass(class).basisFunc(1).sol.uy = uy1;

    
    % second basis solution:
    [ux2, uy2, ~, x_grid2, y_grid2]=fwdSolution3D(xrange ,yrange, z , E, v, [] , [] , [] , [] , forceMesh.basisClass(class).basisFunc(2).f_intp_x, forceMesh.basisClass(class).basisFunc(2).f_intp_y, [], 'fft', 'noIntp');
    
    % These field are good for intp2('*cubic'):
    forceMesh.basisClass(class).basisFunc(2).sol.x  = x_grid2;
    forceMesh.basisClass(class).basisFunc(2).sol.y  = y_grid2;
    forceMesh.basisClass(class).basisFunc(2).sol.ux = ux2;
    forceMesh.basisClass(class).basisFunc(2).sol.uy = uy2;
end


for j=1:numBasis
    display([num2str(j),' of: ',num2str(numBasis)]);
    % limits for integration: integrate basis function only over their
    % respective support:
    
    % Interpolate the basis-solution:
    class  = forceMesh.basis(j).class;
    xShift = forceMesh.basis(j).node(1);
    yShift = forceMesh.basis(j).node(2);
    
    x1  = forceMesh.basisClass(class).basisFunc(1).sol.x;
    y1  = forceMesh.basisClass(class).basisFunc(1).sol.y;
    ux1 = forceMesh.basisClass(class).basisFunc(1).sol.ux;
    uy1 = forceMesh.basisClass(class).basisFunc(1).sol.uy;
    
    x2  = forceMesh.basisClass(class).basisFunc(2).sol.x;
    y2  = forceMesh.basisClass(class).basisFunc(2).sol.y;
    ux2 = forceMesh.basisClass(class).basisFunc(2).sol.ux;
    uy2 = forceMesh.basisClass(class).basisFunc(2).sol.uy;
    
    % Then the interpolants of the first function are:    
    ux(:,j) = interp2(x1+xShift, y1+yShift, ux1, x_vec_u, y_vec_u,'*cubic');
    uy(:,j) = interp2(x1+xShift, y1+yShift, uy1, x_vec_u, y_vec_u,'*cubic');  
    
    % Then the interpolants of the second function are:  
    ux(:,j+numBasis) = interp2(x2+xShift, y2+yShift, ux2, x_vec_u, y_vec_u, '*cubic');
    uy(:,j+numBasis) = interp2(x2+xShift, y2+yShift, uy2, x_vec_u, y_vec_u, '*cubic');
end

M=vertcat(ux,uy);

% plot an example to see if it works correctly
ind=1;
if forceMesh.numBasis>ind-1
    xmin=min(x_vec_u);
    ymin=min(y_vec_u);
    xmax=max(x_vec_u);
    ymax=max(y_vec_u);
    figure(11)
    quiver(x_vec_u,y_vec_u,ux(:,ind),uy(:,ind))
    hold on
    quiver(x_vec_u,y_vec_u,ux(:,ind+forceMesh.numBasis),uy(:,ind+forceMesh.numBasis))
    xlim([xmin xmax])
    ylim([ymin ymax])
    hold off
end