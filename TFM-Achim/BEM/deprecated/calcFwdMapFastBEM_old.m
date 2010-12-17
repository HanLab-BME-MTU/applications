function [M]=calcFwdMapFastBEM_old(x_vec_u, y_vec_u, forceMesh, E, doPlot)

if nargin < 5 || isempty(doPlot)
    doPlot=0;
end

numBasis=length(forceMesh.basis);

% here actually each basis function hast to be evaluated twice, which means
% it actually give two values:
ux=zeros(length(x_vec_u),2*numBasis);
uy=zeros(length(y_vec_u),2*numBasis);

% calculate the basis solution for all basis classes:
numClass=length(forceMesh.basisClass);
for class=1:numClass
    % To make sure that the range over which the solution is calculated,
    % take the double of the initial x and y ranges:
    xmin=min(x_vec_u); xmax=max(x_vec_u);
    ymin=min(y_vec_u); ymax=max(y_vec_u);
    dx=xmax-xmin;
    dy=ymax-ymin;
    
    xrange=[-dx dx]';
    yrange=[-dy dy]';
    
    % Integration bounds used in the refine step:
    xbd_min=min(forceMesh.basisClass(class).neighPos(:,1));
    xbd_max=max(forceMesh.basisClass(class).neighPos(:,1));
    ybd_min=min(forceMesh.basisClass(class).neighPos(:,2));
    ybd_max=max(forceMesh.basisClass(class).neighPos(:,2));
    
    % first basis solution:
    [ux1, uy1, x_grid1, y_grid1]=fwdSolution(xrange,yrange,E,xbd_min,xbd_max,ybd_min,ybd_max,forceMesh.basisClass(class).basisFunc(1).f_intp_x,forceMesh.basisClass(class).basisFunc(1).f_intp_y,'fft','noIntp');

    % These field are good for intp2('*cubic'):
    forceMesh.basisClass(class).basisFunc(1).sol.x  = x_grid1;
    forceMesh.basisClass(class).basisFunc(1).sol.y  = y_grid1;
    forceMesh.basisClass(class).basisFunc(1).sol.ux = ux1;
    forceMesh.basisClass(class).basisFunc(1).sol.uy = uy1;

    
    % second basis solution:
    [ux2, uy2, x_grid2, y_grid2]=fwdSolution(xrange,yrange,E,xbd_min,xbd_max,ybd_min,ybd_max,forceMesh.basisClass(class).basisFunc(2).f_intp_x,forceMesh.basisClass(class).basisFunc(2).f_intp_y,'fft','noIntp');
    
    % These field are good for intp2('*cubic'):
    forceMesh.basisClass(class).basisFunc(2).sol.x  = x_grid2;
    forceMesh.basisClass(class).basisFunc(2).sol.y  = y_grid2;
    forceMesh.basisClass(class).basisFunc(2).sol.ux = ux2;
    forceMesh.basisClass(class).basisFunc(2).sol.uy = uy2;
end

for j=1:numBasis
    displayText=[num2str(j),' of ',num2str(numBasis)];
    progressText(j/numBasis,displayText)
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
if doPlot==1
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
end