function [M]=refineFwdMapFastBEM_old(x_vec_u, y_vec_u, M_old, refinedForceMesh, E, doPlot)

if nargin < 6 || isempty(doPlot)
    doPlot=0;
end

% The structure of M is:

%       sol. for fx                sol. for fy
%    ( ux(:,1:numBasis)   ux(:,numBasis+1:2*numBasis) )
% M= (                                                )
%    ( uy(:,1:numBasis)   uy(:,numBasis+1:2*numBasis) )


% get the original number of basis functions which is colsM/2:
[rowsM, colsM]=size(M_old);
numPts=rowsM/2;

% In order to extend M take the max number of the basis functions:
oldNumBasis=colsM/2;
newNumBasis=refinedForceMesh.numBasis;
if newNumBasis>oldNumBasis
    extendMto=newNumBasis;
else
    % nothing has to be done.
    extendMto=[];
end

% The left half contains the displacements in the x-y-directions caused by
% the basis functions fx: 
M_old_fx=M_old(:,1:oldNumBasis);
%extend it to the new number of basis functions:
if ~isempty(extendMto)
    M_old_fx(rowsM,extendMto)=0;
end

% The right half contains the displacements in the  x-y-directions caused by
% the basis functions fy: 
M_old_fy=M_old(:,(oldNumBasis+1):end);
%extend it to the new number of basis functions:
if ~isempty(extendMto)
    M_old_fy(rowsM,extendMto)=0;
end

% stich the two pieces together again:
M=horzcat(M_old_fx,M_old_fy);

% now clear M to free up memory:
clear M_old M_old_fx M_old_fy;

% here actually each basis function hast to be evaluated twice, which means
% it actually give two values:
% ux=zeros(length(x_vec_u),2*numBasis);
% uy=zeros(length(y_vec_u),2*numBasis);

% calculate the basis solution for all basis classes. If those would be
% stored in the force mesh, then one would not need to recalculate the ones
% which are already known. However, refinement usually vastly increase the
% number of basis classes anyway. So recalculating e.g. 1 old vs. 20 new
% ones make the calculation only 5% slower.
numClass=length(refinedForceMesh.basisClass);
for class=1:numClass
    % To make sure that the range over which the solution is calculated,
    % take the double of the initial x and y ranges:
    xmin=min(x_vec_u); xmax=max(x_vec_u);
    ymin=min(y_vec_u); ymax=max(y_vec_u);
    dx=xmax-xmin;
    dy=ymax-ymin;
    
    xrange=[-dx dx]';
    yrange=[-dy dy]';
    
    % first basis solution:
    [ux1, uy1, x_grid1, y_grid1]=fwdSolution(xrange,yrange,E,[],[],[],[],refinedForceMesh.basisClass(class).basisFunc(1).f_intp_x,refinedForceMesh.basisClass(class).basisFunc(1).f_intp_y,'fft','noIntp');

    % These field are good for intp2('*cubic'):
    refinedForceMesh.basisClass(class).basisFunc(1).sol.x  = x_grid1;
    refinedForceMesh.basisClass(class).basisFunc(1).sol.y  = y_grid1;
    refinedForceMesh.basisClass(class).basisFunc(1).sol.ux = ux1;
    refinedForceMesh.basisClass(class).basisFunc(1).sol.uy = uy1;

    
    % second basis solution:
    [ux2, uy2, x_grid2, y_grid2]=fwdSolution(xrange,yrange,E,[],[],[],[],refinedForceMesh.basisClass(class).basisFunc(2).f_intp_x,refinedForceMesh.basisClass(class).basisFunc(2).f_intp_y,'fft','noIntp');
    
    % These field are good for intp2('*cubic'):
    refinedForceMesh.basisClass(class).basisFunc(2).sol.x  = x_grid2;
    refinedForceMesh.basisClass(class).basisFunc(2).sol.y  = y_grid2;
    refinedForceMesh.basisClass(class).basisFunc(2).sol.ux = ux2;
    refinedForceMesh.basisClass(class).basisFunc(2).sol.uy = uy2;
end


% It would make sense to now run through all new basis functions of this
% class, calculating the solution and then freeing up the memory! This
% should be done in the previous loop and not in an extra loop as follows
% here:
k=0;
toDo=length(refinedForceMesh.IDsOfBasisSolToBeUpdated);
for j=refinedForceMesh.IDsOfBasisSolToBeUpdated'
    k=k+1;
    displayText=[num2str(j),' this is the: ',num2str(k),' of: ',num2str(toDo)];
    progressText(k/toDo,displayText)
    % limits for integration: integrate basis function only over their
    % respective support:
    
    % Interpolate the basis-solution:
    class  = refinedForceMesh.basis(j).class;
    xShift = refinedForceMesh.basis(j).node(1);
    yShift = refinedForceMesh.basis(j).node(2);
    
    x1  = refinedForceMesh.basisClass(class).basisFunc(1).sol.x;
    y1  = refinedForceMesh.basisClass(class).basisFunc(1).sol.y;
    ux1 = refinedForceMesh.basisClass(class).basisFunc(1).sol.ux;
    uy1 = refinedForceMesh.basisClass(class).basisFunc(1).sol.uy;
    
    x2  = refinedForceMesh.basisClass(class).basisFunc(2).sol.x;
    y2  = refinedForceMesh.basisClass(class).basisFunc(2).sol.y;
    ux2 = refinedForceMesh.basisClass(class).basisFunc(2).sol.ux;
    uy2 = refinedForceMesh.basisClass(class).basisFunc(2).sol.uy;
    
    % Then the interpolants of the first function are:    
    M(1:numPts    ,j) = interp2(x1+xShift, y1+yShift, ux1, x_vec_u, y_vec_u,'*cubic');  %This is ux(:,j)
    M(numPts+1:end,j) = interp2(x1+xShift, y1+yShift, uy1, x_vec_u, y_vec_u,'*cubic');  %This is uy(:,j)
    
    % Then the interpolants of the second function are:  (:,j+numBasis)
    M(1:numPts    ,j+newNumBasis) = interp2(x2+xShift, y2+yShift, ux2, x_vec_u, y_vec_u, '*cubic');                %This is ux(:,j+numBasis)
    M(numPts+1:end,j+newNumBasis) = interp2(x2+xShift, y2+yShift, uy2, x_vec_u, y_vec_u, '*cubic');  %This is uy(:,j+numBasis)
end


% plot an example to see if it works correctly
% if doPlot==1
%     ind=1;
%     if refinedForceMesh.numBasis>ind-1
%         xmin=min(x_vec_u);
%         ymin=min(y_vec_u);
%         xmax=max(x_vec_u);
%         ymax=max(y_vec_u);
%         figure(11)
%         quiver(x_vec_u,y_vec_u,ux(:,ind),uy(:,ind))
%         hold on
%         quiver(x_vec_u,y_vec_u,ux(:,ind+refinedForceMesh.numBasis),uy(:,ind+refinedForceMesh.numBasis))
%         xlim([xmin xmax])
%         ylim([ymin ymax])
%         hold off
%     end
% end