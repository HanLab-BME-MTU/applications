function [M]=refineFwdMapFastBEM(x_vec_u, y_vec_u, M_old, refinedForceMesh, E, meshPtsFwdSol, doPlot)
% this isn't really needed anymore. calcFwdMapFastBEM should be generalized
% that it can handle this case. Or, since the fwdCalculation is now so fast,
% just recalculate the whole thing again!

if nargin < 6 || isempty(meshPtsFwdSol)
    meshPtsFwdSol=[];
end

if nargin < 7 || isempty(doPlot)
    doPlot=0;
end

% The structure of M is:

%       sol. for fx                sol. for fy
%    ( ux(:,1:numBasis)   ux(:,numBasis+1:2*numBasis) )
% M= (                                                )
%    ( uy(:,1:numBasis)   uy(:,numBasis+1:2*numBasis) )


% get the original number of basis functions which is colsM/2:
[rowsM, colsM]=size(M_old);
numPts  =rowsM/2;
numBasis=colsM/2;

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

% Only calculate the basis solution for the basis classes that are needed.
% Find these classes: 
basisToBeUpdated = refinedForceMesh.IDsOfBasisSolToBeUpdated;
classToBeCalc    = vertcat(refinedForceMesh.basis(basisToBeUpdated).class);
% sort out doubled values:
classToBeCalcUnique = unique(classToBeCalc);
  
q=0;
numClassToBeCalc=length(classToBeCalcUnique);
for class=classToBeCalcUnique'  
    q=q+1;
    displayText=[num2str(class),' this is the: ',num2str(q),' of: ',num2str(numClassToBeCalc)];
    progressText(q/numClassToBeCalc,displayText)
    
    % First find all basis functions that are of this class and that we
    % have to re-calculate:
    basisToDoList=basisToBeUpdated(classToBeCalc==class);
    
    % To make sure that the range over which the solution is calculated,
    % take the double of the initial x and y ranges:
    xmin=min(x_vec_u); xmax=max(x_vec_u);
    ymin=min(y_vec_u); ymax=max(y_vec_u);
    dx=xmax-xmin;
    dy=ymax-ymin;
    
    xrange=[-dx dx]';
    yrange=[-dy dy]';
    
    for oneORtwo=1:2
        % calculate basis solution:
        [ux_model, uy_model, x_model, y_model]=fwdSolution(xrange,yrange,E,[],[],[],[],refinedForceMesh.basisClass(class).basisFunc(oneORtwo).f_intp_x,refinedForceMesh.basisClass(class).basisFunc(oneORtwo).f_intp_y,'fft','noIntp',meshPtsFwdSol);
        % The upper fields are good for intp2('*cubic')!        
        
        method='*cubic';
        for basisID=basisToDoList'            
            
            % Interpolate the basis-solution:
            xShift = refinedForceMesh.basis(basisID).node(1);
            yShift = refinedForceMesh.basis(basisID).node(2);
            
            if oneORtwo==1
                % Then the interpolants of the first function are:
                M(1:numPts    ,basisID)          = interp2(x_model+xShift, y_model+yShift, ux_model, x_vec_u, y_vec_u, method);  %This is ux(:,j)
                M(numPts+1:end,basisID)          = interp2(x_model+xShift, y_model+yShift, uy_model, x_vec_u, y_vec_u, method);  %This is uy(:,j)
            elseif oneORtwo==2
                % Then the interpolants of the second function are:  (:,j+numBasis)
                M(1:numPts    ,basisID+numBasis) = interp2(x_model+xShift, y_model+yShift, ux_model, x_vec_u, y_vec_u, method);  %This is ux(:,j+numBasis)
                M(numPts+1:end,basisID+numBasis) = interp2(x_model+xShift, y_model+yShift, uy_model, x_vec_u, y_vec_u, method);  %This is uy(:,j+numBasis)
            end
        end
    end
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