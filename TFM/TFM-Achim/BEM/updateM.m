function [M]=updateM(M,forceMesh,basisID,oneORtwo,x_model,y_model,ux_model,uy_model,x_data,y_data,method)
% Get the dimensions of M:
[rowsM, colsM]=size(M);
numPts  =rowsM/2;
numBasis=colsM/2;

% Interpolate the basis-solution:
xShift = forceMesh.basis(basisID).node(1);
yShift = forceMesh.basis(basisID).node(2);

if oneORtwo==1
    % Then the interpolants of the first function are:
    M(1:numPts    ,basisID)             = interp2(x_model+xShift, y_model+yShift, ux_model, x_data, y_data,method);  %This is ux(:,j)
    M(numPts+1:end,basisID)             = interp2(x_model+xShift, y_model+yShift, uy_model, x_data, y_data,method);  %This is uy(:,j)
elseif oneORtwo==2
    % Then the interpolants of the second function are:  (:,j+numBasis)
    M(1:numPts    ,basisID+numBasis) = interp2(x_model+xShift, y_model+yShift, ux_model, x_data, y_data,method);  %This is ux(:,j+numBasis)
    M(numPts+1:end,basisID+numBasis) = interp2(x_model+xShift, y_model+yShift, uy_model, x_data, y_data,method);  %This is uy(:,j+numBasis)
end