function [res,jac]=distTestError(parms,data,gIdxList,mskDataSize,dataProperties)
% callback for nonlinear fitter in function 'fitTest'

% INPUT  parms          : list of gaussian parms and background value
%             data             : vector with intensity data
%             gIdxList        : indices to valid entries (mask)
%             mskDataSize: vector with 3D size of mask 
%
% OUTPUT res  : residual of model
%                jac  : jacobian of model

if(nargin==2)
    dataSize=size(data);
    idxList=[];
%    mask=ones(dataSize);
end;
% create model with current parms
[gaussFit, mgrad] = multiGaussFit(mskDataSize,parms,dataProperties);
% return residual and jacobian
res=data-gaussFit(gIdxList);
jac=-mgrad(gIdxList,:);
