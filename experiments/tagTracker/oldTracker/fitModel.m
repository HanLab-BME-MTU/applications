function [res,jac]=fitModel(parms,data,gIdxList,mskDataSize)
% callback for nonlinear fitter in function 'trackSpots'

% INPUT  parms          : list of model parms
%             data             : vector with relevant raw data
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
[gaussFit, mgrad] = mapData(mskDataSize,parms);
% return residual and jacobian
res=data-gaussFit(gIdxList);
jac=-mgrad(gIdxList,:);
