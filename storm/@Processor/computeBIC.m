function bic = computeBIC(logL,nParam,nPoints)
% function bic = computeBIC(logL,nParam,nPoints)
% SYNOPSIS:
% Computes the BIC.
%
% REQUIRED INPUTS:         
% - logL
% The joint log-likelihood.
% 
% - nParam
% The number of parameters of the model.
% 
% - nPoints
% The number of points represented by the model.
% 
% OPTIONAL INPUTS:
%
% NEEDED PROPERTIES: 
%
% MODIFIED PROPERTIES:
%
% OUTPUTS:
% - bic
% The value of the bic.
%
% Pascal Bérard, October 2011

bic = -2*logL+nParam*log(nPoints);

end

