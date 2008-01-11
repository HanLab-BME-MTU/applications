function [paramV,paramIdent,errFlag] = vectorFromParams(TOPO,nExoIn,...
        arPARAM,maPARAM,tryCONN,identify)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vectorFromParams takes a set of CARMA parameters and concatenates them
% into a single row vector for use by TOMLAB
%%%%% INPUT DESCRIPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       TOPOLOGY  : 3D Concatenation of adjacency matrices at each time lag 
%                   for connections between ARMA processes        
%   
%                   Structure:
%                        *Index in 1st dimension is node connection
%                        originates from
%                        *Index in 2nd dimension is node connection goes
%                        to
%                        *Index in 3rd dimension is time lag for each
%                        connection parameter
%                        *Diagonals should therefore be zero as this is
%                        equivalent to an autoregressive coefficient.
%--------------------------------------------------------------------------
%       nExoIn   : Number of exogenous inputs.
%--------------------------------------------------------------------------
%       arPARAM   : Matrix of autoregressive coefficients.
%       maPARAM   : Matrix of moving average coefficients.
%
%                    ARMA parameter matrix structure:
%                         *1st dimension is time lag for parameter
%                         *2nd dimension is node number
%                         *Both are optional, however they
%                         should be the same number of columns - if
%                         some nodes have AR but not MA then use NaN in
%                         place of MA coef and vice versa.
%--------------------------------------------------------------------------
%       tryCONN   : Optional. A 2D binary adjacency matrix indicating 
%                   the connections from TOPO to include in the vector. 
%
%%%%% OUTPUT DESCRIPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       paramV    : Vector of CARMA parameters.
%
%                      Structure (order of parameters):
%                           AR Params node 1, lag 1
%                           AR Params node 1, lag 2
%                                   ...
%                           AR Params node 1, lag p
%                           MA Params node 1, lag 1
%                           MA Params node 1, lag 2
%                                   ...
%                           MA Params node 1, lag q
%                           AR Params node 2, lag 1
%                           AR Params node 2, lag 2
%                                   ...
%                           AR Params node 2, lag p
%                           MA Params node 2, lag 1
%                           MA Params node 2, lag 2
%                                   ...
%                           MA Params node 2, lag q 
%                                   ...
%                                   ...
%                           ARMA params node n, lag p,q
%                           X parameter, connection 1, lag 0
%                           X parameter, connection 1, lag 1
%                                   ...
%                           X parameter, connection 1, lag m
%                           X parameter, connection 2, lag 0
%                           X parameter, connection 2, lag 1
%                                   ...
%                           X parameter, connection 2, lag m
%                                   ...
%                                   ...
%                           X parameter, connection r, lag m
%
%--------------------------------------------------------------------------
%       paramIdent : Vector listing identity of each parameter in paramV.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Hunter Elliott %%%   Last updated: November, 2007  %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
% INPUT CHECKS %
%%%%%%%%%%%%%%%%
    
    
%check how many inputs are included and set defaults accordingly    
  
if nargin < 2
    disp('--vectorFromParams: Require at least topography and external connection # as input!')
    errFlag = 1;
    return
end

%Check if identities of parameters were requested
if nargin < 6
    identify = false;
end

%Check if connectivity matrix was input
if nargin < 5
    tryCONN = [];
end

%Check if ARMA parameters included
if nargin < 3
    arPARAM = [];
    maPARAM = [];
end
if nargin < 4
    maPARAM = [];
end

%determine number of nodes, external inputs and ARMAX orders
[arOrder,nARp] = size(arPARAM);
[maOrder,nMAp] = size(maPARAM);
[nTot,nTotChk,xOrder] = size(TOPO);
xOrder = xOrder - 1;
nNodes = nTot - nExoIn;

%check if binary connectivity was provided and if not use all connections
%with non-zero x parameters in topology (TOPO)
if isempty(tryCONN)
    %disp('vectorFromParams: WARNING: No connectivity provided, using all non-zero connections in topology.');
    [connFrom,connTo] = find(sum(abs(TOPO),3));
    nTotChk2 = nTot;
    nTotChk3 = nTot;
else
   [connFrom,connTo] = find(tryCONN); 
   [nTotChk2,nTotChk3] = size(tryCONN);
end

nConn = length(connFrom);

ARMAtmp = cat(1,arPARAM,maPARAM);
armaParamT = [];
paramIdentT.iNode = [];
paramIdentT.fNode = [];
paramIdentT.paramLag = [];
paramIdentT.type = [];
ar = 'AR';
ma = 'MA';
X =  'XX';

if ~isempty(arPARAM) && ~isempty(maPARAM)
    if nARp ~= nNodes || nMAp ~= nNodes
        disp('vectorFromParams: incorrect number of ARMA params!');
        errFlag = 1;
        return
    end
end

%verify node number agreement in input parameters
if  (nTotChk2 ~= nTot) || (nTotChk3 ~= nTot)
    disp('Node # disagreement in parameters!');
    return
end

%if there are ARMA 
%if ~isempty(arPARAM) && ~isempty(maPARAM)
    for j = 1:nNodes
        armaParamT = cat(1, armaParamT, ARMAtmp(:,j));
        paramIdentT.iNode = cat(1, paramIdentT.iNode, ones(arOrder+maOrder,1) .* j); 
        paramIdentT.paramLag = cat(1, paramIdentT.paramLag, [1:arOrder 1:maOrder]');
        paramIdentT.type = cat(1, paramIdentT.type, repmat(ar,arOrder,1), repmat(ma,maOrder,1));
    end
%end %if arma

xParamT = [];
paramIdentT.fNode = nan(nNodes*(arOrder+maOrder),1);

for k = 1:nConn
    xParamT = cat(1, xParamT, squeeze(TOPO(connFrom(k),connTo(k),:)));
    paramIdentT.iNode = cat(1,paramIdentT.iNode,ones(xOrder+1,1) .* (connTo(k)-nExoIn));
    paramIdentT.fNode = cat(1,paramIdentT.fNode,ones(xOrder+1,1) .* (connFrom(k)-nExoIn));
    paramIdentT.paramLag = cat(1,paramIdentT.paramLag, [0:xOrder]');
    paramIdentT.type = cat(1, paramIdentT.type, repmat(X,xOrder+1,1));
end

paramV = cat(1,armaParamT,xParamT);

if identify
    paramIdent = paramIdentT;
end