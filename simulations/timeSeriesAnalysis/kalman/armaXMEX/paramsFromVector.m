function [arPARAMS,maPARAMS,TOPO] = paramsFromVector(paramV,nNodes,arOrder,...
    maOrder,xOrder,tryCONN)

if nargin < 6
    disp('Must include all input arguments!');
    return
end

[nTot,nTotChk] = size(tryCONN);
[connFrom,connTo] = find(tryCONN);
nConn = length(connFrom);

arPARAMS = zeros(arOrder,nNodes);
maPARAMS = zeros(maOrder,nNodes);
TOPO = zeros(nTot,nTot,xOrder+1);
armaOrder = arOrder+maOrder;

armaPARAMS = [];
for j = 1:nNodes
    armaPARAMS = cat(2,armaPARAMS,paramV( (j-1)*armaOrder+1 : j*armaOrder ));
end

arPARAMS = armaPARAMS(1:arOrder,:);
maPARAMS = armaPARAMS(arOrder+1:arOrder+maOrder,:);

for k = 1:nConn
    TOPO(connFrom(k),connTo(k),:) = paramV(nNodes*(arOrder+maOrder) + (k-1)*(xOrder+1) + 1 : nNodes*(arOrder+maOrder) + k*(xOrder+1));
end