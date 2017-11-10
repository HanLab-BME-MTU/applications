function [IVals SVals SMagVals Pos]=pairSIprofiles(SPosIn,SValsIn,IPosIn,IValsIn)
% IPosIn=groupedNetworks.cluster{1}.trackedNet{1}.edge{1}.intf_internal;
% SPosIn=groupedNetworks.cluster{1}.trackedNet{1}.edge{1}.cntrs;
% IValsIn=groupedNetworks.cluster{1}.trackedNet{1}.edge{1}.int.val;
% we have to take the stress here!
% SValsIn=groupedNetworks.cluster{1}.trackedNet{1}.edge{1}.s_vec;

% calculate the distance Matrix:
distMat=createDistanceMatrix(SPosIn,IPosIn);

% find the nearest force node for each position along the intensity curve:
[~,fPtId]=min(distMat,[],1);

% sort out multiple Ids (the result is ordered!):
uniqueIds=unique(fPtId);
numPts=length(uniqueIds);

%initialize:
IVals   =zeros(numPts,1);
SVals   =zeros(numPts,2);
SMagVals=zeros(numPts,1);
Pos     =zeros(numPts,2);

k=1;
for iPt=uniqueIds
    checkVec   =(fPtId==iPt);
    IVals(k)   =mean(IValsIn(checkVec));
    SVals(k,:) =SValsIn(iPt,:);
    SMagVals(k)=norm(SValsIn(iPt,:));
    Pos(k,:)   =SPosIn(iPt ,:);
    k=k+1;
end

