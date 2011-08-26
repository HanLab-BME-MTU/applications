function [groupedNetworksSelected,listTaken,listDism]=preSelectForExp(groupedNetworks,expList)
% expList={'2010_11_18_TFM_shMYOIIA93_stiffness','2010_08_18_TFM_shMYOIIa','2010_12_19_TFM_shMYOIIA93'}
% The function runs through the clusterList and selects all clusters whose
% filename of the first bead image contains the pattern in expList

% These are all experiments used till 08/26/2011:

% talin-1 experiments:
% 2011_02_22_TFL_siTLN1_siTLN2

% myo-IIA experiments mixed with only control (only 8kPa?):
% 2010_08_18_TFM_shMYOIIa
% 2010_11_18_TFM_shMYOIIA93_stiffness % there is __, _1, _2, _3 % big exp.
% 2010_12_18_TFM_shMYOIIA93
% 2010_12_19_TFM_shMYOIIA93

% myo-IIB experiments (only 35kPa?):
% 2010_09_23_TFM_shMYOIIB103
% 2010_08_12_TFM_10AshMYOIIB
% 2010_09_23_TFM_103_94

% only control:
% 2010_08_01_TFM_10AEcadGFP_3vs35

%--------------------------------------------------------------------------
[~,numClusters]=size(groupedNetworks.clusterList);
cv=zeros(numClusters,1);
for iEntry=1:length(expList)
    for iCluster=1:numClusters
        foundPos=strfind(groupedNetworks.clusterList{iCluster}, expList{iEntry});
        if ~isempty(foundPos)
            cv(iCluster)=true;
        end            
    end
end

k=1;
d=1;
groupedNetworksSelected.numClusters=sum(cv);
for iCluster=1:numClusters
    if cv(iCluster)
        groupedNetworksSelected.clusterList{k}=groupedNetworks.clusterList{iCluster};
        groupedNetworksSelected.cluster{k}    =groupedNetworks.cluster{iCluster};
        k=k+1;
    else
        listDism{d}=groupedNetworks.clusterList{iCluster};
        d=d+1;
    end
end
listTaken=groupedNetworksSelected.clusterList;