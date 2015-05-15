function [ output_args ] = GCAVisualsPlotFilopodiaPerBranchGroup(filoInfo,imgSize)
%GCAVisualsPlotFilopodiaPerBranchGroup

groupLabels = unique(vertcat(filoInfo.groupCount),'stable');
n = length(groupLabels);
c = linspecer(n);
% c = colormap(jet(n));
reOrder = randperm(n);
c = c(reOrder,:);



hold on
% for each group plot the xy coords of all the filo in that group
for iGroup = 1:length(groupLabels)
    groupLabelsAll = vertcat(filoInfo.groupCount);
    clusterC =  filoInfo(groupLabelsAll == iGroup);
    GCAVisualsMakeOverlaysFilopodia(clusterC,imgSize,1,1,c(iGroup,:),0);
    %arrayfun(@(x) plot(x.Ext_coordsXY(:,1),x.Ext_coordsXY(:,2),'color', c(iGroup,:),'Linewidth',2),clusterC);
    % clusterInt = clusterC;
    %clusterInt = clusterInt(arrayfun(@(x) ~isnan(x.Int_coordsXY(1)),clusterC));
    %if ~isempty(clusterInt)
    %arrayfun(@(x) plot(x.Int_coordsXY(:,1),x.Int_coordsXY(:,2),'color',c(iGroup,:),'Linewidth',2),clusterInt);
    %end
end

