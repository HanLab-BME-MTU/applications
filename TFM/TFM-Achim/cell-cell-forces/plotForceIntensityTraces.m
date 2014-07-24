function []=plotForceIntensityTraces(groupedNetworks)
doPrint=1;
if nargin<1 || isempty(groupedNetworks)
    try
        load('groupedNetworks.mat')
    catch
        display('Couldn''t find any groupedNetworks... and stoped!')
        return;
    end    
end

relErrF_val_corr=Inf;
goodEdgeSet=findEdges(groupedNetworks,'myo',[0],'relErrF',relErrF_val_corr,'errs',0);
% goodEdgeSet=findEdges(groupedNetworks,'myo',[0],'asmbly',[1],'relErrF',relErrF_val_corr,'errs',0);
if ~isempty(goodEdgeSet) && ~isempty(goodEdgeSet(1).edgeId)
    [corrSets]=collectEdgeValues(groupedNetworks,goodEdgeSet,'corr');
else
    display('No myosin cells of this type found!')
end

%**************************************************************************
% Plot the edges from above:
%**************************************************************************
plotForceIntensityTracesCorrSet(corrSets,doPrint)