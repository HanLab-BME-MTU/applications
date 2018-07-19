function plotData = slist2plotData(slist)
% this is a quick hack to transform the slist/cordList into plotData for
% visualization with imarisPlot3. A proper display function would take a
% LOT more effort than I can manage on a Sunday night at 1 am.

% initialization is not fun here, so I'm just glad it's a hack
plotData = struct('XYZ',[],'time',[]);
for t = 1:length(slist)
    newCoords = catStruct(1,'slist.sp.cord');
    nNewCoords = size(newCoords,1);
    plotData.XYZ = [plotData.XYZ;...
        newCoords];
    plotData.time = repmat(t,nNewCoords,1);
end
   
    