function [ output_args ] = GCAVisualsPlotFilopodiaIDs( filoInfo,numbers,color )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
IDs = 1:length(filoInfo); 
filoInfoToPlot = filoInfo(numbers); 
IDsToPlot = IDs(numbers); 


arrayfun(@(i,j) text(i.Ext_coordsXY(1,1),i.Ext_coordsXY(1,2),num2str(j),'color',color),filoInfoToPlot,IDsToPlot);


end

