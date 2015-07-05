function [ output_args ] = GCAVisualizeBranchPointsDistance(filoInfo )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

xyDist = vertcat(filoInfo(:).conXYCoords); 

scatter(xyDist(:,1),xyDist(:,2),5,'r','filled'); 
text(xyDist(:,1),xyDist(:,2),num2str(xyDist(:,3).*0.216,2)); 

end

