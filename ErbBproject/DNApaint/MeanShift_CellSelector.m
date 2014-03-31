function [CompStruct,CI]=MeanShift_CellSelector(clusterInfo)
%MeanShift_CellSelector takes the clusterInfo output of meanShift
%clustering analysis of multicolor paint imaging and prompts you to select
%the cell of interest and returns various information about each cluster
%formatted for machine learning
%

CI = CellSelector(clusterInfo);

comp = vertcat(CI.composition);
area = vertcat(CI.area);
den = comp./repmat(area,[1,size(comp,2)]);

CompStruct = struct('comp',comp,'den',den,'area',area,'name',[] );

end
