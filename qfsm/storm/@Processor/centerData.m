function centerData(obj)
% function centerData(obj)
% SYNOPSIS:
% Centers the data around [0,0,0]. 
%
% REQUIRED INPUTS:         
% 
% OPTIONAL INPUTS:
%
% NEEDED PROPERTIES: 
% - obj.data.points
%
% MODIFIED PROPERTIES: 
% - obj.data.points
%
% OUTPUTS:
%
% Pascal Bérard, October 2011

% Compute the center of the point cloud
center = (max(obj.data.points,[],1)+min(obj.data.points,[],1))/2;

% Translate all the points
obj.data.points = obj.data.points+repmat(-center,size(obj.data.points,1),1);

disp('Process: Data centered!');

end