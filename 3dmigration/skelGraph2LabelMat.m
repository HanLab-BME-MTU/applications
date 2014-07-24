function labelMat = skelGraph2LabelMat(edgePaths,imSize)
%SKELGRAPH2LABELMAT returns a labelled 3D skeleton matrix (indexed image) corresponding to the input skeleton graph.
%
% labelMat = skelGraph2LabelMat(edgePaths,imSize,pixAspect)
%
% Returns a 3D matrix containing a labelled skeleton where each skeleton
% voxel is is set to an integer corresponding to that skeleton element's
% index in the input edgePaths array.
%
%   Input:
%
%       edgePaths - cell array with XYZ coordinates of each edge in the
%       skeleton, as output by skel2Graph. 
%
%       imSize - The image size which corresponds to the input edgePaths.
%       If the edge paths have been scaled to have symmetric voxels, this
%       should be the image size after scaling. 



% Hunter Elliott, 4/2013

labelMat = zeros(imSize,'uint16'); %I'm assuming we'll never have more than 65k skeletal elements..........


nEdge = numel(edgePaths);

%In case this is the scaled image and edgePaths
edgePaths = cellfun(@round,edgePaths,'Unif',0);

edgeInt = cellfun(@(x)(sub2ind(imSize,x(:,1),x(:,2),x(:,3))),edgePaths,'Unif',0);

for j = 1:nEdge
    labelMat(edgeInt{j}) = j;
end
    
