function labelMat = labelMaskVoxelsBySkeleton(maskIn,skelLabel)
%LABELMASKVOXELSBYSKELETON labels each mask voxel according to the nearest skeleton edge
%
% labelMat = labelMaskVoxelsBySkeleton(maskIn,skelLabel)
%   
%   Assigns a label to each voxel in the input mask which associates it
%   with the closest element of the input labelled skeleton.
%
%   Input:
%       maskIn - the mask to be labelled, which the skeleton was derived
%       from.
%
%       skelLabel - The labelled skeleton matrix, e.g. as produced by
%       skelGraph2LabelMat.m
%
%

% Hunter Elliott, 4/2013


%Get distance transform and label matrix. We have to use bwdist_old because
%of the bug in the labelmatrix in bwdist.m
[~,distL] = bwdist_old(skelLabel > 0);
distL(~maskIn) = 0;

iSkelLab = unique(skelLabel(skelLabel(:)>0));
nSkelLab = numel(iSkelLab);

%Get the indices of the skeleton voxels for each element
distLforSkel = arrayfun(@(x)(distL(skelLabel == x)),iSkelLab,'Unif',0);

labelMat = zeros(size(maskIn),'uint16');

for j = 1:nSkelLab
    
    %For each skeleton element, find all the mask voxels for which this is
    %the closest element (as indicated by matching indices in the label
    %matrix) and label them accordingly
    labelMat(any(bsxfun(@eq,distL(:),distLforSkel{j}'),2)) = j;    
    
end