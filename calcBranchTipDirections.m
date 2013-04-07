function [branchDirVec,branchDirAng] = calcBranchTipDirections(edges,edgePaths,nVerts)


%Find only those which are tips
[~,iTipEdge] = findTips(edges,nVerts);

if size(edges,1) > 1
    nTips = numel(iTipEdge);
else
    %Handle the case where there's only one edge. This is bad but will be
    %handled elsewhere.
    nTips = 1;
end
branchDirVec = nan(nTips,3);

for j = 1:nTips
   
    %Get the branch direction vector
    if size(edgePaths{j},1) > 1
        branchDirVec(j,:) = curveDirection3D(edgePaths{j});
    end
            
end

%Convert these to polar coordinates
[branchDirAng(:,1),branchDirAng(:,2)] = cart2sph(branchDirVec(:,1),branchDirVec(:,2),branchDirVec(:,3));    


    



