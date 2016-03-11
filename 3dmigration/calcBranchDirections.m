function [branchDirVec,branchDirAng] = calcBranchDirections(edges,edgePaths,nVerts,edgeDistances)

nBranch = numel(edgePaths);

branchDirVec = nan(nBranch,3);

for j = 1:nBranch
       
    %Get the branch direction vector
    if size(edgePaths{j},1) > 1 && size(edgeDistances{j},1) > 1
        
        %Check the orientation direction of this branch relative to cell
        %center and flip so we get direction outward from cell center        
        if edgeDistances{j}(1) < edgeDistances{j}(end)
            currPath = edgePaths{j};
        else
            currPath = edgePaths{j}(end:-1:1,:);
        end
        
        branchDirVec(j,:) = curveDirection3D(currPath(:,[2 1 3]));%Convert from matrix coord to XYZ before getting direction

    end
            
end

%Convert these to polar coordinates
[branchDirAng(:,1),branchDirAng(:,2)] = cart2sph(branchDirVec(:,1),branchDirVec(:,2),branchDirVec(:,3));    


    



