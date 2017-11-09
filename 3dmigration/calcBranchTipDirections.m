function [branchDirVec,branchDirAng] = calcBranchTipDirections(edges,edgePaths,nVerts)


%Find only those which are tips
[~,iTipEdge] = findTips(edges,nVerts);

if size(edges,1) > 1
    nTips = numel(iTipEdge);
    branchDirVec = nan(nTips,3);

    for j = 1:nTips
       
        %Get the branch direction vector
        if size(edgePaths{iTipEdge(j)},1) > 1
            %Tip edges start at tip, so we reverse these to get direction
            %branch is "pointing" - from base to tip. They are also in YXZ
            %matrix coord so we swap dim before passing as well
            branchDirVec(j,:) = curveDirection3D(edgePaths{iTipEdge(j)}(end:-1:1,[2 1 3]));

        end
            
    end
    
else
    %Handle the case where there's only one edge. This is bad but will be
    %handled elsewhere.    
    
    %In this special case we use both ends of the one edge
    branchDirVec(1,:) = curveDirection3D(edgePaths{1}(end:-1:1,[2 1 3]));
    branchDirVec(2,:) = curveDirection3D(edgePaths{1}(:,[2 1 3]));
    
end





%Convert these to polar coordinates
[branchDirAng(:,1),branchDirAng(:,2)] = cart2sph(branchDirVec(:,1),branchDirVec(:,2),branchDirVec(:,3));    


    



