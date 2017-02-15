%% 
% LBP + deltaLBP
% take all cells with trajectories with frame time = t (frame)
function [feats, cellID] = getLchSpaceTimeT(cellData,t)
n = length(cellData);
feats = [];
cellID = [];
for i = 1 : n
    if sum(cellData{i}.ts==t)
        feats = [feats,[median(cellData{i}.lbp,1),median(cellData{i}.dLbp,1)]'];
        cellID = [cellID, i];
    end
end
end
