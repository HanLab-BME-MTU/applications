%% 
% LBP + deltaLBP
% take all cells with trajectories with frame time = t (frame)
function feats = getLchSpaceTimeT(cellData,t)
n = length(cellData);
feats = [];
for i = 1 : n
    if sum(cellData{i}.ts==t)
        feats = [feats,[median(cellData{i}.lbp,1),median(cellData{i}.dLbp,1)]'];
    end
end
end
