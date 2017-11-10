%% 
% LBP + deltaLBP
% take all cells with trajectories at t = 5 (frame)
function feats = getLchSpaceTime5(cellData)
feats = getLchSpaceTimeT(cellData,5);
end