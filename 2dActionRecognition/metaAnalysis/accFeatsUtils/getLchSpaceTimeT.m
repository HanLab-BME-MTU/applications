%% 
% LBP + deltaLBP
% take all cells with trajectories with frame time = t (frame)
function [feats, cellID, TXY] = getLchSpaceTimeT(cellData,t)
n = length(cellData);
feats = [];
cellID = [];
TXY = {};
ncells = 0;
for i = 1 : n
    if sum(cellData{i}.ts==t)
        ncells = ncells + 1;
        feats = [feats,[median(cellData{i}.lbp,1),median(cellData{i}.dLbp,1)]'];
        cellID = [cellID, i];        
        TXY{ncells}.ts = cellData{i}.ts;
        TXY{ncells}.xs = cellData{i}.xs;
        TXY{ncells}.ys = cellData{i}.ys;
    end
end
end
