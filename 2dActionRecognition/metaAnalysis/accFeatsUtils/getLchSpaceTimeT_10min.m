%% 
% LBP + deltaLBP
% take all cells with trajectories with frame time = t (frame)
function [feats, cellID, TXY] = getLchSpaceTimeT_10min(cellData,t)
n = length(cellData);
feats = [];
cellID = [];
TXY = {};
ncells = 0;
for i = 1 : n
    if sum(cellData{i}.ts==t)
        ncells = ncells + 1;
        sind = find(cellData{i}.ts == t);
        eind = sind + 9;
        feats = [feats,[median(cellData{i}.lbp(sind:eind,:),1),median(cellData{i}.dLbp(sind:eind,:),1)]'];
        cellID = [cellID, i];        
        TXY{ncells}.ts = cellData{i}.ts(sind:eind);
        TXY{ncells}.xs = cellData{i}.xs(sind:eind);
        TXY{ncells}.ys = cellData{i}.ys(sind:eind);
    end
end
end
