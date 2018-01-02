%% LCH_orginizeFeatures2018
% Orginize in conivnient data structures + normalize
function [] = LCH_orginizeFeatures2018(featsDname,featsStrID)
featsFname = [featsDname filesep featsStrID '_all.mat'];
load(featsFname); % allCells

outdir = [featsDname filesep featsStrID];
if ~exist(outdir,'dir')
    mkdir(outdir);
end

nCellTypeDate = length(allCells.strs);

end