function [feats,inds] = LCH_feats2Matrix(cellFeats)
nCells = length(cellFeats);
tmp = cellFeats{1};
nFeats = length(tmp);

inds = 1:nCells;

feats = nan(nFeats,nCells);

for icell = 1 : nCells
    curFeats = cellFeats{icell};    
    feats(:,icell) = curFeats;
end

end
