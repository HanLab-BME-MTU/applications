%% Generic function for accumulating features (at the location level)
% fGetFeats - specific feature
% cellData - where to collect the data from
% Output: features of cells + their IDs
function [feats,ids] = getLchFeats(fGetFeats,cellData)
[feats,ids] = fGetFeats(cellData);
end