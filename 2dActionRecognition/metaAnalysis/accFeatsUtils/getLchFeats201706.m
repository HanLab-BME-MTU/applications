%% Generic function for accumulating features (at the location level)
% fGetFeats - specific feature
% featsInFname - where to collect the data from
% Output: features of cells + their IDs
function [feats,ids,txy] = getLchFeats201706(fGetFeats,featsInFname)
[feats,ids,txy] = fGetFeats(featsInFname);
end