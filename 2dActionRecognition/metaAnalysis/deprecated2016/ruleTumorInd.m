%% Classifies based on Tumor index
function [label, labelStr] = ruleTumorInd(expName,metaData)

[expDetails] = parseExpName(expName,metaData);

if ~isfield(expDetails,'isHighMetPot') || isempty(expDetails.isHighMetPot)
    label = nan; labelStr = nan;
    return;
end

label = expDetails.tumorInd;

labelStr = expDetails.tumor;

end