%% Classifies based on a specific job and tumor
function [label, labelStr] = ruleJobTumor(expName,metaData)

[expDetails] = parseExpName(expName,metaData);

if ~isfield(expDetails,'isHighMetPot') || isempty(expDetails.isHighMetPot)
    label = nan; labelStr = nan;
    return;
end

label = expDetails.jobTumorPairInd;

labelStr = expDetails.jobTumorPair;

end