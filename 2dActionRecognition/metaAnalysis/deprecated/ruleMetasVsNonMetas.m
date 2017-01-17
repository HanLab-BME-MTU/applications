%% Classifies based on Metastatic potential
function [label, labelStr] = ruleMetasVsNonMetas(expName,metaData)

[expDetails] = parseExpName(expName,metaData);

if ~isfield(expDetails,'isHighMetPot')
    label = nan; labelStr = nan;
    return;
end

if expDetails.isHighMetPot 
    label = 2;
    labelStr = 'High metastatic potential';
else
    label = 1;
    labelStr = 'Low metastatic potential';
end

end