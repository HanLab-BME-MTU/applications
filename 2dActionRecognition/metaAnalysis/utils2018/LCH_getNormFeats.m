function [meanVals,stdVals] = LCH_getNormFeats(featsAll)
meanVals = mean(featsAll,2);
stdVals = std(featsAll,0,2);
end