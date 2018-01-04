function [featsNorm] = LCH_setNormFeats(feats,meanVals,stdVals)
meanMat = repmat(meanVals,[1,size(feats,2)]);
stdMat = repmat(stdVals,[1,size(feats,2)]);
featsNorm = (feats - meanMat)./ stdMat;
end