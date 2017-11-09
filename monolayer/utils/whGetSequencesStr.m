function [uniqueSeqStr, uniqueSeqInds] = whGetSequencesStr(strLabels,geneStr)
indsGene = strcmp(whToPrefix(strLabels),geneStr);
interestLabels = strLabels(indsGene);
uniqueSeqStr =  whGetUniqueGeneLabels(interestLabels,false);
N = length(uniqueSeqStr);
uniqueSeqInds = cell(1,N);
for i = 1 : N
    uniqueSeqInds{i} = strcmp(strLabels,[geneStr '_' uniqueSeqStr{i}]);
    
    if sum(uniqueSeqInds{i}) == 0
        uniqueSeqInds{i} = strcmp(strLabels,geneStr);
    end
end
end