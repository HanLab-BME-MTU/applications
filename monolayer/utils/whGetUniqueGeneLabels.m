% also gives the unique sequences if "first" is equal to "false"
% assumes the label is of the format gene_sequence
function [uniqueGenes] = whGetUniqueGeneLabels(strLabels,first)
if nargin < 2
    first = true;
end
prefixArray = whToPrefix(strLabels,first);
uniqueGenes = {};
n = 1;
while ~isempty(prefixArray)    
    uniqueGenes{n} = prefixArray{1};
    prefixArray(strcmp(prefixArray,prefixArray{1})) = '';
    n = n + 1;
end
end