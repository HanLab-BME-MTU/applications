function [prefixArray] = whToPrefix(strLabels,first)
if nargin < 2
    first = true;
end
n = length(strLabels);
prefixArray = cell(1,n);
for i = 1 : n
    tmp = strsplit(strLabels{i},'_');
    if first
        prefixArray{i} = tmp{1};
    else
        if length(tmp) == 2
            prefixArray{i} = tmp{2};
        else
            prefixArray{i} = '';
        end
    end
end
end