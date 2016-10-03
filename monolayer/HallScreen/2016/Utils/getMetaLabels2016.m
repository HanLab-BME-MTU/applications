% ??
function [labels, strLabels] = getMetaLabels2016(metaData)
nLabels = length(metaData.groupsByTreatments);
nExps = metaData.N;

labels = zeros(1,nExps);
strLabels = cell(1,nExps);

for curLabel = 1 : nLabels
    labels(metaData.groupsByTreatments{curLabel}.inds) = curLabel;
    strLabels(metaData.groupsByTreatments{curLabel}.inds) = metaData.groupsByTreatments{curLabel}.treatment;
end


end