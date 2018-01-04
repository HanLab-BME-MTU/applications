function [annLabels,indsLabels] = LCH_getAnnLabels(annotationInfo,curLabel)

n = length(annotationInfo);
indsLabels = [];
annLabels = cell(1,n);


for i = 1 : n
    curAnnots = annotationInfo{i}.annotation;
    for j = 1 : length(curAnnots)
        if strcmp(curAnnots{j},curLabel)
            annLabels{i} = curLabel;
        end
    end
    if length(curAnnots) == 0
        continue;
    else
        if isempty(annLabels{i})
            annLabels{i} = 'Rest';
        end
    end
    indsLabels = [indsLabels i];
end
end