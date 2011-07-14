function mergedField=mergeDisplFields(basicField,extraField)
mergedField=basicField;

for iframe=1:min(length(basicField),length(extraField))
    % find the bounds of the field:
    allPos=vertcat(basicField(iframe).pos,extraField(iframe).pos);
    [maxX,maxY]=max(allPos);

    maxX=maxX+1;
    maxY=maxY+1;

    [basicInd]=sub2ind([maxX maxY],basicField(iframe).pos(:,1),basicField(iframe).pos(:,2));
    [extraInd]=sub2ind([maxX maxY],extraField(iframe).pos(:,1),extraField(iframe).pos(:,2));

    addInd=setdiff(extraInd,basicInd);

    mergedField(iframe).pos=vertcat(basicField(iframe).pos,extraField(iframe).pos(addInd,:));
    mergedField(iframe).vec=vertcat(basicField(iframe).vec,extraField(iframe).vec(addInd,:));
end