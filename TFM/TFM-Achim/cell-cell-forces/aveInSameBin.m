function vec=aveInSameBin(vec,binList,aveType)
for bin=min(binList):max(binList)
    binGrp=find(binList==bin);
    if strcmpi(aveType,'mean') 
        vec(binGrp,:)=repmat(mean(vec(binGrp,:),1),length(binGrp),1);
    elseif strcmpi(aveType,'nanmean')
        vec(binGrp,:)=repmat(nanmean(vec(binGrp,:),1),length(binGrp),1);
    end
end