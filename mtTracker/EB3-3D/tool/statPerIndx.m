function [means,meds,stds,orderedIndex,counts,sums] = statPerIndx(values,index)
means=zeros(1,max(index));
meds=zeros(1,max(index));
stds=zeros(1,max(index));
counts=zeros(1,max(index));
sums=zeros(1,max(index));

orderedIndex=1:max(index);
for t=orderedIndex
    if(t>0)
    biasAtTime=(values(index==t));
    means(t)=nanmean(biasAtTime);
    meds(t)=nanmedian(biasAtTime);
    stds(t)=nanstd(biasAtTime);
    counts(t)=nansum(index==t);
    sums(t)=nansum(biasAtTime);
    end
end
