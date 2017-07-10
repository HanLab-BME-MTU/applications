function [means,meds,stds,orderedIndex,counts] = statPerIndx(values,index)
means=zeros(1,max(index));
meds=zeros(1,max(index));
stds=zeros(1,max(index));
counts=zeros(1,max(index));

orderedIndex=1:max(index);
for t=unique(index)'
    biasAtTime=(values(index==t));
    means(t)=mean(biasAtTime);
    meds(t)=median(biasAtTime);
    stds(t)=std(biasAtTime);
    counts(t)=sum(index==t);
end
