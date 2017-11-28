function [means,meds,stds,orderedIndex,counts,sums] = statPerIndx(values,index,indexValues)
means=zeros(1,max(index));
meds=zeros(1,max(index));
stds=zeros(1,max(index));
counts=zeros(1,max(index));
sums=zeros(1,max(index));

if(nargin<3)
	indexValues=1:max(index);
end
orderedIndex=indexValues;
for t=indexValues
    if(t>0)
    biasAtTime=(values(index==t));
    means(t)=mean(biasAtTime);
    meds(t)=median(biasAtTime);
    stds(t)=std(biasAtTime);
    counts(t)=sum(index==t);
    sums(t)=sum(biasAtTime);
    end
end