function result=getci(y1Aligned)
[nWindow nTime]=size(y1Aligned);
for k=1:nTime
    l=1;
    for j=1:nWindow
        if isnan(y1Aligned(j,k))==0
            y1nonan(l)=y1Aligned(j,k);
            l=l+1;
        end
    end     
    if l<3
        ci_y1Aligned(:,k)=[NaN; NaN];
        mean_y1Aligned(k)=NaN;
     else
        ci_y1Aligned(:,k)=bootci(1000,@mean,y1nonan);
        mean_y1Aligned(k)=mean(y1nonan);
    end
    clear y1nonan;
end
result.mean=mean_y1Aligned;
result.ci=ci_y1Aligned;


    