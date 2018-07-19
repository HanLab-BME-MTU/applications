function c=corrAnalysis(x,y,lag,minLen) 

[nWin, n]=size(x);

for j=1:nWin
    
    xLen=sum(~isnan(x(j,:)));
    yLen=sum(~isnan(y(j,:)));
    
    if xLen>yLen
        Len=yLen;
    else
        Len=xLen;
    end
    
    if Len>minLen
        if lag>(Len/2)
            thisLag=floor(Len/2);
        else
            thisLag=lag;
        end
        tmp=crossCorr(x(j,:)',y(j,:)',thisLag);
        c(j,:)=padNaN2(tmp(:,1)',2*lag+1);
    else
        c(j,:)=NaN*ones(1,2*lag+1);
    end
end
