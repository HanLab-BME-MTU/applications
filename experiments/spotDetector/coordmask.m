function li=coordmask(k,mask)

d=floor(mask/2);
ct=1;
idx=1:size(k,1)';
for i=1:size(k,1)-1
    if(idx(i)~=0)
        li=k(i+1:end,:);
        dif= abs(li-ones(size(li,1),1)*k(i,:)) - ones(size(li,1),1)*d;
        idx(i+1:end)=idx(i+1:end)'.* all(dif>0,2);
    end
end;
li=k(nonzeros(idx),:);
