function sp=findSigSeparationChange(idlist,dataProperties)

TAGS=[1 2];
t1=1;
while(isempty(idlist(t1).linklist))
    t1=t1+1;
end;
sp=[];
for t=t1+1:length(idlist)
    if ~isempty(idlist(t).linklist)
        if testSigSeparationChange(idlist,dataProperties,[t1 t],TAGS)
            sp=[sp; t1 t];
        end;
        t1=t;
    end
end;