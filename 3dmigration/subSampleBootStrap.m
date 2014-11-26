function samps = subSampleBootStrap(nSamp,func,data,nSub)


%NOTE: Vectorize this at some point. In too much of a hurry now.
n=numel(data);

samps = nan(nSamp,1);
for j = 1:nSamp    
    samps(j) = func(data(randsample(n,nSub)));   
end



