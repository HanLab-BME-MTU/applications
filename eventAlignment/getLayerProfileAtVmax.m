function layerProfile = getLayerProfileAtVmax(v,sample,shiftFromVmax)

[nWin nLayer nTime]=size(sample);

pEvents=getProtrusionEvents(v);
iLayerProfile=1;
for j=1:nWin
    nCycle=pEvents.events{j}.nCycle;
    for k=1:nCycle
        vmaxTime=pEvents.events{j}.cycle{k}.pMaxVelocityIndex;
        if isfinite(vmaxTime)
            layerProfile(iLayerProfile,:)=sample(j,:,vmaxTime+shiftFromVmax);
            iLayerProfile=iLayerProfile+1;
        end
    end
end


