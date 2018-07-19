function activityLabFrame = sampleMovieWindowsLabFrame(velocity,activity,movieData,channelIndex, row)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%dThreshold=0;

[nWin nTime]=size(velocity);
pEvents=getProtrusionEvents(velocity);
activityLabFrame=activity;
velocityLabFrame=activity;


for j=1:nWin
    j
    nCycle=pEvents.events{j}.nCycle;
    
    for k=1:nCycle

        iIndex=pEvents.events{j}.cycle{k}.protrusionOnset; %reference index
        if isfinite(iIndex)
            if isfinite(pEvents.events{j}.cycle{k}.protrusionEnd)
                fIndex=pEvents.events{j}.cycle{k}.protrusionEnd;
            else
                fIndex=nTime;
            end
            tmp=sampleMovieLabFrame(movieData,iIndex,j,row,channelIndex,[iIndex fIndex]);
            if nnz(isnan(tmp))>0
                activityLabFrame(j,iIndex:fIndex)=NaN;
            else
                activityLabFrame(j,iIndex:fIndex)=tmp;
            end
            velocityLabFrame(j,iIndex:fIndex)=velocity(j,iIndex:fIndex);
        end

        if isfinite(pEvents.events{j}.cycle{k}.retractionOnset)
            iIndex=pEvents.events{j}.cycle{k}.retractionOnset;
            if isfinite(pEvents.events{j}.cycle{k}.protrusionOnset)
                fIndex=pEvents.events{j}.cycle{k}.protrusionOnset;
            else
                fIndex=nTime;
            end
            
        elseif isfinite(pEvents.events{j}.cycle{k}.protrusionOnset)
            iIndex=1;
            fIndex=pEvents.events{j}.cycle{k}.protrusionOnset;
       
        else
            iIndex=1;
            fIndex=nTime;
        end
        activityLabFrame(j,iIndex:fIndex)=NaN;
        velocityLabFrame(j,iIndex:fIndex)=NaN;
    end

end

