function events = detectProtrusionEvents(v,dThreshold)
%UNTITLED Summary of this function goes here
%   

%check number of input arguments
if nargin < 1
    disp('detectProtrusionEvents: Please enter arguments');
    return
end

%dThresold
if nargin < 2 || isempty(dThreshold)
    dThreshold=10;
end

tThreshold=10;
maxvThreshold=0.2;
minvThreshold=inf;
splineParam=0.01;
splineParam2=0.1;

%   calculate protrusion distance from protrusion velocity
d=getProtrusionDistance(v,1);

%   perform spline filter for protrusion distance
[nWindow nTime]=size(v);
sd_spline= csaps(linspace(1,nTime,nTime),d,splineParam);
sd=ppval(sd_spline,linspace(1,nTime,nTime));

%   calculate protrusion velocity from the filtered distance using analytic
%   expression
sd_spline2= csaps(linspace(1,nTime,nTime),d,splineParam2);
sd2=ppval(sd_spline2,linspace(1,nTime,nTime));

sv_spline=sd_spline2;
sv_spline.order=3;
sv_spline.coefs(:,1)=3*sd_spline2.coefs(:,1);
sv_spline.coefs(:,2)=2*sd_spline2.coefs(:,2);
sv_spline.coefs(:,3)=1*sd_spline2.coefs(:,3);
sv=ppval(sv_spline,linspace(1,nTime,nTime));


%   find local maxima of distance
[dmaxPks dmaxLocs]=findpeaks(sd,'MINPEAKDISTANCE',tThreshold);
[m nDmaxPks]=size(dmaxPks);

%   find minia between maxima
iDminPks=1;dminPks=[];dminLocs=[];
iMinCand=0;minCandPks=[];minCandLocs=[];
decreasing=0;
oldLeftPeak=NaN;
for j=1:(nDmaxPks+1)

    if j==1
        range(1)=1;
    else
        range(1)=dmaxLocs(j-1);
    end
    if j==(nDmaxPks+1)
        range(2)=nTime;
    else 
        range(2)=dmaxLocs(j);
    end
    if (range(2)-range(1))>3
        [pks locs]=findpeaks(-sd(range(1):range(2)),'SORTSTR','descend');
        [m nPks]=size(pks);
    else
        nPks=0;
    end
    if nPks>0
        leftPeak=sd(range(1));
        rightPeak=sd(range(2));
        if leftPeak<oldLeftPeak
            leftPeak=oldLeftPeak;
        else
            oldLeftPeak=leftPeak;
        end
        leftHeight=leftPeak+pks(1);
        rightHeight=sd(range(2))+pks(1);
        
        if (leftHeight>dThreshold)&&(rightHeight>dThreshold) % significant local minimum
            iMinCand=iMinCand+1;
            minCandPks(iMinCand)=-pks(1);minCandLocs(iMinCand)=range(1)+locs(1)-1;
            [minVal minIndex]=min(minCandPks);
            dminPks(iDminPks)=minVal;
            dminLocs(iDminPks)=minCandLocs(minIndex);
            iDminPks=iDminPks+1;
            
            iMinCand=0;minCandPks=[];minCandLocs=[];
            decreasing=0;
            oldLeftPeak=rightPeak;
        elseif  (leftHeight>dThreshold)&&(rightHeight<dThreshold)
            iMinCand=iMinCand+1;minCandPks(iMinCand)=-pks(1);minCandLocs(iMinCand)=range(1)+locs(1)-1;
            decreasing=1;
        elseif (leftHeight<dThreshold)&&(rightHeight>dThreshold)    
            if decreasing==1
                iMinCand=iMinCand+1;
                minCandPks(iMinCand)=-pks(1);minCandLocs(iMinCand)=range(1)+locs(1)-1;
                [minVal minIndex]=min(minCandPks);
                dminPks(iDminPks)=minVal;
                dminLocs(iDminPks)=minCandLocs(minIndex);
                iDminPks=iDminPks+1;
                oldLeftPeak=rightPeak;
                decreasing=0;
            end
            iMinCand=0;minCandPks=[];minCandLocs=[];
        else 
            if decreasing==1
                iMinCand=iMinCand+1;
                minCandPks(iMinCand)=-pks(1);minCandLocs(iMinCand)=range(1)+locs(1)-1;
            end
        end
    end
    
end

%   find maxia between minima
iDmaxPks=1;
dmaxPks=[];
dmaxLocs=[];
[m nDminPks]=size(dminPks);
for j=1:(nDminPks+1)

    if j==1
        range(1)=1;
    else
        range(1)=dminLocs(j-1);
    end
    if j==(nDminPks+1)
        range(2)=nTime;
    else 
        range(2)=dminLocs(j);
    end
    if (range(2)-range(1))>tThreshold
        [pks locs]=findpeaks(sd(range(1):range(2)),'SORTSTR','descend');
        [m nPks]=size(pks);
        if (nPks>0)&&(nDminPks>0)
             if (j==1)||(j==(nDminPks+1))
                if j==1
                    [minPks minLocs]=min(sd((range(1):range(1)+locs(1)-1)));
                else
                    [minPks minLocs]=min(sd((range(1)+locs(1)-1):range(2)));
                end
                [m nMinPks]=size(minPks);
                if nMinPks>1
                    sideMinPks=minPks(1);
                else
                    sideMinPks=minPks;
                end
                if (pks(1)-sideMinPks)>dThreshold
                    dmaxPks(iDmaxPks)=pks(1);
                    dmaxLocs(iDmaxPks)=range(1)+locs(1)-1;
                    iDmaxPks=iDmaxPks+1;                   
                end
             else
                dmaxPks(iDmaxPks)=pks(1);
                dmaxLocs(iDmaxPks)=range(1)+locs(1)-1;
                iDmaxPks=iDmaxPks+1;
             end
        end
    end
    
end


% store minima and maxima of distance
dmin=[dminLocs' dminPks'];
dmax=[dmaxLocs' dmaxPks'];

% remove insignificant max and min
[ndmax m]=size(dmax);[ndmin m]=size(dmin);
%for j=1:ndmax
%    for k=1:ndmin
%        if isfinite(dmin(k,1))&isfinite(dmax(j,1))
%            if abs(dmax(j,1)-dmin(k,1))<0.5*tThreshold
%                dmax(j,1)=NaN;dmax(j,2)=NaN;
%                dmin(k,1)=NaN;dmin(k,2)=NaN;
%            end
%        end
%    end
%end

% remove NaNs
if ndmin>0
    dmin=dmin(isfinite(dmin(:,1)),:);
end
if ndmax>0
    dmax=dmax(isfinite(dmax(:,1)),:);
end

[ndmax m]=size(dmax);[ndmin m]=size(dmin);
n = max([ndmax ndmin]);

% find cycles
if isempty(dmax)&&isempty(dmin)
    cycle{1}.retractionOnset=NaN;
    cycle{1}.protrusionOnset=NaN;
    cycle{1}.protrusionEnd=NaN;   
elseif isempty(dmax)&&(~isempty(dmin))
    cycle{1}.retractionOnset=NaN;
    cycle{1}.protrusionOnset=dmin(1,1);
    cycle{1}.protrusionEnd=NaN;
elseif (~isempty(dmax))&&isempty(dmin)
    cycle{1}.retractionOnset=NaN;
    cycle{1}.protrusionOnset=NaN;
    cycle{1}.protrusionEnd=dmax(1,1);     
    cycle{2}.retractionOnset=dmax(1,1);
    cycle{2}.protrusionOnset=NaN;
    cycle{2}.protrusionEnd=NaN;   
elseif dmax(1,1)<dmin(1,1)
    cycle{1}.retractionOnset=NaN;
    cycle{1}.protrusionOnset=NaN;
    cycle{1}.protrusionEnd=dmax(1,1);
    for j=1:n
        cycle{j+1}.retractionOnset=dmax(j,1);
        if j==n
            if j > ndmin 
                cycle{j+1}.protrusionOnset=NaN;
            else
                cycle{j+1}.protrusionOnset=dmin(j,1);
            end
            cycle{j+1}.protrusionEnd=NaN;
        else
            cycle{j+1}.protrusionOnset=dmin(j,1);
            cycle{j+1}.protrusionEnd=dmax(j+1,1);
        end
    end
else
    cycle{1}.retractionOnset=NaN;
    cycle{1}.protrusionOnset=dmin(1,1);
    cycle{1}.protrusionEnd=dmax(1,1);
    for j=2:n
        cycle{j}.retractionOnset=dmax(j-1,1);
        cycle{j}.protrusionOnset=dmin(j,1);
        if j==n
            if j > ndmax 
                cycle{j}.protrusionEnd=NaN; 
            else
                cycle{j}.protrusionEnd=dmax(j,1);
            end
        else
             cycle{j}.protrusionEnd=dmax(j,1);
        end
    end
end

[m nCycle]=size(cycle);

% final quantifications
for j=1:nCycle
    cycle{j}.protrusionMaxVelocityIndex=NaN;
    cycle{j}.retractionMaxVelocityIndex=NaN;
    cycle{j}.maxPVelocity=NaN;
    if isnan(cycle{j}.retractionOnset)||isnan(cycle{j}.protrusionOnset)
        cycle{j}.retractionTime=NaN;
        cycle{j}.retractionDistance=NaN;
        cycle{j}.rDist=NaN;
        cycle{j}.avgRetractionVelocity=NaN;
        cycle{j}.maxRetractionVelocity=NaN;
        cycle{j}.maxRetractionVelocityTime=NaN;
        cycle{j}.v=NaN;
        cycle{j}.sv=NaN;
        cycle{j}.retractVProfile=[];
        cycle{j}.nRetractStep=NaN;
        if isfinite(cycle{j}.retractionOnset)&&isnan(cycle{j}.protrusionOnset)
            cycle{j}.rDist=d(cycle{j}.retractionOnset)-d(nTime);
        end
        
    else
        fIndex=cycle{j}.protrusionOnset;
        iIndex=cycle{j}.retractionOnset;
        cycle{j}.retractionTime=fIndex-iIndex;
        cycle{j}.retractionDistance=d(iIndex)-d(fIndex);
        cycle{j}.rDist=d(iIndex)-d(fIndex);
        cycle{j}.avgRetractionVelocity=nanmean(v(iIndex:fIndex));
        cycle{j}.maxRetractionVelocity=min(v(iIndex:fIndex));
        cycle{j}.retractV=v(iIndex:fIndex);
        cycle{j}.retractSV=sv(iIndex:fIndex);
        vProfile=getVelocityProfile(-sv(iIndex:fIndex),0.5*tThreshold,maxvThreshold,minvThreshold);
        vProfile.maxVelocityTime=vProfile.maxVelocityTime+iIndex-1;
        vProfile.minVelocityTime=vProfile.minVelocityTime+iIndex-1;
        cycle{j}.retractVProfile=vProfile;
        cycle{j}.nRetractStep=vProfile.nMaxVelocity;

       
    end
    if isnan(cycle{j}.protrusionOnset)||isnan(cycle{j}.protrusionEnd)
        cycle{j}.protrusionTime=NaN;
        cycle{j}.protrusionDistance=NaN;
        cycle{j}.pDist=NaN;
        cycle{j}.avgProtrusionVelocity=NaN;
        cycle{j}.maxProtrusionVelocity=NaN;
        cycle{j}.maxProtrusionVelocityReal=NaN;
        cycle{j}.pMaxVelocityIndex=NaN;
        cycle{j}.avgProtrusionVelocityBeforeMaxV=NaN;
        cycle{j}.avgProtrusionVelocityAfterMaxV=NaN;
        
        cycle{j}.v=NaN;
        cycle{j}.sv=NaN;
        cycle{j}.nStep=NaN;
        cycle{j}.protVProfile=[];
        if isfinite(cycle{j}.protrusionOnset)&&isnan(cycle{j}.protrusionEnd)
            cycle{j}.pDist=d(nTime)-d(cycle{j}.protrusionOnset);
            [maxV maxVIndex]=max(sv(cycle{j}.protrusionOnset:nTime));
            cycle{j}.pMaxVelocityIndex=maxVIndex+cycle{j}.protrusionOnset-1;
            cycle{j}.maxPVelocity=v(maxVIndex+cycle{j}.protrusionOnset-1);
            if (nTime-cycle{j}.pMaxVelocityIndex)<5
                cycle{j}.pMaxVelocityIndex=NaN;
                cycle{j}.maxPVelocity=NaN;
            end
        end
 
    else
        fIndex=cycle{j}.protrusionEnd;
        iIndex=cycle{j}.protrusionOnset;
        cycle{j}.protrusionTime=fIndex-iIndex;
        cycle{j}.protrusionDistance=d(fIndex)-d(iIndex);
        cycle{j}.pDist=d(fIndex)-d(iIndex);
        cycle{j}.avgProtrusionVelocity=nanmean(v(iIndex:fIndex));
        [maxV maxVIndex]=max(sv(iIndex:fIndex));
        [maxVReal maxVRealIndex]=max(v(iIndex:fIndex));
        cycle{j}.protrusionMaxVelocityIndex=maxVIndex+iIndex-1;
        cycle{j}.pMaxVelocityIndex=maxVIndex+iIndex-1;
        cycle{j}.maxProtrusionVelocity=v(maxVIndex+iIndex-1);
        cycle{j}.maxProtrusionVelocityReal=maxVReal;
        cycle{j}.maxPVelocity=v(maxVIndex+iIndex-1);
  
        cycle{j}.avgProtrusionVelocityBeforeMaxV=nanmean(v(iIndex:(maxVIndex+iIndex-2)));
        cycle{j}.avgProtrusionVelocityAfterMaxV=nanmean(v((maxVIndex+iIndex-1):fIndex));
        
        
        cycle{j}.protV=v(iIndex:fIndex);
        cycle{j}.protSV=sv(iIndex:fIndex);
      
        vProfile=getVelocityProfile(sv(iIndex:fIndex),0.5*tThreshold,maxvThreshold,minvThreshold);
        vProfile.maxVelocityTime=vProfile.maxVelocityTime+iIndex-1;
        vProfile.minVelocityTime=vProfile.minVelocityTime+iIndex-1;
        cycle{j}.protVProfile=vProfile;
        cycle{j}.nStep=vProfile.nMaxVelocity;


    end
    
    % for max velocity without having protrusion boundary
    if isnan(cycle{j}.protrusionOnset)&&isfinite(cycle{j}.protrusionEnd)
        fIndex=cycle{j}.protrusionEnd;
        iIndex=1;
        [maxV maxVIndex]=max(sv(iIndex:fIndex));
        cycle{j}.pMaxVelocityIndex=maxVIndex+iIndex-1;
        cycle{j}.maxPVelocity=v(maxVIndex+iIndex-1);
        cycle{j}.pDist=sd(fIndex)-sd(iIndex);
    elseif isnan(cycle{j}.protrusionOnset)&&isnan(cycle{j}.protrusionEnd)&&((d(nTime)-d(1))>30)
        fIndex=nTime;
        iIndex=1;
        [maxV maxVIndex]=max(sv(iIndex:fIndex));
        cycle{j}.pMaxVelocityIndex=maxVIndex+iIndex-1;    
        cycle{j}.maxPVelocity=v(maxVIndex+iIndex-1);
        if (nTime-cycle{j}.pMaxVelocityIndex)<5
            cycle{j}.pMaxVelocityIndex=NaN;
            cycle{j}.maxPVelocity=NaN;
        else
            cycle{j}.pDist=sd(fIndex)-sd(iIndex);
        end
    end
        
end


events.tThreshold=tThreshold;
events.splineParam=splineParam;
events.d=d;
events.v=v;
events.sd=sd;
events.advancement=sd(nTime)-sd(1);
events.sv=sv;
events.dmin=dmin;
events.dmax=dmax;
events.nCycle=nCycle;
events.cycle=cycle;

end

