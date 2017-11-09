function alignments = alignProtrusionEvents(velocity,activity,dThreshold)
%alignProtrusionEvents detects retraction onsets, protrusion onsets, and
%maximal protrusion velocity from protrusion velocity map, and align
%velocity and activity with respect to retraction onsets, protrusion onsets, and
%maximal protrusion velocity

%  Inputs
%  velocity: protrusion velocity map
%  activity: activity map
%  dThreshold: distance fluctuation less than dThreshld is ignored.
%  (Default: 10)



warning off;

%check number of input arguments
if nargin < 1
    disp('alignProtrusionEvents: Please enter arguments');
    return
end

%dThresold
if nargin < 3 || isempty(dThreshold)
    dThreshold=10;
end


vThreshold=3;

[nWin nTime]=size(velocity);
pEvents=getProtrusionEvents(velocity,dThreshold);
alignments.pEvents=pEvents;
alignments.protrusionOnset.velocity=NaN*ones(1,2*nTime+1);
alignments.protrusionOnset.activity=NaN*ones(1,2*nTime+1);
alignments.betweenPOnsetPMaxVelocity.velocity=NaN*ones(1,2*nTime+1);
alignments.betweenPOnsetPMaxVelocity.activity=NaN*ones(1,2*nTime+1);
alignments.retractionOnset.velocity=NaN*ones(1,2*nTime+1);
alignments.retractionOnset.activity=NaN*ones(1,2*nTime+1);
alignments.protrusionMaxVelocity.velocity=NaN*ones(1,2*nTime+1);
alignments.protrusionMaxVelocity.activity=NaN*ones(1,2*nTime+1);
alignments.protrusionMaxVelocity2.velocity=NaN*ones(1,2*nTime+1);
alignments.protrusionMaxVelocity2.activity=NaN*ones(1,2*nTime+1);
alignments.protrusionMaxVelocityFromProtrusionOnsetData.velocity=NaN*ones(1,2*nTime+1);
alignments.protrusionMaxVelocityFromProtrusionOnsetData.activity=NaN*ones(1,2*nTime+1);
alignments.normlaizedPOnsetPMaxVelocity.velocity=NaN*ones(1,101);
alignments.normlaizedPOnsetPMaxVelocity.activity=NaN*ones(1,101);
alignments.activityFromRO=NaN*ones(1,200);
alignments.velocityFromRO=NaN*ones(1,200);
alignments.maxV=NaN;

splineParam=0.5;

% for normalized time from protrusion onset to protrusion max velocity
%for j=1:nWin
%    if sum(isfinite(velocity(j,:)))>2
%        prot_spline(j)= csaps(linspace(1,nTime,nTime),velocity(j,:),splineParam);
        %sProt(j,:)=ppval(act1_spline,linspace(1,nTime,nTime));
    
%        act_spline(j)= csaps(linspace(1,nTime,nTime),activity(j,:),splineParam);
        %sActivity(j,:)=ppval(act1_spline,linspace(1,nTime,nTime));
%    end
%end
% for normalized time from protrusion onset to protrusion max velocity


nProt=0;
nRet=0;
nProtV=0;
nProtV2=0;
nProtRV=0;
nPose=0;
nStep=0;
nMidProtProtV=0;

nRetProt=0;
for j=1:nWin
    nCycle=pEvents.events{j}.nCycle;
    for k=1:nCycle
        % align protrusion onsets
        rIndex=pEvents.events{j}.cycle{k}.protrusionOnset; %reference index
        if isfinite(rIndex)
            if isfinite(pEvents.events{j}.cycle{k}.retractionOnset)
                iIndex=pEvents.events{j}.cycle{k}.retractionOnset;
            else
                iIndex=1;
            end
            if isfinite(pEvents.events{j}.cycle{k}.protrusionEnd)
                fIndex=pEvents.events{j}.cycle{k}.protrusionEnd;
            else
                fIndex=nTime;
            end

            %if isfinite(pEvents.events{j}.cycle{k}.protrusionEnd)
                nProt=nProt+1;
                alignments.protrusionOnset.velocity(nProt,:)=NaN*ones(1,2*nTime+1);
                alignments.protrusionOnset.activity(nProt,:)=NaN*ones(1,2*nTime+1);
                alignments.protrusionOnset.velocity(nProt,(nTime+1):(nTime+fIndex-rIndex+1))=velocity(j,rIndex:fIndex);
                alignments.protrusionOnset.activity(nProt,(nTime+1):(nTime+fIndex-rIndex+1))=activity(j,rIndex:fIndex);
                alignments.protrusionOnset.velocity(nProt,(nTime-rIndex+iIndex+1):nTime)=velocity(j,iIndex:(rIndex-1));
                alignments.protrusionOnset.activity(nProt,(nTime-rIndex+iIndex+1):nTime)=activity(j,iIndex:(rIndex-1));
                
                alignments.protrusionOnset.window(nProt)=j;
                alignments.protrusionOnset.time(nProt)=rIndex;
                alignments.protrusionOnset.protrusionDistance(nProt)=pEvents.events{j}.cycle{k}.protrusionDistance;
                
                iIndex=pEvents.events{j}.cycle{k}.protrusionOnset;
                rIndex=pEvents.events{j}.cycle{k}.pMaxVelocityIndex;
                alignments.protrusionMaxVelocityFromProtrusionOnsetData.velocity(nProt,:)=NaN*ones(1,2*nTime+1);
                alignments.protrusionMaxVelocityFromProtrusionOnsetData.activity(nProt,:)=NaN*ones(1,2*nTime+1);
                if (rIndex>iIndex) && (rIndex<fIndex)
                    alignments.protrusionMaxVelocityFromProtrusionOnsetData.velocity(nProt,(nTime+1):(nTime+fIndex-rIndex+1))=velocity(j,rIndex:fIndex);
                    alignments.protrusionMaxVelocityFromProtrusionOnsetData.activity(nProt,(nTime+1):(nTime+fIndex-rIndex+1))=activity(j,rIndex:fIndex);
                    alignments.protrusionMaxVelocityFromProtrusionOnsetData.velocity(nProt,(nTime-rIndex+iIndex+1):nTime)=velocity(j,iIndex:(rIndex-1));
                    alignments.protrusionMaxVelocityFromProtrusionOnsetData.activity(nProt,(nTime-rIndex+iIndex+1):nTime)=activity(j,iIndex:(rIndex-1));
                end
            %end
            %   rIndex= round(0.25*pEvents.events{j}.cycle{k}.protrusionOnset+0.75*pEvents.events{j}.cycle{k}.pMaxVelocityIndex);
            %   if isfinite(rIndex)
            %        nMidProtProtV=nMidProtProtV+1;
            %        alignments.betweenPOnsetPMaxVelocity.velocity(nMidProtProtV,:)=NaN*ones(1,2*nTime+1);
            %        alignments.betweenPOnsetPMaxVelocity.activity(nMidProtProtV,:)=NaN*ones(1,2*nTime+1);
            %        alignments.betweenPOnsetPMaxVelocity.velocity(nMidProtProtV,(nTime+1):(nTime+fIndex-rIndex+1))=velocity(j,rIndex:fIndex);
            %        alignments.betweenPOnsetPMaxVelocity.activity(nMidProtProtV,(nTime+1):(nTime+fIndex-rIndex+1))=activity(j,rIndex:fIndex);
            %        alignments.betweenPOnsetPMaxVelocity.velocity(nMidProtProtV,(nTime-rIndex+iIndex+1):nTime)=velocity(j,iIndex:(rIndex-1));
            %        alignments.betweenPOnsetPMaxVelocity.activity(nMidProtProtV,(nTime-rIndex+iIndex+1):nTime)=activity(j,iIndex:(rIndex-1));
                
            %        alignments.betweenPOnsetPMaxVelocity.window(nMidProtProtV)=j;
            %        alignments.betweenPOnsetPMaxVelocity.time(nMidProtProtV)=rIndex;
            %   end
            
            
            % for normalized time from protrusion onset to protrusion max velocity
            %iIndex=pEvents.events{j}.cycle{k}.protrusionOnset;
            %fIndex=pEvents.events{j}.cycle{k}.pMaxVelocityIndex;
                        
            %if (fIndex-iIndex)>10
            %    nMidProtProtV=nMidProtProtV+1;
            %    alignments.normlaizedPOnsetPMaxVelocity.velocity(nMidProtProtV,:)=NaN*ones(1,101);
            %    alignments.normlaizedPOnsetPMaxVelocity.activity(nMidProtProtV,:)=NaN*ones(1,101);             
            %    normTime=iIndex+(0:0.01:1)*(fIndex-iIndex);
            %    alignments.normlaizedPOnsetPMaxVelocity.velocity(nMidProtProtV,:)=ppval(prot_spline(j),normTime);  
            %    alignments.normlaizedPOnsetPMaxVelocity.activity(nMidProtProtV,:)=ppval(act_spline(j),normTime);   
            %end
            % for normalized time from protrusion onset to protrusion max velocity
            
            
        end
        % align retraction onsets
        rIndex=pEvents.events{j}.cycle{k}.retractionOnset;
        if isfinite(rIndex)
            if k==1;
                iIndex=1;
            elseif isnan(pEvents.events{j}.cycle{k-1}.protrusionOnset)
                iIndex=1;
            else
                iIndex=pEvents.events{j}.cycle{k-1}.protrusionOnset;
            end
            if isfinite(pEvents.events{j}.cycle{k}.protrusionOnset)
                fIndex=pEvents.events{j}.cycle{k}.protrusionOnset;
            else
                fIndex=nTime;
            end

            nRet=nRet+1;
            alignments.retractionOnset.velocity(nRet,:)=NaN*ones(1,2*nTime+1);
            alignments.retractionOnset.activity(nRet,:)=NaN*ones(1,2*nTime+1);
            alignments.retractionOnset.velocity(nRet,(nTime+1):(nTime+fIndex-rIndex+1))=velocity(j,rIndex:fIndex);
            alignments.retractionOnset.activity(nRet,(nTime+1):(nTime+fIndex-rIndex+1))=activity(j,rIndex:fIndex);
            alignments.retractionOnset.velocity(nRet,(nTime-rIndex+iIndex+1):nTime)=velocity(j,iIndex:(rIndex-1));
            alignments.retractionOnset.activity(nRet,(nTime-rIndex+iIndex+1):nTime)=activity(j,iIndex:(rIndex-1));
                
            alignments.retractionOnset.window(nRet)=j;
            alignments.retractionOnset.time(nRet)=rIndex;
        end
        
        
        
        % align max protrusion velocity in each protrusion event
        rIndex=pEvents.events{j}.cycle{k}.pMaxVelocityIndex;
        if isfinite(rIndex)
            if isfinite(pEvents.events{j}.cycle{k}.protrusionOnset)
                iIndex=pEvents.events{j}.cycle{k}.protrusionOnset;
            else
                iIndex=1;
            end
                        
            if k==nCycle
                fIndex=nTime;
            elseif isfinite(pEvents.events{j}.cycle{k+1}.protrusionOnset)
               fIndex=pEvents.events{j}.cycle{k+1}.protrusionOnset;
            else
                fIndex=nTime;
            end
            
           
            nProtV=nProtV+1;
            alignments.protrusionMaxVelocity.velocity(nProtV,:)=NaN*ones(1,2*nTime+1);
            alignments.protrusionMaxVelocity.activity(nProtV,:)=NaN*ones(1,2*nTime+1);
            alignments.protrusionMaxVelocity.velocity(nProtV,(nTime+1):(nTime+fIndex-rIndex+1))=velocity(j,rIndex:fIndex);
            alignments.protrusionMaxVelocity.activity(nProtV,(nTime+1):(nTime+fIndex-rIndex+1))=activity(j,rIndex:fIndex);
            alignments.protrusionMaxVelocity.velocity(nProtV,(nTime-rIndex+iIndex+1):nTime)=velocity(j,iIndex:(rIndex-1));
            alignments.protrusionMaxVelocity.activity(nProtV,(nTime-rIndex+iIndex+1):nTime)=activity(j,iIndex:(rIndex-1));
                
            alignments.protrusionMaxVelocity.window(nProtV)=j;
            alignments.protrusionMaxVelocity.time(nProtV)=rIndex;
        end       

        % align max protrusion velocity in each protrusion event
        rIndex=pEvents.events{j}.cycle{k}.pMaxVelocityIndex;
        if isfinite(rIndex)
            if isfinite(pEvents.events{j}.cycle{k}.protrusionOnset)
                iIndex=pEvents.events{j}.cycle{k}.protrusionOnset;
            else
                iIndex=1;
            end
                        
            
            if isfinite(pEvents.events{j}.cycle{k}.protrusionEnd)
                fIndex=pEvents.events{j}.cycle{k}.protrusionEnd;
            else
                fIndex=nTime;
            end
            
            nProtV2=nProtV2+1;
            alignments.protrusionMaxVelocity2.velocity(nProtV2,:)=NaN*ones(1,2*nTime+1);
            alignments.protrusionMaxVelocity2.activity(nProtV2,:)=NaN*ones(1,2*nTime+1);
            alignments.protrusionMaxVelocity2.velocity(nProtV2,(nTime+1):(nTime+fIndex-rIndex+1))=velocity(j,rIndex:fIndex);
            alignments.protrusionMaxVelocity2.activity(nProtV2,(nTime+1):(nTime+fIndex-rIndex+1))=activity(j,rIndex:fIndex);
            alignments.protrusionMaxVelocity2.velocity(nProtV2,(nTime-rIndex+iIndex+1):nTime)=velocity(j,iIndex:(rIndex-1));
            alignments.protrusionMaxVelocity2.activity(nProtV2,(nTime-rIndex+iIndex+1):nTime)=activity(j,iIndex:(rIndex-1));
                
            alignments.protrusionMaxVelocity2.window(nProtV2)=j;
            alignments.protrusionMaxVelocity2.time(nProtV2)=rIndex;
        end   
        
        
        
        pOIndex=pEvents.events{j}.cycle{k}.protrusionOnset;
        rOIndex=pEvents.events{j}.cycle{k}.retractionOnset;
        if (~isnan(pOIndex))&&(~isnan(rOIndex))
            nRetProt=nRetProt+1;
            
            %velocitySpline= csaps(linspace(1,pOIndex-rOIndex+1,pOIndex-rOIndex+1),velocity(j,rOIndex:pOIndex),0.5);
            %activitySpline= csaps(linspace(1,pOIndex-rOIndex+1,pOIndex-rOIndex+1),activity(j,rOIndex:pOIndex),0.5);
            %activityFilteredSpline= csaps(linspace(1,nTime,nTime),activity(j,:),0.5);
            %activityFiltered=ppval(activityFilteredSpline,linspace(1,nTime,nTime));
            %velocityFilteredSpline= csaps(linspace(1,nTime,nTime),velocity(j,:),0.5);
            %velocityFiltered=ppval(velocityFilteredSpline,linspace(1,nTime,nTime));
            %timeArray=(pOIndex-rOIndex)/100*linspace(0,100,101);
            %alignments.velocityTimeNormalized(nRetProt,:)=ppval(velocitySpline,timeArray);
            %alignments.activityTimeNormalized(nRetProt,:)=ppval(activitySpline,timeArray);
            
            alignments.maxV(nRetProt,1)=max(velocity(j,pOIndex:end));%max(velocityFiltered(pOIndex:end)); 
            alignments.retractionTime(nRetProt,1)=pOIndex-rOIndex;
            
            alignments.activityFromRO(nRetProt,:)=NaN*ones(1,200);
            alignments.velocityFromRO(nRetProt,:)=NaN*ones(1,200);
            alignments.activityFromRO(nRetProt,1:(pOIndex-rOIndex+1))=activity(j,rOIndex:pOIndex);%activityFiltered(rOIndex:pOIndex); 
            alignments.velocityFromRO(nRetProt,1:(pOIndex-rOIndex+1))=velocity(j,rOIndex:pOIndex);%velocityFiltered(rOIndex:pOIndex); 
            %alignments.activityToPO(nRetProt,(200-(pOIndex-rOIndex)):200)=activity(j,rOIndex:pOIndex);
            %alignments.velocityToPO(nRetProt,(200-(pOIndex-rOIndex):200))=velocity(j,rOIndex:pOIndex);
        end
        
        % align max retraction velocity in each retraction event
        %rIndex=pEvents.events{j}.cycle{k}.retractionMaxVelocityIndex;
        %if isfinite(rIndex)
        %    if isfinite(pEvents.events{j}.cycle{k}.retractionOnset)
        %        iIndex=pEvents.events{j}.cycle{k}.retractionOnset;
        %    else
        %        iIndex=1;
        %    end
        %    if isfinite(pEvents.events{j}.cycle{k}.protrusionOnset)
        %        fIndex=pEvents.events{j}.cycle{k}.protrusionOnset;
        %    else
        %        fIndex=nTime;
        %    end
        %    nProtRV=nProtRV+1;
        %    alignments.retractionMaxVelocity.velocity(nProtRV,:)=NaN*ones(1,2*nTime+1);
        %    alignments.retractionMaxVelocity.activity(nProtRV,:)=NaN*ones(1,2*nTime+1);
        %    alignments.retractionMaxVelocity.velocity(nProtRV,(nTime+1):(nTime+fIndex-rIndex+1))=velocity(j,rIndex:fIndex);
        %    alignments.retractionMaxVelocity.activity(nProtRV,(nTime+1):(nTime+fIndex-rIndex+1))=activity(j,rIndex:fIndex);
        %    alignments.retractionMaxVelocity.velocity(nProtRV,(nTime-rIndex+iIndex+1):nTime)=velocity(j,iIndex:(rIndex-1));
        %    alignments.retractionMaxVelocity.activity(nProtRV,(nTime-rIndex+iIndex+1):nTime)=activity(j,iIndex:(rIndex-1));
        %end  
        
        
        % align protrusion pose events
        %[m nPoseCycle]=size(pEvents.events{j}.cycle{k}.protrusionPose);
        %for l=1:nPoseCycle
        %    rIndex=pEvents.events{j}.cycle{k}.protrusionPose(l);
        %    if l==1
        %       iIndex=pEvents.events{j}.cycle{k}.protrusionOnset;             
        %    else
        %       iIndex=pEvents.events{j}.cycle{k}.protrusionPose(l-1); 
        %    end
        %    if l==nPoseCycle
        %        fIndex=pEvents.events{j}.cycle{k}.protrusionEnd;
        %        if isnan(fIndex)
        %            fIndex=nTime;
        %        end
        %    else
        %        fIndex=pEvents.events{j}.cycle{k}.protrusionPose(l+1); 
        %    end
        %    nPose=nPose+1;
        %    alignments.protrusionPose.velocity(nPose,:)=NaN*ones(1,2*nTime+1);
        %    alignments.protrusionPose.activity(nPose,:)=NaN*ones(1,2*nTime+1);
        %    alignments.protrusionPose.velocity(nPose,(nTime+1):(nTime+fIndex-rIndex+1))=velocity(j,rIndex:fIndex);
        %    alignments.protrusionPose.activity(nPose,(nTime+1):(nTime+fIndex-rIndex+1))=activity(j,rIndex:fIndex);
        %    alignments.protrusionPose.velocity(nPose,(nTime-rIndex+iIndex+1):nTime)=velocity(j,iIndex:(rIndex-1));
        %    alignments.protrusionPose.activity(nPose,(nTime-rIndex+iIndex+1):nTime)=activity(j,iIndex:(rIndex-1));            
            
        %end
        
        % align max protrusion velocity events in each protrusion step
        %nStepCycle=pEvents.events{j}.cycle{k}.nStep;
        %if ~isnan(nStepCycle) 
        %    for l=1:nStepCycle
        %    
        %        iIndex=pEvents.events{j}.cycle{k}.protrusionStep{l}.onset;
        %        fIndex=pEvents.events{j}.cycle{k}.protrusionStep{l}.end;
        %        if isnan(fIndex)
        %            fIndex=nTime;
        %        end
        %        rIndex=pEvents.events{j}.cycle{k}.protrusionStep{l}.maxVelocity;
        %        nStep=nStep+1;
        %        alignments.protrusionStepMaxVelocity.velocity(nStep,:)=NaN*ones(1,2*nTime+1);
        %        alignments.protrusionStepMaxVelocity.activity(nStep,:)=NaN*ones(1,2*nTime+1);
        %        alignments.protrusionStepMaxVelocity.velocity(nStep,(nTime+1):(nTime+fIndex-rIndex+1))=velocity(j,rIndex:fIndex);
        %        alignments.protrusionStepMaxVelocity.activity(nStep,(nTime+1):(nTime+fIndex-rIndex+1))=activity(j,rIndex:fIndex);
        %        alignments.protrusionStepMaxVelocity.velocity(nStep,(nTime-rIndex+iIndex+1):nTime)=velocity(j,iIndex:(rIndex-1));
        %        alignments.protrusionStepMaxVelocity.activity(nStep,(nTime-rIndex+iIndex+1):nTime)=activity(j,iIndex:(rIndex-1));    
        %    end
        %end
    end

end

