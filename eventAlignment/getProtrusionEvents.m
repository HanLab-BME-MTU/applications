function protrusion = getProtrusionEvents( protVelocity, dThreshold)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%check number of input arguments
if nargin < 1
    disp('getProtrusionEvents: Please enter arguments');
    return
end

%dThresold
if nargin < 2 || isempty(dThreshold)
    dThreshold=10;
end


[nWindow nTime]=size(protVelocity);

iCycle=0;
iStep=0;
poseVelocity=[];
poseRetractVelocity=[];

for j =1:nWindow

    events{j}=detectProtrusionEvents(protVelocity(j,:),dThreshold);
    winData(j,1)=events{j}.nCycle;
    winData(j,2)=events{j}.advancement;
    [m nCycle]=size(events{j}.cycle);
    totalStep=0;
    maxProtDistance=NaN;
    maxRetDistance=NaN;
    for k=1:nCycle
        iCycle=iCycle+1;
        cycleData(iCycle,1)=j;
        cycleData(iCycle,2)=k;
        cycleData(iCycle,3)=events{j}.cycle{k}.protrusionDistance;
        cycleData(iCycle,4)=events{j}.cycle{k}.protrusionTime;
        cycleData(iCycle,5)=events{j}.cycle{k}.avgProtrusionVelocity;
        cycleData(iCycle,6)=events{j}.cycle{k}.maxProtrusionVelocity;
        cycleData(iCycle,7)=events{j}.cycle{k}.retractionDistance;
        cycleData(iCycle,8)=events{j}.cycle{k}.retractionTime; 
        cycleData(iCycle,9)=events{j}.cycle{k}.avgRetractionVelocity; 
        cycleData(iCycle,10)=events{j}.cycle{k}.maxRetractionVelocity; 
        cycleData(iCycle,11)=events{j}.cycle{k}.nStep; %totalStep=totalStep+events{j}.cycle{k}.nStep;
        cycleData(iCycle,12)=events{j}.cycle{k}.nRetractStep;
        cycleData(iCycle,13)=events{j}.cycle{k}.protrusionMaxVelocityIndex-events{j}.cycle{k}.protrusionOnset;
        cycleData(iCycle,14)=events{j}.cycle{k}.avgProtrusionVelocityBeforeMaxV;
        cycleData(iCycle,15)=events{j}.cycle{k}.avgProtrusionVelocityAfterMaxV;
        cycleData(iCycle,16)=events{j}.cycle{k}.maxProtrusionVelocityReal;
        maxProtDistance=max([maxProtDistance events{j}.cycle{k}.protrusionDistance]);
        maxRetDistance=max([maxRetDistance events{j}.cycle{k}.retractionDistance]);
        %      cycleData(iCycle,12)=events{j}.cycle{k}.protrusionOnsetAcceleration; 
  %      cycleData(iCycle,13)=events{j}.cycle{k}.maxRetractionVelocityTime; 
        if ~isempty(events{j}.cycle{k}.protVProfile)
            poseVelocity=[poseVelocity events{j}.cycle{k}.protVProfile.minVelocity];
        end
        if ~isempty(events{j}.cycle{k}.retractVProfile)
            poseRetractVelocity=[poseRetractVelocity -events{j}.cycle{k}.retractVProfile.minVelocity];
        end
  
  %      if ~isempty(events{j}.cycle{k}.protrusionStep)
  %          [m nStep]=size(events{j}.cycle{k}.protrusionStep);

  %          for l=1:nStep
  %              iStep=iStep+1;
  %              stepData(iStep,1)=j;
  %              stepData(iStep,2)=k;
  %              stepData(iStep,3)=l;
  %              stepData(iStep,4)=events{j}.cycle{k}.protrusionStep{l}.distance;
  %              stepData(iStep,5)=events{j}.cycle{k}.protrusionStep{l}.time;
  %              stepData(iStep,6)=events{j}.cycle{k}.protrusionStep{l}.maxVelocityValue;
  %              stepData(iStep,7)=events{j}.cycle{k}.protrusionStep{l}.maxVelocityTime;
  %          end
  %      end
    end  
  % winData(j,2)=totalStep;
  winData(j,3)=maxProtDistance;
  winData(j,4)=maxRetDistance;
  sv(j,:)=events{j}.sv;
end

protrusion.events=events;
protrusion.sv=sv;
protrusion.winData=winData;
protrusion.cycleData=cycleData;
%protrusion.stepData=stepData;
protrusion.poseVelocity=poseVelocity;
protrusion.poseRetractVelocity=poseRetractVelocity;

protrusion.nCycle=winData(:,1);
protrusion.advancement=winData(:,2);
protrusion.maxProtDistance=winData(:,3);
protrusion.maxRetDistance=winData(:,4);
protrusion.protrusionDistance=cycleData(:,3);
protrusion.protrusionTime=cycleData(:,4);
protrusion.avgProtrusionVelocity=cycleData(:,5);
protrusion.maxProtrusionVelocity=cycleData(:,6);
protrusion.retractionDistance=cycleData(:,7);
protrusion.retractionTime=cycleData(:,8);
protrusion.avgRetractionVelocity=cycleData(:,9);
protrusion.maxRetractionVelocity=cycleData(:,10);
protrusion.nStep=cycleData(:,11);
protrusion.nRetractStep=cycleData(:,12);
protrusion.maxVelocityTime=cycleData(:,13);
protrusion.avgProtrusionVelocityBeforeMaxV=cycleData(:,14);
protrusion.avgProtrusionVelocityAfterMaxV=cycleData(:,15);
protrusion.maxProtrusionVelocityReal=cycleData(:,16);
%protrusion.protrusionOnsetAcceleration=cycleData(:,12);
%protrusion.maxRetractionVelocityTime=cycleData(:,13);
%protrusion.stepDistance=stepData(:,4);
%protrusion.stepTime=stepData(:,5);
%protrusion.stepMaxVelocity=stepData(:,6);

protrusion.nCycleMean=nanmean(winData(:,1));
protrusion.advancementMean=nanmean(winData(:,2));
protrusion.protrusionDistanceMean=nanmean(cycleData(:,3));
protrusion.protrusionTimeMean=nanmean(cycleData(:,4));
protrusion.avgProtrusionVelocityMean=nanmean(cycleData(:,5));
protrusion.maxProtrusionVelocityMean=nanmean(cycleData(:,6));
protrusion.retractionDistanceMean=nanmean(cycleData(:,7));
protrusion.retractionTimeMean=nanmean(cycleData(:,8));
protrusion.avgRetractionVelocityMean=nanmean(cycleData(:,9));
protrusion.maxRetractionVelocityMean=nanmean(cycleData(:,10));
protrusion.nStepMean=nanmean(cycleData(:,11));
protrusion.nRetractStepMean=nanmean(cycleData(:,12));
protrusion.avgProtrusionVelocityBeforeMaxVMean=nanmean(cycleData(:,14));
protrusion.avgProtrusionVelocityAfterMaxVMean=nanmean(cycleData(:,15));
protrusion.maxProtrusionVelocityRealMean=nanmean(cycleData(:,16));

protrusion.poseVelocityMean=nanmean(poseVelocity);
protrusion.poseRetractVelocityMean=nanmean(poseRetractVelocity);
%protrusion.protrusionOnsetAccelerationMean=nanmean(cycleData(:,12));
%protrusion.maxRetractionVelocityTimeMean=nanmean(cycleData(:,13));
%protrusion.stepDistanceMean=nanmean(stepData(:,4));
%protrusion.stepTimeMean=nanmean(stepData(:,5));
%protrusion.stepMaxVelocityMean=nanmean(stepData(:,6));
%protrusion.stepMaxVelocityTimeMean=nanmean(stepData(:,7));




