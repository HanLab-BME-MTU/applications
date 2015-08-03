function [clusterSize,rateAssocComp,rateDissocComp,type1,typeN1,transitionType] = testChemicalKinetics_Mod(...
    rateAssoc,rateDissoc,timeStep,totalTime)

%get probabilities from rates
probAssoc = rateAssoc * timeStep;
probDissoc = rateDissoc * timeStep;

%get number of steps
numSteps = ceil(totalTime/timeStep);

%initialize output
clusterSize = ones(numSteps,1);

%determine when there are association or dissociation events
associationEvent = rand(numSteps,1) < probAssoc;
dissociationEvent = rand(numSteps,1) < probDissoc;

%simulate cluster size over time
for iStep = 2 : numSteps
    clusterSize(iStep) = clusterSize(iStep-1) + associationEvent(iStep) - dissociationEvent(iStep);
end

%recover association and dissociation rates ...

%find association and dissociation (i.e. transition) time points, and thus cluster lifetimes
clusterSizeDiff = diff(clusterSize);
transitionPoints = find(clusterSizeDiff~=0);
clusterLifetime = diff(transitionPoints);
lifetimeMean2Any = mean(clusterLifetime)*timeStep;
transitionType = clusterSizeDiff(transitionPoints(2:end));
%
probAssocComp = length(find(transitionType==1))/length(transitionType);
probDissocComp = length(find(transitionType==-1))/length(transitionType);
%
%thus calculate back association and dissociation rates
rateAssocComp = probAssocComp / lifetimeMean2Any;
rateDissocComp = probDissocComp / lifetimeMean2Any;

type1 = length(find(transitionType == 1));
typeN1 = length(find(transitionType == -1));




