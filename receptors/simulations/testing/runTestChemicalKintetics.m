%   Script to run testChemicalKinetics.
%

rateAssoc = 10;
rateDissoc = 5;
timeStep = 0.01;
totalTime = 100;

simNum = 1;
type1 = zeros(simNum,1);

typeN1 = zeros(simNum,1);
n1To1 = zeros(simNum,1);

for simIter=1:simNum
    [clusterSize,rateAssocComp,rateDissocComp,type1(simIter,1),typeN1(simIter,1),transitionType] = testChemicalKinetics_Mod(...
    rateAssoc,rateDissoc,timeStep,totalTime);
    n1To1(simIter,1) = typeN1(simIter,1)/type1(simIter,1);
end

csumNumAssocEvents = cumsum(transitionType == 1);
csumNumDissocEvents = cumsum(transitionType == -1);
csumProbAssoc = csumNumAssocEvents./(csumNumAssocEvents+csumNumDissocEvents);
csumProbDissoc = csumNumDissocEvents./(csumNumAssocEvents+csumNumDissocEvents);
figure()
semilogx(1:length(csumProbDissoc),csumProbDissoc,'k:',1:length(csumProbAssoc),csumProbAssoc,'r:');