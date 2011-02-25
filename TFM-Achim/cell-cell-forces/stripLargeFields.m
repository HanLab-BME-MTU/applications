function constrForceField=stripLargeFields(constrForceField,doSave)
if nargin<2 || isempty(doSave)
    doSave=0;
end
for frame=1:length(constrForceField)
    if ~isempty(constrForceField{frame}) && ~isempty(constrForceField{frame}.clusterAnalysis)
        constrForceField{frame}.clusterAnalysis.sol.u =[];    %These take too much space!
        constrForceField{frame}.clusterAnalysis.mesh.p=[];   %These take too much space!
        constrForceField{frame}.clusterAnalysis.mesh.e=[];   %These take too much space!
        constrForceField{frame}.clusterAnalysis.mesh.t=[];   %These take too much space!
        constrForceField{frame}.clusterAnalysis.par.globForce=[];
        constrForceField{frame}.clusterAnalysis.par.globYoung=[];
    end
end
if doSave
    save('cellCellForcesStripped.mat',constrForceField,'-v7.3')
end