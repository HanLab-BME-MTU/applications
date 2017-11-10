function [result]=compareBasisFunctions(ForceMesh1,ForceMesh2,basisID)
test1 = (ForceMesh1.basis(basisID).class==ForceMesh2.basis(basisID).class);
test2 = compPts(ForceMesh1.basis(basisID).node,ForceMesh2.basis(basisID).node);
test3 = (ForceMesh1.basis(basisID).nodeID==ForceMesh2.basis(basisID).nodeID);

result = test1 && test2 && test3;