function display(obj)
% Display selected object properties
disp('Config: ===================');
fprintf('| configName: %s\n', obj.configName);
fprintf('| snapshotsEnabled: %f\n', obj.snapshotsEnabled);
fprintf('| displayEnabled: %f\n', obj.displayEnabled);
fprintf('| snapshotsPath: %s\n', obj.snapshotsPath);
fprintf('| path: %s\n', obj.path);
fprintf('| fileName: %s\n', obj.fileName);
fprintf('| errorX: %f\n', obj.errorX);
fprintf('| errorY: %f\n', obj.errorY);
fprintf('| errorZ: %f\n', obj.errorZ);
fprintf('| roiPosition: %f %f %f\n', obj.roiPosition);
fprintf('| roiSize: %.0f %.0f %.0f\n', obj.roiSize);
fprintf('| dataReductionEnabled: %u\n', obj.dataReductionEnabled);
fprintf('| nReductionRun: %u\n', obj.nReductionRun);
fprintf('| reductionEdgeRadius: %f\n', obj.reductionEdgeRadius);
fprintf('| densityFilteringEnabled: %u\n', obj.densityFilteringEnabled);
fprintf('| nNeighborsThreshold: %u\n', obj.nNeighborsThreshold);
fprintf('| neighborBallRadius: %f\n', obj.neighborBallRadius);
fprintf('| edgeWidthInitFree: %f\n', obj.edgeWidthInitFree);
fprintf('| filterLength: %f\n', obj.filterLength);
fprintf('| angularSampling: %f\n', obj.angularSampling);
fprintf('| maxDegreeBezier: %f\n', obj.maxDegreeBezier);
fprintf('| maxCurvature: %f\n', obj.maxCurvature);
fprintf('| betaVar: %f\n', obj.betaVar);
fprintf('| modeVar: %f\n', obj.modeVar);
fprintf('| fitMethod: %u\n', obj.fitMethod);
fprintf('| initialEdgeRadiusGeom: %f\n', obj.initialEdgeRadiusGeom);
fprintf('| initialEdgeRadius: %f\n', obj.initialEdgeRadius);
fprintf('| nIterGeomMatching: %f\n', obj.nIterGeomMatching);
fprintf('| modelLength: %f\n', obj.modelLength);
fprintf('| angleThreshold: %f\n', obj.angleThreshold);
fprintf('| nIterEM: %f\n', obj.nIterEM);
fprintf('| maxIterMerge: %f\n', obj.maxIterMerge);
fprintf('| nSigmaThreshold: %f\n', obj.nSigmaThreshold);
disp('===========================');
end
