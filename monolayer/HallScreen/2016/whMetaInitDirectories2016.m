
%% Analysis of the different treatments
function [dirs] = whMetaInitDirectories2016(mainDirname)
% ??
% allDataDirname = [mainDirname '../../allData/'];

dirs.hairpinEfficiencyDir = [mainDirname 'hairpinEfficiency/']; %

dirs.pcsDir = [mainDirname 'pcs/'];
% dirs.kymographDir = [mainDirname 'kymographs/'];
% dirs.pcaDir = [mainDirname 'pca/'];
% dirs.pcaCoeff = [dirs.pcaDir 'coeff/'];%
% dirs.pcaDayDir = [dirs.pcaDir 'day/'];
% dirs.pcaGeneDir = [dirs.pcaDir 'gene/'];
% dirs.pcaDayGeneDir = [dirs.pcaDir 'dayGene/'];
% dirs.associationDir = [mainDirname 'association/'];
dirs.healingRateDir = [mainDirname 'healingRate/'];
dirs.varianceDir = [mainDirname 'variance/'];
% dirs.dayDir = [mainDirname 'day/'];
dirs.dayScreenDir = [mainDirname 'screen/'];
% dirs.dayClusterDir = [mainDirname 'clusters/'];
% dirs.dayGeneDir = [mainDirname 'day/gene/'];
% dirs.micrscopeRepeatErrorStatsDir = [mainDirname 'micrscopeError/'];
dirs.targetsDir = [mainDirname 'targets/'];
dirs.dayGeneControlKymographDir = [mainDirname 'dayGeneControlKymograph/'];
dirs.dayGeneControlKymographStdDir = [mainDirname 'dayGeneControlKymographStd/'];
% dirs.dayGenePlotsDir = [mainDirname 'dayGenePlots/'];
dirs.dayWellReplicatesKymographDir = [mainDirname 'dayWellReplicatesKymographs/'];
% dirs.combinedPropertiesDir = [mainDirname 'combinedProperties/'];
% dirs.combinedPropertiesAnglesPcsDir = [mainDirname 'combinedProperties/anglesPC/'];
% dirs.controlInterdayAssessmentDir = [mainDirname 'controlInterdayAssessment/'];
% dirs.plithotaxisDir = [mainDirname 'plithotaxis/'];
% dirs.plithotaxisOutDir = [mainDirname 'plithotaxisOut/'];
% dirs.RhoGTPasesDir = [mainDirname 'RhoGTPases/'];

dirs.dayGeneControlPCA = [mainDirname 'dayGeneControlPCA'];
dirs.dayGeneControlDirect = [mainDirname 'dayGeneControlDirectMeasure'];% magnitude, space and time derivatives
dirs.dayGeneControlFollowup = [mainDirname 'dayGeneControlFollowup'];

% input (no output/writing)
% dirs.correctMotionDir = [allDataDirname 'correctMotion/']; %

if ~exist(dirs.pcsDir,'dir')
    mkdir(dirs.pcsDir);
end

if ~exist(dirs.hairpinEfficiencyDir,'dir')
    mkdir(dirs.hairpinEfficiencyDir);
end

% if ~exist(dirs.kymographDir,'dir')
%     mkdir(dirs.kymographDir);
% end
% 
% if ~exist(dirs.pcaDir,'dir')
%     mkdir(dirs.pcaDir);
% end
% 
% if ~exist(dirs.pcaCoeff,'dir')
%     mkdir(dirs.pcaCoeff);
% end
% 
% if ~exist(dirs.pcaDayDir,'dir')
%     mkdir(dirs.pcaDayDir);
% end
% 
% if ~exist(dirs.pcaDayGeneDir,'dir')
%     mkdir(dirs.pcaDayGeneDir);
% end
% 
% if ~exist(dirs.pcaGeneDir,'dir')
%     mkdir(dirs.pcaGeneDir);
% end
% 
% if ~exist(dirs.associationDir,'dir')
%     mkdir(dirs.associationDir);
% end
% 
if ~exist(dirs.healingRateDir,'dir')
    mkdir(dirs.healingRateDir);
end
% 
% if ~exist(dirs.varianceDir,'dir')
%     mkdir(dirs.varianceDir);
% end
% 
% if ~exist(dirs.dayDir,'dir')
%     mkdir(dirs.dayDir);
% end
% 
if ~exist(dirs.dayScreenDir,'dir')
    mkdir(dirs.dayScreenDir);
end
% 
% if ~exist(dirs.dayClusterDir,'dir')
%     mkdir(dirs.dayClusterDir);
% end
% 
% if ~exist(dirs.dayGeneDir,'dir')
%     mkdir(dirs.dayGeneDir);
% end
% 
% if ~exist(dirs.micrscopeRepeatErrorStatsDir,'dir')
%     mkdir(dirs.micrscopeRepeatErrorStatsDir);
% end
% 
% if ~exist(dirs.targetsDir,'dir')
%     mkdir(dirs.targetsDir);
% end
% 
if ~exist(dirs.dayGeneControlKymographDir,'dir')
    mkdir(dirs.dayGeneControlKymographDir);
end

if ~exist(dirs.dayGeneControlKymographStdDir,'dir')
    mkdir(dirs.dayGeneControlKymographStdDir);
end
% 
% if ~exist(dirs.dayGenePlotsDir,'dir')
%     mkdir(dirs.dayGenePlotsDir);
% end
% 
if ~exist(dirs.dayWellReplicatesKymographDir,'dir')
    mkdir(dirs.dayWellReplicatesKymographDir);
end
% 
% if ~exist(dirs.combinedPropertiesDir,'dir')
%     mkdir(dirs.combinedPropertiesDir);
% end
% 
% if ~exist(dirs.combinedPropertiesAnglesPcsDir,'dir')
%     mkdir(dirs.combinedPropertiesAnglesPcsDir);
% end
% 
% if ~exist(dirs.controlInterdayAssessmentDir,'dir')
%     mkdir(dirs.controlInterdayAssessmentDir);
% end
% 
% if ~exist(dirs.plithotaxisOutDir,'dir')
%     mkdir(dirs.plithotaxisOutDir);
% end
% 
% if ~exist(dirs.RhoGTPasesDir,'dir')
%     mkdir(dirs.RhoGTPasesDir);
% end

propertiesStr = {'Speed','Directionality','Coordination'};

if ~exist(dirs.dayGeneControlPCA,'dir')
    mkdir(dirs.dayGeneControlPCA);
    for iprop = 1 : 3
        curProp = propertiesStr{iprop};
        mkdir([dirs.dayGeneControlPCA filesep curProp]);
    end
end

if ~exist(dirs.dayGeneControlDirect,'dir')
    mkdir(dirs.dayGeneControlDirect);
    for iprop = 1 : 3
        curProp = propertiesStr{iprop};
        mkdir([dirs.dayGeneControlDirect filesep curProp]);
    end
end

if ~exist(dirs.dayGeneControlFollowup,'dir')
    mkdir(dirs.dayGeneControlFollowup);
    for iprop = 1 : 3
        curProp = propertiesStr{iprop};
        mkdir([dirs.dayGeneControlDirect filesep curProp]);
    end
end

end