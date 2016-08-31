function [] = whMetaDayAnalysis2016(allFeatures,healingRate,strLabels,metaData,mainDirname,plotBasicPCAFname)
warning('off','all');

addpath(genpath('/apps/MATLAB/R2015a/toolbox/stats/stats'));

[dayGeneDataSpeed,speedDiffs,healingRateSpeed] = ...
    dayAnalysisProperty(allFeatures.speedFeats,healingRate,strLabels,metaData,mainDirname,'Speed',plotBasicPCAFname);

[dayGeneDataDirectionality,directionalityDiffs,healingRateDirectionality] = ...
    dayAnalysisProperty(allFeatures.directionalityFeats,healingRate,strLabels,metaData,mainDirname,'Directionality',plotBasicPCAFname);

[dayGeneDataCoordination,coordinationDiffs,healingRateCoordination] = ...
    dayAnalysisProperty(allFeatures.coordinationFeats,healingRate,strLabels,metaData,mainDirname,'Coordination',plotBasicPCAFname);

% targetGenesStr = {'RHOA','TRIO','SOS1','ARHGEF11','TUBA','ARHGEF28'};
% whAssociatePropertyDirection(dayGeneDataSpeed,healingRateSpeed,speedDiffs,directionalityDiffs,coordinationDiffs,targetGenesStr,mainDirname);

close all;
end

%% dayAnalysisProperty
% TODO: update by day
%   1. Construct dayGeneData structure that includes the atomic unit of a daily experiment, print kymographs, etc. (till line 196)
%   DISCARDED   2. PCA of (gene KD - control) - line 200
%   DISCARDED   3. Assess variance
%   DISCARDED   4. dayAnalysis - PCA analysis: printing PCs, mean PCA per day, and for each condition
%   5. dayClassification - detect targtes (hits) via cluster seperation measures (in svm directory)
%   DISCARDED   6. dayClustering - clustering GEFs by phenotype (not active)
%   7. dayVisualizeKymographs - mean kymographs of control, KD and the subtraction KD - control (in dayGeneControlKymograph directory)
%   DISCARDED   8. whVisualizeTargets - PCA of (KD - control), highlighting targets. (in targets directory)
%   DISCARDED   9. whControlInterdayAssessment - assesseing variability of control
%   between days experiments (specifically targets), (in controlInterdayAssessment directory)
%   10. dayVisPCA - PCA of gene vs. control
%   11. dayVisDirectMeasures - direct meausres of gene vs. control

% TODO: Z-scores for step 3

function [dayGeneData,allDiffMeanGeneToMeanControl,healingRateOut] = dayAnalysisProperty(data,healingRate,strLabels,metaData,mainDirname,propertyStr,plotBasicPCAFname)
close all;

% flags.variance = 0;
% flags.dayAnalysisPCA = 0;
flags.dayScreen = 1;
flags.dayVisKymographs = 0;
flags.dayVisKymographsStd = 0;
% flags.dayVisTargets = 0;
% flags.controlInterdayAssessment = 0;

%% New analyses
flags.dayVisPCA = 0;
flags.dayVisDirectMeasures = 0;
flags.dayFollowupHits = 0;

[dayGeneData,allDiffMeanGeneToMeanControl,healingRateOut] = whCollectDayGeneData2016(data,healingRate,strLabels,metaData,mainDirname,propertyStr,plotBasicPCAFname);

save([mainDirname 'dayGeneData_' propertyStr '.mat'],'dayGeneData','allDiffMeanGeneToMeanControl','healingRateOut');


% %% 3. Assess variance
% if flags.variance
%     dayGeneData = whVariance(features,dayGeneData,indsPSup,mainDirname,propertyStr);
%     save([mainDirname propertyStr '_dayGeneData.mat'],'dayGeneData','allDiffMeanGeneToMeanControl','healingRateOut');
% end

% %% 4. dayAnalysis: PCA analysis: printing PCs, mean PCA per day, and for each condition KD - Control
% if flags.dayAnalysisPCA
%     dayAnalysis(allDiffVectors,allDiffToMeanControl,allDiffMeanGeneToMeanControl,dayGeneData,mainDirname,propertyStr);
% end

%% 5. whScreenbyDays - hit detection via cluster seperation measures
if flags.dayScreen
    whDayScreen2016(dayGeneData,mainDirname,propertyStr);
end

% %% 6. dayClustering: clustering GEFs by phenotype (not active)
% if flags.dayClustering
%     dayClustering(dayGeneData,allDiffVectors,allDiffMeanGeneToMeanControl,mainDirname,propertyStr);
% end

%% 7. dayVisualizeKymographs: mean kymographs of control, KD and the subtraction KD - control (in dayGeneControlKymograph directory)
if flags.dayVisKymographs
    dayVisualizeKymographs(dayGeneData,mainDirname,propertyStr,metaData);
end

if flags.dayVisKymographsStd && ~strcmp(propertyStr,'Coordination')
    dayVisualizeKymographsStd(dayGeneData,mainDirname,propertyStr,metaData);
end

% %% 8. whVisualizeTargets: PCA of (KD - control), highlighting targets. Output in targets directory.
% if flags.dayVisTargets
%     whVisualizeTargets(dayGeneData,allDiffMeanGeneToMeanControl,targetGenesStr,mainDirname,propertyStr);
% end
% 
% %% 9. whControlInterdayAssessment: 
% if flags.controlInterdayAssessment
%     whControlInterdayAssessment(dayGeneData,mainDirname,propertyStr,metaData,targetGenesStr,healingRateControl,healingRateGene);
% end

%% 10. day PCA
if flags.dayVisPCA
    whDayVisPCA2016(dayGeneData,mainDirname,propertyStr);
end
%% 11 direct meausres of gene vs. control
if flags.dayVisDirectMeasures
    whDayVisDirectMeasures(dayGeneData,mainDirname,propertyStr);
end

%% Followup hits
if flags.dayFollowupHits            
    %     validateGenes = {'Lcf'};
    
    % Screen hits, RhoGTPAses, SOS1 pathway, contractility inhibition
    validateGenes = {'beta-PIX','CDC42','RAC1',...
        'TRIO','SOS1',...
        'ARHGEF18','RHOA',...
        'ARHGEF11','ARHGEF28','ARHGEF3',...
        'PLEKHG2','PLEKHG6','BCR',...
        'GSK','PD','ERKi','PIK90',...
        'Y2763210uM','Y2763210uMp24',...
        'Blebbistatin10uM','Blebbistatin25uM','Blebbistatin50uM',...
        'Blebbistatin10uMp24','Blebbistatin25uMp24',...
        'Rhosin12dot5uM','Rhosin25uM','Rhosin50uM','Rhosin100uM','Rhosin200uM',...
        'FAK14dose1uM','FAK14dose5uM','FAK14dose10uM',...
        'SMIFH2dose5uM','SMIFH2dose10uM','SMIFH2dose15uM','SMIFH2dose25uM',...
        'Blebbistatin10uMp48','Y15uMp48','Y20uMp48','Y25uMp48',...
        'GSKp96','ERKip96','PDp96'...
        }; % 'Lcf' is out becuase the code asks for KD efficiency >= 50%...
    
    
    
    %     %% Followup 20160523: longer treatment times (SOS1 pathway, contractility inhibition), fishing for other phenotypes
    %     validateGenes = {'Rhosin12dot5uM','Rhosin25uM','Rhosin50uM','Rhosin100uM','Rhosin200uM',...
    %         'FAK14dose1uM','FAK14dose5uM','FAK14dose10uM',...
    %         'SMIFH2dose5uM','SMIFH2dose10uM','SMIFH2dose15uM','SMIFH2dose25uM',...
    %         'Blebbistatin10uMp48','Y15uMp48','Y20uMp48','Y25uMp48',...
    %         'GSKp96','ERKip96','PDp96',...
    %         };
        
    
    %% Screen hits, RhoGTPAses, SOS1 pathway, contractility inhibition
    %     validateGenes = {'beta-PIX','CDC42','RAC1',...
    %         'TRIO','SOS1',...
    %         'ARHGEF18','RHOA',...
    %         'ARHGEF11','ARHGEF28','ARHGEF3',...
    %         'PLEKHG2','PLEKHG6','BCR',...
    %         'GSK','PD','ERKi','PIK90',...
    %         'Y2763210uM','Y2763210uMp24',...
    %         'Blebbistatin10uM','Blebbistatin25uM','Blebbistatin50uM',...
    %         'Blebbistatin10uMp24','Blebbistatin25uMp24'}; % 'Lcf' is out becuase the code asks for KD efficiency >= 50%...
    
    %      validateGenes = {'Lcf'}; % when checking Lcf, make sure to
    %      change the threshold on the KD efficiency!
    %     if strcmp(propertyStr,'Speed')
    %         validateGenes = {'ARHGEF18','RHOA','TRIO','SOS1'}; % 'GSK','PD','ERKi','PIK90'
    %     else if strcmp(propertyStr,'Directionality') || strcmp(propertyStr,'Coordination')
     %             validateGenes = {'ARHGEF11','ARHGEF28','ARHGEF3','PLEKHG2','PLEKHG6','BCR','Y27632','Blebbistatin'}; % Blebbistatin - will include all dosages!!
     %         end
     %     end        
    whDayFollowupHits2016(dayGeneData,mainDirname,propertyStr,validateGenes);    
end


end