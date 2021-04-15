% Sue's code for saving all publishable figures
% This is to compare the PIVs and PTVs with ground-truth filed 
% written by Sangyoon Han, March 2021
cd('/Volumes/GoogleDrive/My Drive/Documents/Faculty/Projects/PTV-Retracking/Manuscript/Data')
%% Load ground truth ux and uy
uxStruct = load('ux_sim.mat');
ux = uxStruct.ux;
uyStruct = load('uy_sim.mat');
uy = uyStruct.uy;
x_mat_uStruct = load('x_mat_u.mat');
x_mat_u=x_mat_uStruct.x_mat_u;
y_mat_uStruct = load('y_mat_u.mat');
y_mat_u = y_mat_uStruct.y_mat_u;
xR = 32:8:480; yR = xR;
% imshow(ones(512,512)); hold on
quiver(x_mat_u(xR,yR),y_mat_u(xR,yR), ux(xR,yR), uy(xR,yR), 0,'Color','b')
set(gca, 'YDir','reverse')
set(gca, 'XLim',[32 480], 'YLim',[32 480])
set(gca, 'PlotBoxAspectRatio',[1 1 1])
%% map for org field
generateHeatmapFromGridData(x_mat_u,y_mat_u,ux,uy,'./orgField',32,0,60,0)
%% Mean Squared Error for PIVSuite
displFieldObj = load('./analysisWithoutSDC/TFMPackage/displacementField_PIVSuite/displField.mat');
displFieldSuite = displFieldObj.displField;
quiver(displFieldSuite(1).pos(:,1),displFieldSuite(1).pos(:,2),...
    displFieldSuite(1).vec(:,1),displFieldSuite(1).vec(:,2),0,'Color','r')
set(gca, 'YDir','reverse')
set(gca, 'XLim',[32 480], 'YLim',[32 480])
set(gca, 'PlotBoxAspectRatio',[1 1 1])
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])
%% Compare for PIVSuite
displFieldOrg.x_mat_u = x_mat_u;
displFieldOrg.y_mat_u = y_mat_u;
displFieldOrg.ux = ux;
displFieldOrg.uy = uy;
plotRes = true;

[deviationLevelsSuite,indWrongVectorsSuite,indMissingSuite, MSESuite] =  ...
    compareFields(displFieldSuite,displFieldOrg,plotRes);
intactVectorsSuite = ~(indWrongVectorsSuite | indMissingSuite);
figure;
quiver(displFieldSuite(1).pos(intactVectorsSuite,1),...
    displFieldSuite(1).pos(intactVectorsSuite,2),...
    displFieldSuite(1).vec(intactVectorsSuite,1),...
    displFieldSuite(1).vec(intactVectorsSuite,2),0,'Color','b')
hold on
quiver(displFieldSuite(1).pos(indWrongVectorsSuite,1),...
    displFieldSuite(1).pos(indWrongVectorsSuite,2),...
    displFieldSuite(1).vec(indWrongVectorsSuite,1),...
    displFieldSuite(1).vec(indWrongVectorsSuite,2),0,'Color','g')
plot(displFieldSuite(1).pos(indMissingSuite,1),...
    displFieldSuite(1).pos(indMissingSuite,2),'ro')
set(gca, 'YDir','reverse')
set(gca, 'XLim',[32 480], 'YLim',[32 480])
set(gca, 'PlotBoxAspectRatio',[1 1 1])
%% map for PIVSuite
dataPath='./analysisWithoutSDC/TFMPackage/displacementField_PIVSuite';
generateHeatmapFromField(displFieldSuite,dataPath,0,60); % step 
%% Tseng's PIV
%% Displacement filed for ImageJ PIV Results
% displOriginalImageJstruct = load('/storage/network/TFM_Development/TFM2D/PIVimprovement/simulation/TFMPackage/displacementFieldImageJ/dispFieldImageJ.mat');
dispPositionImageJstruct = load('dispPositionImageJ.mat');
dispVectorImageJstruct = load('dispVectorImageJ.mat');

% quiver plot - ImageJ result
figure(3), quiver(dispPositionImageJstruct.positionImageJ(:,1),dispPositionImageJstruct.positionImageJ(:,2),...
    dispVectorImageJstruct.vectorImageJ(:,1),dispVectorImageJstruct.vectorImageJ(:,2),0,'Color','b')
% hold on
% generateHeatmapFromGridData(dispPositionImageJstruct.positionImageJ(:,1),dispPositionImageJstruct.positionImageJ(:,2),dispVectorImageJstruct.vectorImageJ(:,1),dispVectorImageJstruct.vectorImageJ(:,2),[dataPath '/Figures'],0,0,30,false,500,500);
displFieldTseng.pos = dispPositionImageJstruct.positionImageJ;
displFieldTseng.vec = dispVectorImageJstruct.vectorImageJ;
[deviationLevels_tseng,indWrongVectors_tseng,indMissing_tseng,MSE_tseng] = ...
    compareFields(displFieldTseng,displFieldOrg,plotRes);
%% map for Tseng's PIV
dataPath='./analysisWithoutSDC/TFMPackage/displacementField_TsengPIV';
generateHeatmapFromField(displFieldTseng,dataPath,0,60); % step 
%% mpiv
%load displacement field
iuStruct = load('iu.mat'); % displacement vector calculated by mpiv
ivStruct = load('iv.mat');  % displacement vector calculated by mpiv
iu = iuStruct.iu;
iv = ivStruct.iv;

%load position
xiStruct = load('xi.mat'); % 1x1 struct
yiStruct = load('yi.mat'); % 1x1 struct
xi = xiStruct.xi;
yi = yiStruct.yi;

%make them 30x30
[xgrid,ygrid] = meshgrid(xi,yi);

%quiver plot
figure, quiver(xgrid,ygrid,iu,iv,0,'Color','g')
hold on
%% mpiv - load filtered displacement field
iu_ft_Struct = load('iu_ft.mat');
iv_ft_Struct = load('iv_ft.mat');
iu_ft = iu_ft_Struct.iu_ft;
iv_ft = iv_ft_Struct.iv_ft;

quiver(xgrid,ygrid,iu_ft,iv_ft,0,'Color','b')
hold on
%% mpiv - load interpolated displacement field
load('missingPos_mpiv.mat')
iu_ip_Struct = load('iu_ip.mat');
iv_ip_Struct = load('iv_ip.mat');
iu_ip = iu_ip_Struct.iu_ip;
iv_ip = iv_ip_Struct.iv_ip;

plot(missingPos_mpiv(:,1),missingPos_mpiv(:,2),'ro')
%% MSE for mpiv
displFieldMPIV.pos = [ygrid(:) xgrid(:)];
displFieldMPIV.vec = [iv_ft(:) iu_ft(:)];

[deviationLevels_mpiv,indWrongVectors_mpiv,indMissing_mpiv,MSE_mpiv] = ...
    compareFields(displFieldMPIV,displFieldOrg,plotRes);

%% colormap - mpiv
dataPath='./analysisWithoutSDC/TFMPackage/displacementField_mpiv';
generateHeatmapFromField(displFieldMPIV,dataPath,0,60); % step 
%% Displacement field for our original PTV
displOriginalPTVstruct = load('./analysisWithoutSDC/TFMPackage/displacementField_PTV/displField.mat');
displOriginalPTV = displOriginalPTVstruct.displField;

%% MSE for original PTV
[deviationLevels_ptv,indWrongVectors_ptv,indMissing_ptv,MSE_ptv] = ...
    compareFields(displOriginalPTV,displFieldOrg,plotRes);
%% quiver plot for PTV
intactVectorsPTV = ~(indWrongVectors_ptv | indMissing_ptv);
figure;
quiver(displOriginalPTV(1).pos(intactVectorsPTV,1),displOriginalPTV(1).pos(intactVectorsPTV,2),...
    displOriginalPTV(1).vec(intactVectorsPTV,1),displOriginalPTV(1).vec(intactVectorsPTV,2),0,'Color','b')
hold on
quiver(displOriginalPTV(1).pos(indWrongVectors_ptv,1),displOriginalPTV(1).pos(indWrongVectors_ptv,2),...
    displOriginalPTV(1).vec(indWrongVectors_ptv,1),displOriginalPTV(1).vec(indWrongVectors_ptv,2),0,'Color','g')
plot(displOriginalPTV(1).pos(indMissing_ptv,1),displOriginalPTV(1).pos(indMissing_ptv,2),'ro')
set(gca, 'YDir','reverse')
set(gca, 'XLim',[32 480], 'YLim',[32 480])
set(gca, 'PlotBoxAspectRatio',[1 1 1])
%% map for ptv
dataPath='./analysisWithoutSDC/TFMPackage/displacementField_PTV';
generateHeatmapFromField(displOriginalPTV,dataPath,0,60); % step 

%% Filtered vectors
displOriginalPTVFilteredstruct = load('./analysisWithoutSDC/TFMPackage/correctedDisplacementField_justFiltering/displField.mat');
displOriginalPTVFiltered = displOriginalPTVFilteredstruct.displField;
%% MSE for filtered PTV
[deviationLevels_ptvFilt,indWrongVectors_ptvFilt,indMissing_ptvFilt,MSE_ptvFilt] = ...
    compareFields(displOriginalPTVFiltered,displFieldOrg,plotRes);
%% quiver plot for filtered PTV
intactVectorsPTVFilt = ~(indWrongVectors_ptvFilt | indMissing_ptvFilt);
figure;
quiver(displOriginalPTVFiltered(1).pos(intactVectorsPTVFilt,1),displOriginalPTVFiltered(1).pos(intactVectorsPTVFilt,2),...
    displOriginalPTVFiltered(1).vec(intactVectorsPTVFilt,1),displOriginalPTVFiltered(1).vec(intactVectorsPTVFilt,2),0,'Color','b')
hold on
quiver(displOriginalPTVFiltered(1).pos(indWrongVectors_ptvFilt,1),displOriginalPTVFiltered(1).pos(indWrongVectors_ptvFilt,2),...
    displOriginalPTVFiltered(1).vec(indWrongVectors_ptvFilt,1),displOriginalPTVFiltered(1).vec(indWrongVectors_ptvFilt,2),0,'Color','g')
plot(displOriginalPTVFiltered(1).pos(indMissing_ptvFilt,1),displOriginalPTVFiltered(1).pos(indMissing_ptvFilt,2),'ro')
set(gca, 'YDir','reverse')
set(gca, 'XLim',[32 480], 'YLim',[32 480])
set(gca, 'PlotBoxAspectRatio',[1 1 1])
%% map for ptv filtered
dataPath='./analysisWithoutSDC/TFMPackage/correctedDisplacementField_justFiltering';
generateHeatmapFromField(displOriginalPTVFiltered,dataPath,0,60); % step 

%% MSE plot
MSEGroup = {MSESuite, MSE_tseng, MSE_mpiv, MSE_ptv,MSE_ptvFilt};
namePIVs = {'PIVSuite','TsengPIV','mpiv','cPTV','cPTV-filtered'};
figure, barPlotCellArray(MSEGroup,namePIVs,1)
ylabel('Mean Squared Deviation (1)')
title('MSD in FOV')

%% MSE for small and large adhesions
% We need a mask for large and small adhesions
open('./orgField/heatMap/figs/hMapFig.fig');
h = drawellipse;
maskLargeAdh = createMask(h);
h2 = drawpolygon;
maskSmallAdh = createMask(h2);
% figure, imshow(maskSmallAdh)
%% Large adh first
dF_Suite_LA = filterDisplacementField( displFieldSuite, maskLargeAdh);
[~,~,~, MSELA{1}] =  compareFields(dF_Suite_LA,displFieldOrg);
dF_Tseng_LA = filterDisplacementField( displFieldTseng, maskLargeAdh);
[~,~,~, MSELA{2}] =  compareFields(dF_Tseng_LA,displFieldOrg);
dF_mpiv_LA = filterDisplacementField( displFieldMPIV, maskLargeAdh);
[~,~,~, MSELA{3}] =  compareFields(dF_mpiv_LA,displFieldOrg);
dF_ptv_LA = filterDisplacementField( displOriginalPTV, maskLargeAdh);
[~,~,~, MSELA{4}] =  compareFields(dF_ptv_LA,displFieldOrg);
dF_ptvf_LA = filterDisplacementField( displOriginalPTVFiltered, maskLargeAdh);
[~,~,~, MSELA{5}] =  compareFields(dF_ptvf_LA,displFieldOrg);
figure, barPlotCellArray(MSELA,namePIVs,1)
ylabel('Mean Squared Deviation (1)')
title('MSD at large forces')
%% Small adhs
dF_Suite_SA = filterDisplacementField( displFieldSuite, maskSmallAdh);
[~,~,~, MSESA{1}] =  compareFields(dF_Suite_SA,displFieldOrg);
dF_Tseng_SA = filterDisplacementField( displFieldTseng, maskSmallAdh);
[~,~,~, MSESA{2}] =  compareFields(dF_Tseng_SA,displFieldOrg);
dF_mpiv_SA = filterDisplacementField( displFieldMPIV, maskSmallAdh);
[~,~,~, MSESA{3}] =  compareFields(dF_mpiv_SA,displFieldOrg);
dF_ptv_SA = filterDisplacementField( displOriginalPTV, maskSmallAdh);
[~,~,~, MSESA{4}] =  compareFields(dF_ptv_SA,displFieldOrg);
dF_ptvf_SA = filterDisplacementField( displOriginalPTVFiltered, maskSmallAdh);
[~,~,~, MSESA{5}] =  compareFields(dF_ptvf_SA,displFieldOrg);
figure, barPlotCellArray(MSESA,namePIVs,1)
ylabel('Mean Squared Deviation (1)')
title('MSD at small forces')
%% Fig 4 Enlargement factor
% We want to compare MSE for them
% loading it again
%% Retrack without enlargement factor
displPTVretractedStruct = load('./analysisWithoutSDC/TFMPackage/correctedDisplacementField_retrackedUpperMedian/displField.mat');
displPTVretracted = displPTVretractedStruct.displField;
%% MSE for filtered PTV
[deviationLevels_ptvR,indWrongVectors_ptvR,indMissing_ptvR,MSEptvR] = ...
    compareFields(displPTVretracted,displFieldOrg,plotRes);
figure
intactVectors_ptvR = ~(indWrongVectors_ptvR | indMissing_ptvR);
quiver(displPTVretracted(1).pos(intactVectors_ptvR,1),displPTVretracted(1).pos(intactVectors_ptvR,2),...
    displPTVretracted(1).vec(intactVectors_ptvR,1),displPTVretracted(1).vec(intactVectors_ptvR,2),0,'Color','b')
hold on
quiver(displPTVretracted(1).pos(indWrongVectors_ptvR,1),displPTVretracted(1).pos(indWrongVectors_ptvR,2),...
    displPTVretracted(1).vec(indWrongVectors_ptvR,1),displPTVretracted(1).vec(indWrongVectors_ptvR,2),0,'Color','g')
plot(displPTVretracted(1).pos(indMissing_ptvR,1),displPTVretracted(1).pos(indMissing_ptvR,2),'ro')
set(gca, 'YDir','reverse')
set(gca, 'XLim',[32 480], 'YLim',[32 480])
set(gca, 'PlotBoxAspectRatio',[1 1 1])
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])
%% Enla
displPTVretractedEnlargedStruct = load('./analysisWithoutSDC/TFMPackage/correctedDisplacementField_enlarge7_cor_nei_med/displField.mat');
displPTVretractedEnlarged = displPTVretractedEnlargedStruct.displField;
%% MSE for PTVR enlargementFactor
[deviationLevels_ptvEnla,indWrongVectors_ptvEnla,indMissing_ptvEnla,MSEptvEnla] = ...
    compareFields(displPTVretractedEnlarged,displFieldOrg,plotRes);
%% Quiver for PTVR enlarge
figure
intactVectors_ptvEnla = ~(indWrongVectors_ptvEnla | indMissing_ptvEnla);
quiver(displPTVretractedEnlarged(1).pos(intactVectors_ptvEnla,1),displPTVretractedEnlarged(1).pos(intactVectors_ptvEnla,2),...
    displPTVretractedEnlarged(1).vec(intactVectors_ptvEnla,1),displPTVretractedEnlarged(1).vec(intactVectors_ptvEnla,2),0,'Color','b')
hold on
quiver(displPTVretractedEnlarged(1).pos(indWrongVectors_ptvEnla,1),displPTVretractedEnlarged(1).pos(indWrongVectors_ptvEnla,2),...
    displPTVretractedEnlarged(1).vec(indWrongVectors_ptvEnla,1),displPTVretractedEnlarged(1).vec(indWrongVectors_ptvEnla,2),0,'Color','g')
plot(displPTVretractedEnlarged(1).pos(indMissing_ptvEnla,1),displPTVretractedEnlarged(1).pos(indMissing_ptvEnla,2),'ro')
set(gca, 'YDir','reverse')
set(gca, 'XLim',[32 480], 'YLim',[32 480])
set(gca, 'PlotBoxAspectRatio',[1 1 1])
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])
%% Plotting MSEs for PTVfiltered, PTVretracked, and PTVEnla
MSEGroup2 = {MSE_ptvFilt,MSEptvR,MSEptvEnla};
namePTVs = {'cPTV-filtered','cPTVR-median','cPTVR-with enlargement factor'};
figure, barPlotCellArray(MSEGroup2,namePTVs,1)
ylabel('Mean Squared Deviation (1)')
title('MSD in FOV')
%% map for PTVretracked
dataPath='./analysisWithoutSDC/TFMPackage/correctedDisplacementField_retrackedUpperMedian';
generateHeatmapFromField(displPTVretracted,dataPath,0,60); % step 

%% map for PTVretrackedEnla
dataPath='./analysisWithoutSDC/TFMPackage/correctedDisplacementField_enlarge7_cor_nei_med';
generateHeatmapFromField(displPTVretractedEnlarged,dataPath,0,60); % step 

%% Now Fig 5!!! Force map
%% Org force field
fStruct = load('forceGTlarger.mat');
fx = fStruct.force_x;
fy = fStruct.force_y;
sf = 1/50;
quiver(x_mat_u(xR,yR),y_mat_u(xR,yR), sf*fx(xR,yR), sf*fy(xR,yR), 0,'b')
set(gca, 'YDir','reverse')
set(gca, 'XLim',[32 480], 'YLim',[32 480])
set(gca, 'PlotBoxAspectRatio',[1 1 1])
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])
%% Org forcemap
generateHeatmapFromGridData(x_mat_u,y_mat_u,fx,fy,'./orgForceField',32,0,12000,0)

%% PIV Suite-based force
forceFieldOrg.x_mat_u = x_mat_u;
forceFieldOrg.y_mat_u = y_mat_u;
forceFieldOrg.ux = fx;
forceFieldOrg.uy = fy;
plotRes = true;
%% field for PTVretrackedEnla
forcePIVSuiteStruct= load('./analysisWithoutSDC/TFMPackage/forceField_fromPIVSuite/forceField.mat');
forceFieldSuite = forcePIVSuiteStruct.forceField;
%% force field for suite
[deviationLevelsSuiteF,indWrongVectorsSuiteF,indMissingSuiteF, MSESuiteF,DTMSuiteF] =  ...
    compareFields(forceFieldSuite,forceFieldOrg,plotRes, sf);
% intactVectorsSuiteF = ~(indWrongVectorsSuiteF | indMissingSuiteF);
% figure;
% quiver(forceFieldSuite(1).pos(intactVectorsSuiteF,1),...
%     forceFieldSuite(1).pos(intactVectorsSuiteF,2),...
%     sf*forceFieldSuite(1).vec(intactVectorsSuiteF,1),...
%     sf*forceFieldSuite(1).vec(intactVectorsSuiteF,2),0,'Color','b')
% hold on
% quiver(forceFieldSuite(1).pos(indWrongVectorsSuiteF,1),...
%     forceFieldSuite(1).pos(indWrongVectorsSuiteF,2),...
%     sf*forceFieldSuite(1).vec(indWrongVectorsSuiteF,1),...
%     sf*forceFieldSuite(1).vec(indWrongVectorsSuiteF,2),0,'Color','g')
% plot(forceFieldSuite(1).pos(indMissingSuiteF,1),...
%     forceFieldSuite(1).pos(indMissingSuiteF,2),'ro')
% set(gca, 'YDir','reverse')
% set(gca, 'XLim',[32 480], 'YLim',[32 480])
% set(gca, 'PlotBoxAspectRatio',[1 1 1])
figure;
quiver(forceFieldSuite(1).pos(:,1),...
    forceFieldSuite(1).pos(:,2),...
    sf*forceFieldSuite(1).vec(:,1),...
    sf*forceFieldSuite(1).vec(:,2),0,'Color','k')
set(gca, 'YDir','reverse')
set(gca, 'XLim',[32 480], 'YLim',[32 480])
set(gca, 'PlotBoxAspectRatio',[1 1 1])
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])
%% map for PIVSuite
dataPath='./analysisWithoutSDC/TFMPackage/forceField_fromPIVSuite';
generateHeatmapFromField(forceFieldSuite,dataPath,0,10000); % step 
%% running force calculation for Tseng's PIV result
displField = displFieldTseng;
mkdir('./analysisWithoutSDC/TFMPackage/correctedDisplacementField_Tseng')
save('./analysisWithoutSDC/TFMPackage/correctedDisplacementField_Tseng/displField.mat','displField')
%% run TFM
%% field for Tseng
forceTsengStruct= load('./analysisWithoutSDC/TFMPackage/forceField_Tseng/forceField.mat');
forceFieldTseng = forceTsengStruct.forceField;
%% force field for Tseng
[deviationLevelsTsengF,indWrongVectorsTsengF,indMissingTsengF, MSETsengF, DTMTsengF] =  ...
    compareFields(forceFieldTseng,forceFieldOrg,plotRes, sf);
%% field
figure;
quiver(forceFieldTseng(1).pos(:,1),...
    forceFieldTseng(1).pos(:,2),...
    sf*forceFieldTseng(1).vec(:,1),...
    sf*forceFieldTseng(1).vec(:,2),0,'Color','k')
set(gca, 'YDir','reverse')
set(gca, 'XLim',[32 480], 'YLim',[32 480])
set(gca, 'PlotBoxAspectRatio',[1 1 1])
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])
%% force map for Tseng
dataPath='./analysisWithoutSDC/TFMPackage/forceField_Tseng';
generateHeatmapFromField(forceFieldTseng,dataPath,0,10000); % step 

%% mpiv force
displField=displFieldMPIV;
mkdir('./analysisWithoutSDC/TFMPackage/correctedDisplacementField_mpiv')
save('./analysisWithoutSDC/TFMPackage/correctedDisplacementField_mpiv/displField.mat','displField')
%% run TFM
%% field for mpiv
forceMpivStruct= load('./analysisWithoutSDC/TFMPackage/forceField_mpiv/forceField.mat');
forceFieldMpiv = forceMpivStruct.forceField;
%% force field for mpiv
[deviationLevelsMpivF,indWrongVectorsMpivF,indMissingMpivF, MSEMpivF, DTMMpivF] =  ...
    compareFields(forceFieldMpiv,forceFieldOrg,plotRes, sf);

%% field
figure;
quiver(forceFieldMpiv(1).pos(:,1),...
    forceFieldMpiv(1).pos(:,2),...
    sf*forceFieldMpiv(1).vec(:,1),...
    sf*forceFieldMpiv(1).vec(:,2),0,'Color','k')
set(gca, 'YDir','reverse')
set(gca, 'XLim',[32 480], 'YLim',[32 480])
set(gca, 'PlotBoxAspectRatio',[1 1 1])
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])
%% force map for mpiv
dataPath='./analysisWithoutSDC/TFMPackage/forceField_mpiv';
generateHeatmapFromField(forceFieldMpiv,dataPath,0,10000); % step 
%% cptv filtering
%% field for cptv filtering
forceCPTVFStruct= load('./analysisWithoutSDC/TFMPackage/forceField_cPTVFiltering/forceField.mat');
forceFieldCPTVF = forceCPTVFStruct.forceField;
%% force field for cptv filtering
[deviationLevelsCPTVFF,indWrongVectorsCPTVFF,indMissingCPTVFF, MSECPTVFF, DTMcPTVFF] =  ...
    compareFields(forceFieldCPTVF,forceFieldOrg,plotRes, sf);
%% field
figure;
quiver(forceFieldCPTVF(1).pos(:,1),...
    forceFieldCPTVF(1).pos(:,2),...
    sf*forceFieldCPTVF(1).vec(:,1),...
    sf*forceFieldCPTVF(1).vec(:,2),0,'Color','k')
set(gca, 'YDir','reverse')
set(gca, 'XLim',[32 480], 'YLim',[32 480])
set(gca, 'PlotBoxAspectRatio',[1 1 1])
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])

%% force map for cPTV filtering
dataPath='./analysisWithoutSDC/TFMPackage/forceField_cPTVFiltering';
generateHeatmapFromField(forceFieldCPTVF,dataPath,0,10000); % step 

%% cPTVR with median
%run TFM package
displField = displPTVretracted;
save('./analysisWithoutSDC/TFMPackage/correctedDisplacementField_retrackedUpperMedian/displField.mat','displField')
%% field for cptvr med
forceCPTVR_MedStruct= load('./analysisWithoutSDC/TFMPackage/forceField_cPTVR_median/forceField.mat');
forceFieldCPTVR_Med = forceCPTVR_MedStruct.forceField;
%% force field for cptvr med
[deviationLevelsCPTVR_MedF,indWrongVectorsCPTVR_MedF,indMissingCPTVR_MedF, MSECPTVR_MedF, DTMcPTVR_MedF] =  ...
    compareFields(forceFieldCPTVR_Med,forceFieldOrg,plotRes, sf);
%% field cptvr med
figure;
quiver(forceFieldCPTVR_Med(1).pos(:,1),...
    forceFieldCPTVR_Med(1).pos(:,2),...
    sf*forceFieldCPTVR_Med(1).vec(:,1),...
    sf*forceFieldCPTVR_Med(1).vec(:,2),0,'Color','k')
set(gca, 'YDir','reverse')
set(gca, 'XLim',[32 480], 'YLim',[32 480])
set(gca, 'PlotBoxAspectRatio',[1 1 1])
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])
%% force map for cptvr med
dataPath='./analysisWithoutSDC/TFMPackage/forceField_cPTVR_median';
generateHeatmapFromField(forceFieldCPTVR_Med,dataPath,0,10000); % step 
%% CPTVR with EF
%% field for cptvr ef
forceCPTVR_EFStruct= load('./analysisWithoutSDC/TFMPackage/forceField_cPTVR_EF_reg2e-6/forceField.mat');
forceFieldCPTVR_EF = forceCPTVR_EFStruct.forceField;
%% force field for cptvr ef
[deviationLevelsCPTVR_EF,indWrongVectorsCPTVR_EF,indMissingCPTVR_EF, MSECPTVR_EF,DTMcPTVR_EF] =  ...
    compareFields(forceFieldCPTVR_EF,forceFieldOrg,plotRes, sf);
%% field cptvr ef
figure;
quiver(forceFieldCPTVR_EF(1).pos(:,1),...
    forceFieldCPTVR_EF(1).pos(:,2),...
    sf*forceFieldCPTVR_EF(1).vec(:,1),...
    sf*forceFieldCPTVR_EF(1).vec(:,2),0,'Color','k')
set(gca, 'YDir','reverse')
set(gca, 'XLim',[32 480], 'YLim',[32 480])
set(gca, 'PlotBoxAspectRatio',[1 1 1])
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])
%% force map for cptvr ef
dataPath='./analysisWithoutSDC/TFMPackage/forceField_cPTVR_EF';
generateHeatmapFromField(forceFieldCPTVR_EF,dataPath,0,10000); % step 

%% MSE plot
MSDGroup = {MSESuiteF, MSETsengF, MSEMpivF,MSECPTVFF, MSECPTVR_MedF,MSECPTVR_EF};
namePITVs = {'PIVSuite','TsengPIV','mpiv','cPTV-filtered','cPTVR-med','cPTVR-EF'};
figure, barPlotCellArray(MSDGroup,namePITVs,1)
ylabel('Mean Squared Deviation (1)')
title('MSD in FOV')
%% MSE for small and large adhesions
% We need a mask for large and small adhesions
open('./orgField/heatMap/figs/hMapFigWithMasks.fig');
h = drawellipse;
maskLargeAdh = createMask(h);
h2 = drawpolygon;
maskSmallAdh = createMask(h2);
save('./maskAdhs.mat','maskLargeAdh','maskSmallAdh')
%% Large adh first
fF_Suite_LA = filterDisplacementField( forceFieldSuite, maskLargeAdh);
[~,~,~, MSELA_f{1}] =  compareFields(fF_Suite_LA,forceFieldOrg);
fF_Tseng_LA = filterDisplacementField( forceFieldTseng, maskLargeAdh);
[~,~,~, MSELA_f{2}] =  compareFields(fF_Tseng_LA,forceFieldOrg);
fF_mpiv_LA = filterDisplacementField( forceFieldMpiv, maskLargeAdh);
[~,~,~, MSELA_f{3}] =  compareFields(fF_mpiv_LA,forceFieldOrg);
fF_ptvf_LA = filterDisplacementField( forceFieldCPTVF, maskLargeAdh);
[~,~,~, MSELA_f{4}] =  compareFields(fF_ptvf_LA,forceFieldOrg);
fF_ptvr_med_LA = filterDisplacementField( forceFieldCPTVR_Med, maskLargeAdh);
[~,~,~, MSELA_f{5}] =  compareFields(fF_ptvr_med_LA,forceFieldOrg);
fF_ptvr_ef_LA = filterDisplacementField( forceFieldCPTVR_EF, maskLargeAdh);
[~,~,~, MSELA_f{6}] =  compareFields(fF_ptvr_ef_LA,forceFieldOrg);
figure, barPlotCellArray(MSELA_f,namePITVs,1)
ylabel('Mean squared deviation (1)')
title('Force MSD at large forces')
%% Small adhs
fF_Suite_SA = filterDisplacementField( forceFieldSuite, maskSmallAdh);
[~,~,~, MSESA_f{1}] =  compareFields(fF_Suite_SA,forceFieldOrg);
fF_Tseng_SA = filterDisplacementField( forceFieldTseng, maskSmallAdh);
[~,~,~, MSESA_f{2}] =  compareFields(fF_Tseng_SA,forceFieldOrg);
fF_mpiv_SA = filterDisplacementField( forceFieldMpiv, maskSmallAdh);
[~,~,~, MSESA_f{3}] =  compareFields(fF_mpiv_SA,forceFieldOrg);
fF_ptvf_SA = filterDisplacementField( forceFieldCPTVF, maskSmallAdh);
[~,~,~, MSESA_f{4}] =  compareFields(fF_ptvf_SA,forceFieldOrg);
fF_ptvr_med_SA = filterDisplacementField( forceFieldCPTVR_Med, maskSmallAdh);
[~,~,~, MSESA_f{5}] =  compareFields(fF_ptvr_med_SA,forceFieldOrg);
fF_ptvr_ef_SA = filterDisplacementField( forceFieldCPTVR_EF, maskSmallAdh);
[~,~,~, MSESA_f{6}] =  compareFields(fF_ptvr_ef_SA,forceFieldOrg);
figure, barPlotCellArray(MSESA_f,namePITVs,1)
ylabel('Mean squared deviation (1)')
title('Force MSD at small forces')
%% save all data
save('allData.mat')