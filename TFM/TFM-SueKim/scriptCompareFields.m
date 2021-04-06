% Sue's code for saving all publishable figures
% This is to compare the PIVs and PTVs with ground-truth filed 
% written by Sangyoon Han, March 2021
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
[deviationLevels_ptvFilt,indWrongVectors_ptvFilt,indMissing_ptvFilt] = ...
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

