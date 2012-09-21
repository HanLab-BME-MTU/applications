clear all
close all
load particleBehaviorAdaptiveWindows120921.mat
mkdir figuresSptVsWindows120921
cd figuresSptVsWindows120921

movieName = 'AlphaV 110909 Cs2C1';

numMode = 4;

for iType = [1 2]
    
    %particle density
    plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).directAll.densityParticles,...
        windowDistFromEdge(iType),1,20,[movieName ', Particle density overall, protType ' num2str(iType)],...
        [movieName ' DensityOverall_protType' num2str(iType) '.fig'],[0.111 10 1/(0.111^2)],'Density overall (um^-^2)');
    
    %overall diffusion coefficient
    plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).diffModeAnalysis.diffCoefModeAll,...
        windowDistFromEdge(iType),1,20,[movieName ', Diffusion coefficient overall, protType ' num2str(iType)],...
        [movieName ' DiffCoefOverall_protType' num2str(iType) '.fig'],[0.111 10 (0.111^2)/0.025],'Diffusion coefficient overall (um^2/s)');
    
    %fraction in each diffusion mode
    for iMode = 1 : numMode
        plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).diffModeAnalysis.fracModeClass,...
            windowDistFromEdge(iType),iMode,20,[movieName ', Fraction Mode ' num2str(iMode) ', protType ' num2str(iType)],...
            [movieName ' FractionMode' num2str(iMode) '_protType' num2str(iType) '.fig'],[0.111 10 1],['Fraction Mode ' num2str(iMode)]);
    end
    
    %density of each diffusion mode
    for iMode = 1 : numMode
        plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).diffModeAnalysis.densityModeClass,...
            windowDistFromEdge(iType),iMode,20,[movieName ', Density Mode ' num2str(iMode) ', protType ' num2str(iType)],...
            [movieName ' DensityMode' num2str(iMode) '_protType' num2str(iType) '.fig'],[0.111 10 1/(0.111^2)],['Density Mode ' num2str(iMode) ' (um^-^2) ']);
    end
    
end

%     %cell area
%     plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).directAll.cellArea,...
%         1,20,[movieName ', Contributing cell area per event, protType ' num2str(iType)],...
%         ['cellArea_protType' num2str(iType) '.fig']);

%     %frame-to-frame displacement
%     plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).directAll.f2fDispMag2D,...
%         1,20,[movieName ', Frame-to-frame disp., protType ' num2str(iType)],...
%         ['f2fDisp_protType' num2str(iType) '.fig']);

%     %diffusion coefficient in each diffusion mode
%     for iMode = 1 : numMode
%         plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).diffModeAnalysis.diffCoefModeInd,...
%             iMode,20,[movieName ', Diff. coef. Mode ' num2str(iMode) ', protType ' num2str(iType)],...
%             ['diffCoefMode' num2str(iMode) '_protType' num2str(iType) '.fig']);
%     end

%     %fraction of merges
%     plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).mergesSplits.fracMerge,...
%         1,20,[movieName ', Fraction merges, protType ' num2str(iType)],...
%         ['fracMerge_protType' num2str(iType) '.fig']);

%     %fraction of splits
%     plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).mergesSplits.fracSplit,...
%         1,20,[movieName ', Fraction splits, protType ' num2str(iType)],...
%         ['fracSplit_protType' num2str(iType) '.fig']);

