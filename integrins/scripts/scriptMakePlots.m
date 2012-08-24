clear all
close all
load particleBehaviorAdaptiveWindows120817.mat
mkdir figuresSptVsWindows120817
cd figuresSptVsWindows120817

movieName = 'Talin 120524 Cs1C2';

numMode = 3;

for iType = [1 2 9]
    
    %particle density
    plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).directAll.densityParticles,...
        1,20,[movieName ', Overall particle density, protType ' num2str(iType)],...
        ['densityOverall_protType' num2str(iType) '.fig']);
    
    %cell area
    plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).directAll.cellArea,...
        1,20,[movieName ', Contributing cell area per event, protType ' num2str(iType)],...
        ['cellArea_protType' num2str(iType) '.fig']);
    
    %frame-to-frame displacement
    plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).directAll.f2fDispMag2D,...
        1,20,[movieName ', Frame-to-frame disp., protType ' num2str(iType)],...
        ['f2fDisp_protType' num2str(iType) '.fig']);
    
    %overall diffusion coefficient
    plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).diffModeAnalysis.diffCoefModeAll,...
        1,20,[movieName ', Overall diff. coef., protType ' num2str(iType)],...
        ['diffCoefOverall_protType' num2str(iType) '.fig']);
    
    %fraction in each diffusion mode
    for iMode = 1 : numMode
        plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).diffModeAnalysis.fracModeClass,...
            iMode,20,[movieName ', Fraction Mode ' num2str(iMode) ', protType ' num2str(iType)],...
            ['fractionMode' num2str(iMode) '_protType' num2str(iType) '.fig']);
    end
    
    %density of each diffusion mode
    for iMode = 1 : numMode
        plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).diffModeAnalysis.densityModeClass,...
            iMode,20,[movieName ', Density Mode ' num2str(iMode) ', protType ' num2str(iType)],...
            ['densityMode' num2str(iMode) '_protType' num2str(iType) '.fig']);
    end
    
    %diffusion coefficient in each diffusion mode
    for iMode = 1 : numMode
        plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).diffModeAnalysis.diffCoefModeInd,...
            iMode,20,[movieName ', Diff. coef. Mode ' num2str(iMode) ', protType ' num2str(iType)],...
            ['diffCoefMode' num2str(iMode) '_protType' num2str(iType) '.fig']);
    end
    
    %fraction of merges
    plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).mergesSplits.fracMerge,...
        1,20,[movieName ', Fraction merges, protType ' num2str(iType)],...
        ['fracMerge_protType' num2str(iType) '.fig']);
    
    %fraction of splits
    plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).mergesSplits.fracSplit,...
        1,20,[movieName ', Fraction splits, protType ' num2str(iType)],...
        ['fracSplit_protType' num2str(iType) '.fig']);
    
end

