clear all
close all
load particleBehaviorAdaptiveWindows.mat
cd figuresSptVsWindows

movieName = 'Talin 120516 Cs2C3';

for iType = [1 2 9]
    
    %particle density
    plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).directAll.spDensity,...
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
    for iMode = 1 : 3
        plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).diffModeAnalysis.fracMode,...
            iMode,20,[movieName ', Fraction Mode ' num2str(iMode) ', protType ' num2str(iType)],...
            ['fractionMode' num2str(iMode) '_protType' num2str(iType) '.fig']);
    end
    
    %density of each diffusion mode
    for iMode = 1 : 3
        plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).diffModeAnalysis.modeDensity,...
            iMode,20,[movieName ', Density Mode ' num2str(iMode) ', protType ' num2str(iType)],...
            ['densityMode' num2str(iMode) '_protType' num2str(iType) '.fig']);
    end
    
    %diffusion coefficient in each diffusion mode
    for iMode = 1 : 3
        plotSptRelToActivityOnsetAdaptiveWindows(sptPropInWindow(iType).diffModeAnalysis.diffCoefModeInd,...
            iMode,20,[movieName ', Diff. coef. Mode ' num2str(iMode) ', protType ' num2str(iType)],...
            ['diffCoefMode' num2str(iMode) '_protType' num2str(iType) '.fig']);
    end
    
end

