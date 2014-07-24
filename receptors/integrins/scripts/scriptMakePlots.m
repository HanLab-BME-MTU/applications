clear all
close all
load particleBehaviorAdaptiveWindows130219.mat
mkdir figuresSptVsWindows
cd figuresSptVsWindows

movieName = 'Farn Sim 01';

numMode = 4;

for iType = 1
    
    %particle density
    plotSptRelToActivityOnsetAdaptiveWindowsV2(sptPropInWindow(iType).directAll.densityParticles,...
        windowDistFromEdge(iType),1,20,[movieName ', Molecule density, protType ' num2str(iType)],...
        [movieName ' DensityOverall_protType' num2str(iType) '.fig'],[0.111 10 1/(0.111^2)],...
        'Molecule density (per um^2)',[-100 150 0 11],1,[movieName ', Edge & windows, protType ' num2str(iType)],...
            [movieName ' EdgeWindows_protType' num2str(iType) '.fig']);
    
    %overall diffusion coefficient
    plotSptRelToActivityOnsetAdaptiveWindowsV2(sptPropInWindow(iType).diffModeAnalysis.diffCoefModeAll,...
        windowDistFromEdge(iType),1,20,[movieName ', Diffusion coefficient, protType ' num2str(iType)],...
        [movieName ' DiffCoefOverall_protType' num2str(iType) '.fig'],[0.111 10 (0.111^2)/0.025],'Diffusion coefficient (um^2/s)');
    
    %     %fraction in each diffusion mode
    %     for iMode = 1 : numMode
    %         plotSptRelToActivityOnsetAdaptiveWindowsV2(sptPropInWindow(iType).diffModeAnalysis.fracModeClass,...
    %             windowDistFromEdge(iType),iMode,20,[movieName ', Fraction Mode ' num2str(iMode) ', protType ' num2str(iType)],...
    %             [movieName ' FractionMode' num2str(iMode) '_protType' num2str(iType) '.fig'],[0.111 10 1],['Fraction Mode ' num2str(iMode)]);
    %     end
    %
    %     %density of each diffusion mode
    %     for iMode = 1 : numMode
    %         plotSptRelToActivityOnsetAdaptiveWindowsV2(sptPropInWindow(iType).diffModeAnalysis.densityModeClass,...
    %             windowDistFromEdge(iType),iMode,20,[movieName ', Density Mode ' num2str(iMode) ', protType ' num2str(iType)],...
    %             [movieName ' DensityMode' num2str(iMode) '_protType' num2str(iType) '.fig'],[0.111 10 1/(0.111^2)],['Density Mode ' num2str(iMode) ' (um^-^2) ']);
    %     end
    
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

