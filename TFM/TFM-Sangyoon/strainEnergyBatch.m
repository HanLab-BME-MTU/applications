%% open necessary MLs
[fileSFolders, pathSFolders] = uigetfile('*.mat','Select selectedFolders.');
selectedFolders=load([pathSFolders filesep fileSFolders]);
pathAnalysisAll=selectedFolders.pathAnalysisAll;
numConditions = numel(pathAnalysisAll);
%% Load movieLists for each condition
for k=1:numConditions
    MLAll(k) = MovieList.load([pathAnalysisAll{k} filesep 'movieList.mat']);
end
%% Output
rootAnalysis = fileparts(pathAnalysisAll{1});
figPath = [rootAnalysis '/AnalysisSummaryTFM/Figs'];
mkdir(figPath)
dataPath = [rootAnalysis '/AnalysisSummaryTFM/Data'];
mkdir(dataPath)
%% Now Force analysis: total force, strain energy and SE density
SEGroup = cell(numConditions,1);
SEDenGroup = cell(numConditions,1);
totalForceGroup = cell(numConditions,1);
for ii=1:numConditions
    curML=MLAll(ii);
    curMovies = curML.movies_;
    N(ii) = numel(curMovies);
    SEPerGroup = cell(N(ii),1);
    SEDenPerGroup = cell(N(ii),1);
    totForcePerGroup = cell(N(ii),1);
    for k=1:N(ii)
        % get the tracksNA
        curMovie=curMovies{k};
        % get the strain energy
        iCurSEProc = curMovie.getProcessIndex('StrainEnergyCalculationProcess');
        curSEProc = curMovie.getProcess(iCurSEProc);
        seCellStruct=load(curSEProc.outFilePaths_{3});
        curSE_struct = seCellStruct.SE_Blobs;
        curSE = curSE_struct.SE;
        curSEDen = curSE_struct.SEDensity;
        forceStruct = seCellStruct.totalForceBlobs;
        curTotalForce = forceStruct.force;
%         curSE_struct = seCellStruct.SE_Cell;
%         curSE = curSE_struct.SE;
%         curSEDen = curSE_struct.SEDensity;
%         curTotalForce = seCellStruct.totalForceCell;
        
        SEPerGroup{k}=curSE;
        SEDenPerGroup{k}=curSEDen;
        totForcePerGroup{k}=curTotalForce;
    end
    SEGroup{ii,1}=SEPerGroup;
    SEDenGroup{ii,1}=SEDenPerGroup;
    totalForceGroup{ii,1}=totForcePerGroup;
    clear SEPerGroup SEDenPerGroup totForcePerGroup 
end
disp('Done')
%% setting up group name
for ii=1:numConditions
    [~, finalFolder]=fileparts(pathAnalysisAll{ii});
    groupNames{ii} = finalFolder;
end
%% Plotting each - SE
nameList=groupNames'; %{'pLVX' 'P29S'};
SEGroupCell = cellfun(@(x) cell2mat(x),SEGroup,'unif',false);
h1=figure; 
% barPlotCellArray(SEGroupCell,nameList,1)
boxPlotCellArray(SEGroupCell,nameList,1,false,true)
ylabel('Strain energy (femto-Joule)')
title('Total strain energy')
hgexport(h1,strcat(figPath,'/strainEnergyTotal'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/strainEnergyTotal'),'-v7.3')
print(h1,strcat(figPath,'/strainEnergyTotal.tif'),'-dtiff')
%% Strain energy density
SEDenGroupCell = cellfun(@(x) cell2mat(x),SEDenGroup,'unif',false);
h1=figure; 
% barPlotCellArray(SEGroupCell,nameList,1)
boxPlotCellArray(SEDenGroupCell,nameList,1,false,true)
ylabel('Strain energy density (J/m^2)')
title('Strain energy density')
hgexport(h1,strcat(figPath,'/SEDensity'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/SEDensity'),'-v7.3')
print(h1,strcat(figPath,'/SEDensity.tif'),'-dtiff')
%% Total force
totForceCell = cellfun(@(x) cell2mat(x),totalForceGroup,'unif',false);
h1=figure; 
% barPlotCellArray(SEGroupCell,nameList,1)
boxPlotCellArray(totForceCell,nameList,1,false,true)
ylabel('Total force (nN)')
title('Total force in force blobs')
hgexport(h1,strcat(figPath,'/totForce'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/totForce'),'-v7.3')
print(h1,strcat(figPath,'/totForce.tif'),'-dtiff')
%% save entire workspace for later
save([dataPath filesep 'allData.mat'])
