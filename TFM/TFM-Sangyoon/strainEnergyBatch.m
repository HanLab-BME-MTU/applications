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
N=zeros(numConditions,1);
SE_FB_Group = cell(numConditions,1);
SEDen_FB_Group = cell(numConditions,1);
totalForce_FB_Group = cell(numConditions,1);
SE_Cell_Group = cell(numConditions,1);
SEDen_Cell_Group = cell(numConditions,1);
totalForce_Cell_Group = cell(numConditions,1);
SE_FOV_Group = cell(numConditions,1);
SEDen_FOV_Group = cell(numConditions,1);
totalForce_FOV_Group = cell(numConditions,1);
for ii=1:numConditions
    curML=MLAll(ii);
    curMovies = curML.movies_;
    N(ii) = numel(curMovies);
    SEPerFBGroup = cell(N(ii),1);
    SEDenPerFBGroup = cell(N(ii),1);
    totForcePerFBGroup = cell(N(ii),1);
    SEPerCellGroup = cell(N(ii),1);
    SEDenPerCellGroup = cell(N(ii),1);
    totForcePerCellGroup = cell(N(ii),1);
    SEPerFOVGroup = cell(N(ii),1);
    SEDenPerFOVGroup = cell(N(ii),1);
    totForcePerFOVGroup = cell(N(ii),1);
    for k=1:N(ii)
        % get the tracksNA
        curMovie=curMovies{k};
        % get the strain energy
        iCurSEProc = curMovie.getProcessIndex('StrainEnergyCalculationProcess');
        curSEProc = curMovie.getProcess(iCurSEProc);
        % There are three types: FOV(1). cell(2). and forceBlobs(3)
        % 3. ForceBlob
        seFBStruct=load(curSEProc.outFilePaths_{3});
        curSEFB_struct = seFBStruct.SE_Blobs;
        curSEFB = curSEFB_struct.SE;
        curSEFBDen = curSEFB_struct.SEDensity;
        forceFBStruct = seFBStruct.totalForceBlobs;
        curTotalForceFB = forceFBStruct.force;
        
        SEPerFBGroup{k}=curSEFB;
        SEDenPerFBGroup{k}=curSEFBDen;
        totForcePerFBGroup{k}=curTotalForceFB;
        
        % 2. Cell - It is possible this is zero (when cell segmentation is
        % not there)
        seCellStruct=load(curSEProc.outFilePaths_{2});
        curSECell_struct = seCellStruct.SE_Cell;
        curSECell = curSECell_struct.SE;
        curSEDenCell = curSECell_struct.SEDensity;
        curTotalForceCell = seCellStruct.totalForceCell;

        SEPerCellGroup{k}=curSECell;
        SEDenPerCellGroup{k}=curSEDenCell;
        totForcePerCellGroup{k}=curTotalForceCell;
        
        % 3. FOV
        seFOVStruct=load(curSEProc.outFilePaths_{1});
        curSEFOV_struct = seFOVStruct.SE_FOV;
        curSEFOV = curSEFOV_struct.SE;
        curSEDenFOV = curSEFOV_struct.SEDensity;
        curTotalForceFOV = seFOVStruct.totalForceFOV;

        SEPerFOVGroup{k}=curSEFOV;
        SEDenPerFOVGroup{k}=curSEDenFOV;
        totForcePerFOVGroup{k}=curTotalForceFOV;
    end
    SE_FB_Group{ii,1}=SEPerFBGroup;
    SEDen_FB_Group{ii,1}=SEDenPerFBGroup;
    totalForce_FB_Group{ii,1}=totForcePerFBGroup;
    SE_Cell_Group{ii,1}=SEPerCellGroup;
    SEDen_Cell_Group{ii,1}=SEDenPerCellGroup;
    totalForce_Cell_Group{ii,1}=totForcePerCellGroup;
    SE_FOV_Group{ii,1}=SEPerFOVGroup;
    SEDen_FOV_Group{ii,1}=SEDenPerFOVGroup;
    totalForce_FOV_Group{ii,1}=totForcePerFOVGroup;
end
disp('Done')
%% setting up group name
for ii=1:numConditions
    [~, finalFolder]=fileparts(pathAnalysisAll{ii});
    groupNames{ii} = finalFolder;
end
nameList=groupNames'; 
%% Plotting each - SE-ForceBlob
SE_FB_GroupCellArray = cellfun(@(x) cell2mat(x),SE_FB_Group,'unif',false);
h1=figure; 
% barPlotCellArray(SEGroupCell,nameList,1)
boxPlotCellArray(SE_FB_GroupCellArray,nameList,1,false,true)
ylabel('Strain energy (femto-Joule)')
title('Total strain energy in force blobs')
hgexport(h1,strcat(figPath,'/strainEnergyForceBlobs'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/strainEnergyForceBlobs'),'-v7.3')
print(h1,strcat(figPath,'/strainEnergyForceBlobs.tif'),'-dtiff')

tableSE_FB=table(SE_FB_GroupCellArray,'RowNames',nameList);
writetable(tableSE_FB,strcat(dataPath,'/strainEnergyForceBlobs.csv'))
%% Plotting each - SE-Cell
SE_Cell_GroupCellArray = cellfun(@(x) cell2mat(x),SE_Cell_Group,'unif',false);
h1=figure; 
% barPlotCellArray(SEGroupCell,nameList,1)
boxPlotCellArray(SE_Cell_GroupCellArray,nameList,1,false,true)
ylabel('Strain energy (femto-Joule)')
title('Total strain energy in a cell')
hgexport(h1,strcat(figPath,'/strainEnergyCell'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/strainEnergyCell'),'-v7.3')
print(h1,strcat(figPath,'/strainEnergyCell.tif'),'-dtiff')

tableSE_Cell=table(SE_Cell_GroupCellArray,'RowNames',nameList);
writetable(tableSE_Cell,strcat(dataPath,'/strainEnergyCell.csv'))
%% Plotting each - SE-FOV
SE_FOV_GroupCellArray = cellfun(@(x) cell2mat(x),SE_FOV_Group,'unif',false);
h1=figure; 
% barPlotCellArray(SEGroupCell,nameList,1)
boxPlotCellArray(SE_FOV_GroupCellArray,nameList,1,false,true)
ylabel('Strain energy (femto-Joule)')
title('Total strain energy in a field of view')
hgexport(h1,strcat(figPath,'/strainEnergyFOV'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/strainEnergyFOV'),'-v7.3')
print(h1,strcat(figPath,'/strainEnergyFOV.tif'),'-dtiff')

tableSE_FOV=table(SE_FOV_GroupCellArray,'RowNames',nameList);
writetable(tableSE_FOV,strcat(dataPath,'/strainEnergyFOV.csv'))
%% Strain energy density - ForceBlob
SEDenGroupFB_CellArray = cellfun(@(x) cell2mat(x),SEDen_FB_Group,'unif',false);
h1=figure; 
% barPlotCellArray(SEGroupCell,nameList,1)
boxPlotCellArray(SEDenGroupFB_CellArray,nameList,1,false,true)
ylabel('Strain energy density (J/m^2)')
title('Strain energy density in force blobs')
hgexport(h1,strcat(figPath,'/SEDensityFB'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/SEDensityFB'),'-v7.3')
print(h1,strcat(figPath,'/SEDensityFB.tif'),'-dtiff')

tableSED_FB=table(SEDenGroupFB_CellArray,'RowNames',nameList);
writetable(tableSED_FB,strcat(dataPath,'/strainEnergyDensityForceBlobs.csv'))
%% Strain energy density - Cell
SEDenGroupCell_CellArray = cellfun(@(x) cell2mat(x),SEDen_Cell_Group,'unif',false);
h1=figure; 
% barPlotCellArray(SEGroupCell,nameList,1)
boxPlotCellArray(SEDenGroupCell_CellArray,nameList,1,false,true)
ylabel('Strain energy density (J/m^2)')
title('Strain energy density in a cell')
hgexport(h1,strcat(figPath,'/SEDensityCell'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/SEDensityCell'),'-v7.3')
print(h1,strcat(figPath,'/SEDensityCell.tif'),'-dtiff')

tableSED_Cell=table(SEDenGroupCell_CellArray,'RowNames',nameList);
writetable(tableSED_Cell,strcat(dataPath,'/strainEnergyDensityCell.csv'))
%% Strain energy density - FOV
SEDenGroupFOV_CellArray = cellfun(@(x) cell2mat(x),SEDen_FOV_Group,'unif',false);
h1=figure; 
% barPlotCellArray(SEGroupCell,nameList,1)
boxPlotCellArray(SEDenGroupFOV_CellArray,nameList,1,false,true)
ylabel('Strain energy density (J/m^2)')
title('Strain energy density in a field of view')
hgexport(h1,strcat(figPath,'/SEDensityFOV'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/SEDensityFOV'),'-v7.3')
print(h1,strcat(figPath,'/SEDensityFOV.tif'),'-dtiff')

tableSED_FOV=table(SEDenGroupFOV_CellArray,'RowNames',nameList);
writetable(tableSED_FOV,strcat(dataPath,'/strainEnergyDensityFOV.csv'))
%% Total force - Force Blobs
totForceFBCellArray = cellfun(@(x) cell2mat(x),totalForce_FB_Group,'unif',false);
h1=figure; 
boxPlotCellArray(totForceFBCellArray,nameList,1,false,true)
ylabel('Total force (nN)')
title('Total force in force blobs')
hgexport(h1,strcat(figPath,'/totForceFB'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/totForceFB'),'-v7.3')
print(h1,strcat(figPath,'/totForceFB.tif'),'-dtiff')

tableTotalForce_FB=table(totForceFBCellArray,'RowNames',nameList);
writetable(tableTotalForce_FB,strcat(dataPath,'/totalForce_ForceBlobs.csv'))
%% Total force - Cell
totForceCell_CellArray = cellfun(@(x) cell2mat(x),totalForce_Cell_Group,'unif',false);
h1=figure; 
boxPlotCellArray(totForceCell_CellArray,nameList,1,false,true)
ylabel('Total force (nN)')
title('Total force in a cell')
hgexport(h1,strcat(figPath,'/totForceCell'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/totForceCell'),'-v7.3')
print(h1,strcat(figPath,'/totForceCell.tif'),'-dtiff')

tableTotalForce_Cell=table(totForceCell_CellArray,'RowNames',nameList);
writetable(tableTotalForce_Cell,strcat(dataPath,'/totalForce_Cell.csv'))
%% Total force - FOV
totForceFOV_CellArray = cellfun(@(x) cell2mat(x),totalForce_FOV_Group,'unif',false);
h1=figure; 
boxPlotCellArray(totForceFOV_CellArray,nameList,1,false,true)
ylabel('Total force (nN)')
title('Total force in a field of view')
hgexport(h1,strcat(figPath,'/totForceFOV'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/totForceFOV'),'-v7.3')
print(h1,strcat(figPath,'/totForceFOV.tif'),'-dtiff')

tableTotalForce_FOV=table(totForceFOV_CellArray,'RowNames',nameList);
writetable(tableTotalForce_FOV,strcat(dataPath,'/totalForce_FOV.csv'))
%% save entire workspace for later
save([dataPath filesep 'allData.mat'])
