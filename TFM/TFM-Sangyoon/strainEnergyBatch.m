%% open necessary MLs
[fileSFolders, pathSFolders] = uigetfile('*.mat','Select selectedFolders.mat.  If do not have one, click cancel');
if ~ischar(pathSFolders) && pathSFolders==0
    analysisFolderSelectionDone = false;
    ii=0;
    rootFolder=pwd;
    while ~analysisFolderSelectionDone
        ii=ii+1;
        curPathProject = uigetdir(rootFolder,'Select each analysis folder that contains movieList.mat (Click Cancel when no more)');
        if ~ischar(curPathProject) && curPathProject==0
            analysisFolderSelectionDone=true;
        else
            pathAnalysisAll{ii} = curPathProject;
        end
    end
else
    selectedFolders=load([pathSFolders filesep fileSFolders]);
    pathAnalysisAll=selectedFolders.pathAnalysisAll;
end
%% Load movieLists for each condition
numConditions = numel(pathAnalysisAll);
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

SE_CellPeri_Group = cell(numConditions,1);
SE_CellInside_Group = cell(numConditions,1);
totalForce_CellPeri_Group = cell(numConditions,1);
totalForce_CellInside_Group = cell(numConditions,1);
spreadArea_Group = cell(numConditions,1);
SEDen_CellPeri_Group = cell(numConditions,1);
SEDen_CellInside_Group = cell(numConditions,1);

for ii=1:numConditions
    curML=MLAll(ii);
    curMovies = curML.movies_;
    N(ii) = numel(curMovies);
    curSEPerFBGroup = cell(N(ii),1);
    curSEDenPerFBGroup = cell(N(ii),1);
    curTotForcePerFBGroup = cell(N(ii),1);
    curSECellGroup = cell(N(ii),1);
    curSEDenCellGroup = cell(N(ii),1);
    curTotalForceCellGroup = cell(N(ii),1);
    curSEPerFOVGroup = cell(N(ii),1);
    curSEDenPerFOVGroup = cell(N(ii),1);
    curTotForcePerFOVGroup = cell(N(ii),1);
    
    curSECellPeriGroup = cell(N(ii),1);
    curSECellInsideGroup = cell(N(ii),1);
    curTotalForceCellPeriGroup = cell(N(ii),1);
    curTotalForceCellInsideGroup = cell(N(ii),1);
    curSpreadAreaGroup = cell(N(ii),1);
    curSEDenCellPeriGroup = cell(N(ii),1);
    curSEDenCellInsideGroup = cell(N(ii),1);
    
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
        
        curSEPerFBGroup{k}=curSEFB;
        curSEDenPerFBGroup{k}=curSEFBDen;
        curTotForcePerFBGroup{k}=curTotalForceFB;
        
        % 2. Cell - It is possible this is zero (when cell segmentation is
        % not there)
        seCellStruct=load(curSEProc.outFilePaths_{2});
        curSECell_struct = seCellStruct.SE_Cell;
        curSECellGroup{k} = curSECell_struct.SE; %in femto-Joule=1e15*(N*m)
        curSEDenCellGroup{k} = curSECell_struct.SEDensity; % in J/m2
        curTotalForceCellGroup{k} = seCellStruct.totalForceCell; % in nN
        
        curSECellPeriGroup{k} = curSECell_struct.SE_peri; %in femto-Joule=1e15*(N*m)
        curSECellInsideGroup{k} = curSECell_struct.SE_inside; %in femto-Joule=1e15*(N*m)
        curSEDenCellPeriGroup{k} = curSECell_struct.SEDensityPeri; % in J/m2
        curSEDenCellInsideGroup{k} = curSECell_struct.SEDensityInside; % in J/m2
        curSpreadAreaGroup{k} = curSECell_struct.area; % in um2
        curTotalForceCellPeriGroup{k} = seCellStruct.totalForceCellPeri; % in nN
        curTotalForceCellInsideGroup{k} = seCellStruct.totalForceCellInside; % in nN
        
        % 3. FOV
        seFOVStruct=load(curSEProc.outFilePaths_{1});
        curSEFOV_struct = seFOVStruct.SE_FOV;
        curSEFOV = curSEFOV_struct.SE;
        curSEDenFOV = curSEFOV_struct.SEDensity;
        curTotalForceFOV = seFOVStruct.totalForceFOV;

        curSEPerFOVGroup{k}=curSEFOV;
        curSEDenPerFOVGroup{k}=curSEDenFOV;
        curTotForcePerFOVGroup{k}=curTotalForceFOV;
    end
    SE_FB_Group{ii,1}=curSEPerFBGroup;
    SEDen_FB_Group{ii,1}=curSEDenPerFBGroup;
    totalForce_FB_Group{ii,1}=curTotForcePerFBGroup;
    SE_Cell_Group{ii,1}=curSECellGroup;
    SEDen_Cell_Group{ii,1}=curSEDenCellGroup;
    totalForce_Cell_Group{ii,1}=curTotalForceCellGroup;
    SE_FOV_Group{ii,1}=curSEPerFOVGroup;
    SEDen_FOV_Group{ii,1}=curSEDenPerFOVGroup;
    totalForce_FOV_Group{ii,1}=curTotForcePerFOVGroup;
    
    SE_CellPeri_Group{ii,1}=curSECellPeriGroup;
    SE_CellInside_Group{ii,1}=curSECellInsideGroup;
    totalForce_CellPeri_Group{ii,1}=curTotalForceCellPeriGroup;
    totalForce_CellInside_Group{ii,1}=curTotalForceCellInsideGroup;
    spreadArea_Group{ii,1}=curSpreadAreaGroup;
    SEDen_CellPeri_Group{ii,1}=curSEDenCellPeriGroup;
    SEDen_CellInside_Group{ii,1}=curSEDenCellInsideGroup;
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
writetable(tableSE_FB,strcat(dataPath,'/strainEnergyForceBlobs.csv'),'WriteRowNames',true)
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
writetable(tableSE_Cell,strcat(dataPath,'/strainEnergyCell.csv'),'WriteRowNames',true)
%% Plotting each - SE-Cell-Peri
SE_CellPeri_GroupCellArray = cellfun(@(x) cell2mat(x),SE_CellPeri_Group,'unif',false);
h1=figure; 
% barPlotCellArray(SEGroupCell,nameList,1)
boxPlotCellArray(SE_CellPeri_GroupCellArray,nameList,1,false,true)
ylabel('Strain energy (femto-Joule)')
title('Total strain energy in a cell periphery')
hgexport(h1,strcat(figPath,'/strainEnergyCellPeri'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/strainEnergyCellPeri'),'-v7.3')
print(h1,strcat(figPath,'/strainEnergyCellPeri.tif'),'-dtiff')

tableSE_CellPeri=table(SE_CellPeri_GroupCellArray,'RowNames',nameList);
writetable(tableSE_CellPeri,strcat(dataPath,'/strainEnergyCellPeri.csv'),'WriteRowNames',true)
%% Plotting each - SE-Cell-Inside
SE_CellInside_GroupCellArray = cellfun(@(x) cell2mat(x),SE_CellInside_Group,'unif',false);
h1=figure; 
% barPlotCellArray(SEGroupCell,nameList,1)
boxPlotCellArray(SE_CellInside_GroupCellArray,nameList,1,false,true)
ylabel('Strain energy (femto-Joule)')
title('Total strain energy in a cell interior')
hgexport(h1,strcat(figPath,'/strainEnergyCellInside'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/strainEnergyCellInside'))
print(h1,strcat(figPath,'/strainEnergyCellInside.tif'),'-dtiff')

tableSE_CellInside=table(SE_CellInside_GroupCellArray,'RowNames',nameList);
writetable(tableSE_CellInside,strcat(dataPath,'/strainEnergyCellInside.csv'),'WriteRowNames',true)
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
writetable(tableSE_FOV,strcat(dataPath,'/strainEnergyFOV.csv'),'WriteRowNames',true)
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
writetable(tableSED_FB,strcat(dataPath,'/strainEnergyDensityForceBlobs.csv'),'WriteRowNames',true)
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
writetable(tableSED_Cell,strcat(dataPath,'/strainEnergyDensityCell.csv'),'WriteRowNames',true)
%% Strain energy density - CellPeri
SEDenGroupCellPeri_CellArray = cellfun(@(x) cell2mat(x),SEDen_CellPeri_Group,'unif',false);
h1=figure; 
% barPlotCellArray(SEGroupCell,nameList,1)
boxPlotCellArray(SEDenGroupCellPeri_CellArray,nameList,1,false,true)
ylabel('Strain energy density (J/m^2)')
title('Strain energy density in a cell periphery')
hgexport(h1,strcat(figPath,'/SEDensityCellPeri'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/SEDensityCellPeri'),'-v7.3')
print(h1,strcat(figPath,'/SEDensityCellPeri.tif'),'-dtiff')

tableSED_CellPeri=table(SEDenGroupCellPeri_CellArray,'RowNames',nameList);
writetable(tableSED_CellPeri,strcat(dataPath,'/strainEnergyDensityCellPeri.csv'),'WriteRowNames',true)
%% Strain energy density - CellInside
SEDenGroupCellInside_CellArray = cellfun(@(x) cell2mat(x),SEDen_CellInside_Group,'unif',false);
h1=figure; 
% barPlotCellArray(SEGroupCell,nameList,1)
boxPlotCellArray(SEDenGroupCellInside_CellArray,nameList,1,false,true)
ylabel('Strain energy density (J/m^2)')
title('Strain energy density in a cell interior')
hgexport(h1,strcat(figPath,'/SEDensityCellInside'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/SEDensityCellInside'),'-v7.3')
print(h1,strcat(figPath,'/SEDensityCellInside.tif'),'-dtiff')

tableSED_CellInside=table(SEDenGroupCellInside_CellArray,'RowNames',nameList);
writetable(tableSED_CellInside,strcat(dataPath,'/strainEnergyDensityCellInside.csv'),'WriteRowNames',true)
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
writetable(tableSED_FOV,strcat(dataPath,'/strainEnergyDensityFOV.csv'),'WriteRowNames',true)
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
writetable(tableTotalForce_FB,strcat(dataPath,'/totalForce_ForceBlobs.csv'),'WriteRowNames',true)
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
writetable(tableTotalForce_Cell,strcat(dataPath,'/totalForce_Cell.csv'),'WriteRowNames',true)
%% Total force - CellPeri
totForceCellPeri_CellArray = cellfun(@(x) cell2mat(x),totalForce_CellPeri_Group,'unif',false);
h1=figure; 
boxPlotCellArray(totForceCellPeri_CellArray,nameList,1,false,true)
ylabel('Total force (nN)')
title('Total force in a cell periphery')
hgexport(h1,strcat(figPath,'/totForceCellPeri'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/totForceCellPeri'),'-v7.3')
print(h1,strcat(figPath,'/totForceCellPeri.tif'),'-dtiff')

tableTotalForce_CellPeri=table(totForceCellPeri_CellArray,'RowNames',nameList);
writetable(tableTotalForce_CellPeri,strcat(dataPath,'/totalForce_CellPeri.csv'),'WriteRowNames',true)
%% Total force - CellInside
totForceCellInside_CellArray = cellfun(@(x) cell2mat(x),totalForce_CellInside_Group,'unif',false);
h1=figure; 
boxPlotCellArray(totForceCellInside_CellArray,nameList,1,false,true)
ylabel('Total force (nN)')
title('Total force in a cell interior')
hgexport(h1,strcat(figPath,'/totForceCellInside'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/totForceCellInside'),'-v7.3')
print(h1,strcat(figPath,'/totForceCellInside.tif'),'-dtiff')

tableTotalForce_CellInside=table(totForceCellInside_CellArray,'RowNames',nameList);
writetable(tableTotalForce_CellInside,strcat(dataPath,'/totalForce_CellInside.csv'),'WriteRowNames',true)
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
%% Total force - Spread Area
spreadArea_CellArray = cellfun(@(x) cell2mat(x),spreadArea_Group,'unif',false);
h1=figure; 
boxPlotCellArray(spreadArea_CellArray,nameList,1,false,true)
ylabel('Spread area (um2)')
title('Cell spread area')
hgexport(h1,strcat(figPath,'/spreadArea'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/spreadArea'),'-v7.3')
print(h1,strcat(figPath,'/spreadArea.tif'),'-dtiff')

tableSpreadArea=table(spreadArea_CellArray,'RowNames',nameList);
writetable(tableSpreadArea,strcat(dataPath,'/spreadArea.csv'),'WriteRowNames',true)
%% save entire workspace for later
save([dataPath filesep 'allData.mat'])
