%% open necessary MLs
MLdirect=false;
%% Read selectedFolders or MLs directly
[pathAnalysisAll, MLNames, groupNames,usedSelectedFoldersMat, specificName] = chooseSelectedFolders;
%% Load movieLists for each condition
numConditions = numel(pathAnalysisAll);
for k=1:numConditions
    MLAll(k) = MovieList.load([pathAnalysisAll{k} filesep MLNames{k}]);
end
%% Output
rootAnalysis = pathAnalysisAll{1};
figPath = [rootAnalysis '/AnalysisSummary_TFM' specificName '/Figs'];
mkdir(figPath)
dataPath = [rootAnalysis '/AnalysisSummary_TFM' specificName '/Data'];
mkdir(dataPath)
%% Now Force analysis: total force, strain energy and SE density
N=zeros(numConditions,1);
SE_FB_Group = cell(numConditions,1);
SEDen_FB_Group = cell(numConditions,1);
totalForce_FB_Group = cell(numConditions,1);
SE_Cell_Group = cell(numConditions,1);
SEDen_Cell_Group = cell(numConditions,1);
totalForce_Cell_Group = cell(numConditions,1);
totalDisp_Cell_Group = cell(numConditions,1);
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

isCellSeg = true;

for ii=1:numConditions
    curML=MLAll(ii);
    curMovies = curML.movies_;
    N(ii) = numel(curMovies);
    curSEPerFBGroup = cell(N(ii),1);
    curSEDenPerFBGroup = cell(N(ii),1);
    curTotForceAllFBGroup = cell(N(ii),1);
    curSECellGroup = cell(N(ii),1);
    curSEDenCellGroup = cell(N(ii),1);
    curTotalForceCellGroup = cell(N(ii),1);
    curTotalDispCellGroup = cell(N(ii),1);
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
        if isfield(forceFBStruct,'avgTractionCell')
            curForcePerBlobs = forceFBStruct.avgTractionCell;
            curTotForcePerFBGroup{k}=cell2mat(curForcePerBlobs);
        end
        
        curSEPerFBGroup{k}=curSEFB;
        curSEDenPerFBGroup{k}=curSEFBDen;
        curTotForceAllFBGroup{k}=curTotalForceFB;
        
        % 2. Cell - It is possible this is zero (when cell segmentation is
        % not there)
        if exist(curSEProc.outFilePaths_{2},'file')
            seCellStruct=load(curSEProc.outFilePaths_{2});
            curSECell_struct = seCellStruct.SE_Cell;
            curSECellGroup{k} = curSECell_struct.SE; %in femto-Joule=1e15*(N*m)
            curSEDenCellGroup{k} = curSECell_struct.SEDensity; % in J/m2
            curTotalForceCellGroup{k} = seCellStruct.totalForceCell; % in nN
            curTotalDispCellGroup{k} = seCellStruct.totDispCell; % in um2
        
            curSECellPeriGroup{k} = curSECell_struct.SE_peri; %in femto-Joule=1e15*(N*m)
            curSECellInsideGroup{k} = curSECell_struct.SE_inside; %in femto-Joule=1e15*(N*m)
            curSEDenCellPeriGroup{k} = curSECell_struct.SEDensityPeri; % in J/m2
            curSEDenCellInsideGroup{k} = curSECell_struct.SEDensityInside; % in J/m2
            curSpreadAreaGroup{k} = curSECell_struct.area; % in um2
            curTotalForceCellPeriGroup{k} = seCellStruct.totalForceCellPeri; % in nN
            curTotalForceCellInsideGroup{k} = seCellStruct.totalForceCellInside; % in nN
        else
            isCellSeg = false;
        end
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
    totalForce_FB_Group{ii,1}=curTotForceAllFBGroup;
    if isfield(forceFBStruct,'avgTractionCell')
        avgForce_FBindiv_Group{ii,1}=cell2mat(curTotForcePerFBGroup');
        avgForce_FB_Group{ii,1}=cellfun(@(x) nanmean(x),curTotForcePerFBGroup);
        totForce_FBInCell_Group{ii,1}=cellfun(@(x) nansum(x),curTotForcePerFBGroup);
    else
        avgForce_FBindiv_Group{ii,1}=NaN;
    end
    SE_FOV_Group{ii,1}=curSEPerFOVGroup;
    SEDen_FOV_Group{ii,1}=curSEDenPerFOVGroup;
    totalForce_FOV_Group{ii,1}=curTotForcePerFOVGroup;
    
    if isCellSeg
        SE_Cell_Group{ii,1}=curSECellGroup;
        SEDen_Cell_Group{ii,1}=curSEDenCellGroup;
        totalForce_Cell_Group{ii,1}=curTotalForceCellGroup;
        totalDisp_Cell_Group{ii,1}=curTotalDispCellGroup;

        SE_CellPeri_Group{ii,1}=curSECellPeriGroup;
        SE_CellInside_Group{ii,1}=curSECellInsideGroup;
        totalForce_CellPeri_Group{ii,1}=curTotalForceCellPeriGroup;
        totalForce_CellInside_Group{ii,1}=curTotalForceCellInsideGroup;
        spreadArea_Group{ii,1}=curSpreadAreaGroup;
        SEDen_CellPeri_Group{ii,1}=curSEDenCellPeriGroup;
        SEDen_CellInside_Group{ii,1}=curSEDenCellInsideGroup;
    end
end
disp('Done')
%% setting up group name
% if exist('analysisFolderSelectionDone','var')
    if usedSelectedFoldersMat
        groupNames2=groupNames;
        for ii=1:numConditions
            [~, finalFolder]=fileparts(pathAnalysisAll{ii});
            groupNames{ii} = finalFolder;
        end
        nameList=groupNames'; 
        if any(cellfun(@isempty,nameList))
            nameList = MLNames;
            if strcmp(nameList{1},nameList{2})
                for ii=1:numConditions
                    curPath=fileparts(pathAnalysisAll{ii});
                    [~,nameList{ii}] = fileparts(curPath);
                end
            end
        end
    else
        nameList=groupNames'; 
    end
% else
%     nameList=groupNames'; 
% end
%% Plotting each - SE-ForceBlob
SE_FB_GroupCellArray = cellfun(@(x) cell2mat(x),SE_FB_Group,'unif',false);
h1=figure; 
% barPlotCellArray(SEGroupCell,nameList,1)
try
    boxPlotCellArray(SE_FB_GroupCellArray,nameList,1,false,true)
catch
    nameList = nameList';
    boxPlotCellArray(SE_FB_GroupCellArray,nameList,1,false,true)
end
ylabel('Strain energy (femto-Joule)')
title('Total strain energy in force blobs')
hgexport(h1,strcat(figPath,'/strainEnergyForceBlobs'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/strainEnergyForceBlobs'),'-v7.3')
print(h1,strcat(figPath,'/strainEnergyForceBlobs.tif'),'-dtiff')

tableSE_FB=table(SE_FB_GroupCellArray,'RowNames',nameList);
writetable(tableSE_FB,strcat(dataPath,'/strainEnergyForceBlobs.csv'),'WriteRowNames',true)
%% Plotting each - SE-Cell
if isCellSeg
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
end
%% Plotting each - SE-Cell-Peri
if isCellSeg
    SE_CellPeri_GroupCellArray = cellfun(@(x) cell2mat(x),SE_CellPeri_Group,'unif',false);
    samNum=cellfun(@numel,SE_CellPeri_GroupCellArray);
    if length(samNum)>1 && (samNum(1)>10*samNum(2) || samNum(2)>10*samNum(1))
        SE_CellPeri_GroupCellArray = cellfun(@(x) cellfun(@mean,x),SE_CellPeri_Group,'unif',false);
    end
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
end
%% Plotting each - SE-Cell-Inside
if isCellSeg
    SE_CellInside_GroupCellArray = cellfun(@(x) cell2mat(x),SE_CellInside_Group,'unif',false);
    samNum=cellfun(@numel,SE_CellInside_GroupCellArray);
    if length(samNum)>1 && (samNum(1)>10*samNum(2) || samNum(2)>10*samNum(1))
        SE_CellInside_GroupCellArray = cellfun(@(x) cellfun(@mean,x),SE_CellInside_Group,'unif',false);
    end
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
end
%% Plotting each - SE-FOV
SE_FOV_GroupCellArray = cellfun(@(x) cell2mat(x),SE_FOV_Group,'unif',false);
samNum=cellfun(@numel,SE_FOV_GroupCellArray);
if length(samNum)>1 && (samNum(1)>10*samNum(2) || samNum(2)>10*samNum(1))
    SE_FOV_GroupCellArray = cellfun(@(x) cellfun(@mean,x),SE_FOV_Group,'unif',false);
end
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
samNum=cellfun(@numel,SEDenGroupFB_CellArray);
if length(samNum)>1 && (samNum(1)>10*samNum(2) || samNum(2)>10*samNum(1))
    SEDenGroupFB_CellArray = cellfun(@(x) cellfun(@mean,x),SEDen_FB_Group,'unif',false);
end
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
if isCellSeg
    SEDenGroupCell_CellArray = cellfun(@(x) cell2mat(x),SEDen_Cell_Group,'unif',false);
    samNum=cellfun(@numel,SEDenGroupCell_CellArray);
    if length(samNum)>1 && (samNum(1)>10*samNum(2) || samNum(2)>10*samNum(1))
        SEDenGroupCell_CellArray = cellfun(@(x) cellfun(@mean,x),SEDen_Cell_Group,'unif',false);
    end
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
end
%% Strain energy density - CellPeri
if isCellSeg
    SEDenGroupCellPeri_CellArray = cellfun(@(x) cell2mat(x),SEDen_CellPeri_Group,'unif',false);
    samNum=cellfun(@numel,SEDenGroupCellPeri_CellArray);
    if length(samNum)>1 && (samNum(1)>10*samNum(2) || samNum(2)>10*samNum(1))
        SEDenGroupCellPeri_CellArray = cellfun(@(x) cellfun(@mean,x),SEDen_CellPeri_Group,'unif',false);
    end
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
end
%% Strain energy density - CellInside
if isCellSeg
    SEDenGroupCellInside_CellArray = cellfun(@(x) cell2mat(x),SEDen_CellInside_Group,'unif',false);
    samNum=cellfun(@numel,SEDenGroupCellInside_CellArray);
    if length(samNum)>1 && (samNum(1)>10*samNum(2) || samNum(2)>10*samNum(1))
        SEDenGroupCellInside_CellArray = cellfun(@(x) cellfun(@mean,x),SEDen_CellInside_Group,'unif',false);
    end
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
end
%% Strain energy density - FOV
SEDenGroupFOV_CellArray = cellfun(@(x) cell2mat(x),SEDen_FOV_Group,'unif',false);
    samNum=cellfun(@numel,SEDenGroupFOV_CellArray);
    if length(samNum)>1 && (samNum(1)>10*samNum(2) || samNum(2)>10*samNum(1))
        SEDenGroupFOV_CellArray = cellfun(@(x) cellfun(@mean,x),SEDen_FOV_Group,'unif',false);
    end
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
    samNum=cellfun(@numel,totForceFBCellArray);
    if length(samNum)>1 && (samNum(1)>10*samNum(2) || samNum(2)>10*samNum(1))
        totForceFBCellArray = cellfun(@(x) cellfun(@mean,x),totalForce_FB_Group,'unif',false);
    end
h1=figure; 
boxPlotCellArray(totForceFBCellArray,nameList,1,false,true)
ylabel('Total force (nN)')
title('Total force in force blobs')
hgexport(h1,strcat(figPath,'/totForceFB'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/totForceFB'),'-v7.3')
print(h1,strcat(figPath,'/totForceFB.tif'),'-dtiff')

tableTotalForce_FB=table(totForceFBCellArray,'RowNames',nameList);
writetable(tableTotalForce_FB,strcat(dataPath,'/totalForce_ForceBlobs.csv'),'WriteRowNames',true)
%% avg force in all individual force blob in Cell
forceFBCellCellArray = avgForce_FBindiv_Group;%cellfun(@(x) cell2mat(x),avgForce_FBindiv_Group,'unif',false);
h1=figure; 
boxPlotCellArray(forceFBCellCellArray,nameList,1,false,true)
ylabel('Avereage force (nN)')
title('Average force in force blobs in Cells (per force blob)')
hgexport(h1,strcat(figPath,'/avgForceFBinCell'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/avgForceFBinCell'),'-v7.3')
print(h1,strcat(figPath,'/avgForceFBinCell.tif'),'-dtiff')

tableAvgForceAllFBinCell=table(forceFBCellCellArray,'RowNames',nameList);
writetable(tableAvgForceAllFBinCell,strcat(dataPath,'/forceAllFBs.csv'),'WriteRowNames',true)
%% avg force of average force blob in Cell
avgforceFBCellArray = avgForce_FB_Group;%cellfun(@(x) cell2mat(x),avgForce_FB_Group,'unif',false);
h1=figure; 
boxPlotCellArray(avgforceFBCellArray,nameList,1,false,true)
ylabel('Avereage force (nN)')
title('Average force of force blobs in cell per cell')
hgexport(h1,strcat(figPath,'/avgForceFBinCell'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/avgForceFBinCell'),'-v7.3')
print(h1,strcat(figPath,'/avgForceFBinCell.tif'),'-dtiff')

try
    tableAvgForceFBinCell=table(avgforceFBCellArray,'RowNames',nameList);
    writetable(tableAvgForceFBinCell,strcat(dataPath,'/avgForceFBsCell.csv'),'WriteRowNames',true)
catch
    disp('Not all movieList had processed new force blob analysis')
end
%% total force of average force blob in Cell
totforceFBCellArray = totForce_FBInCell_Group;%cellfun(@(x) cell2mat(x),totForce_FBInCell_Group,'unif',false);
h1=figure; 
boxPlotCellArray(totforceFBCellArray,nameList,1,false,true)
ylabel('Total force (nN)')
title('Total force of force blobs in cell (per cell)')
hgexport(h1,strcat(figPath,'/totForceFBinCell'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/totForceFBinCell'),'-v7.3')
print(h1,strcat(figPath,'/totForceFBinCell.tif'),'-dtiff')

try
    tableTotForceFBinCell=table(totforceFBCellArray,'RowNames',nameList);
    writetable(tableTotForceFBinCell,strcat(dataPath,'/totForceFBsCell.csv'),'WriteRowNames',true)
catch
    disp('Not all movieList had processed new force blob analysis')
end
%% Total force - Cell
if isCellSeg
    totForceCell_CellArray = cellfun(@(x) cell2mat(x),totalForce_Cell_Group,'unif',false);
    samNum=cellfun(@numel,totForceCell_CellArray);
    if length(samNum)>1 && (samNum(1)>10*samNum(2) || samNum(2)>10*samNum(1))
        totForceCell_CellArray = cellfun(@(x) cellfun(@mean,x),totalForce_Cell_Group,'unif',false);
    end
    h1=figure; 
    boxPlotCellArray(totForceCell_CellArray,nameList,1,false,true)
    ylabel('Total force (nN)')
    title('Total force in a cell')
    hgexport(h1,strcat(figPath,'/totForceCell'),hgexport('factorystyle'),'Format','eps')
    hgsave(h1,strcat(figPath,'/totForceCell'),'-v7.3')
    print(h1,strcat(figPath,'/totForceCell.tif'),'-dtiff')

    tableTotalForce_Cell=table(totForceCell_CellArray,'RowNames',nameList);
    writetable(tableTotalForce_Cell,strcat(dataPath,'/totalForce_Cell.csv'),'WriteRowNames',true)
end
%% Integrated Displacement - Cell
if isCellSeg
    totDispCell_CellArray = cellfun(@(x) cell2mat(x),totalDisp_Cell_Group,'unif',false);
    samNum=cellfun(@numel,totDispCell_CellArray);
    if length(samNum)>1 && (samNum(1)>10*samNum(2) || samNum(2)>10*samNum(1))
        totDispCell_CellArray = cellfun(@(x) cellfun(@mean,x),totalForce_Cell_Group,'unif',false);
    end
    h1=figure; 
    boxPlotCellArray(totDispCell_CellArray,nameList,1,false,true)
    ylabel('Displacement, integrated (um^3)')
    title('Integrated displacement over cell area')
    hgexport(h1,strcat(figPath,'/totDispCell'),hgexport('factorystyle'),'Format','eps')
    hgsave(h1,strcat(figPath,'/totDispCell'),'-v7.3')
    print(h1,strcat(figPath,'/totDispCell.tif'),'-dtiff')

    tableTotDispCell_CellArray=table(totDispCell_CellArray,'RowNames',nameList);
    writetable(tableTotDispCell_CellArray,strcat(dataPath,'/totDispCell.csv'),'WriteRowNames',true)
end

%% Total force - CellPeri
if isCellSeg
    totForceCellPeri_CellArray = cellfun(@(x) cell2mat(x'),totalForce_CellPeri_Group,'unif',false);
    samNum=cellfun(@numel,totForceCellPeri_CellArray);
    if length(samNum)>1 && (samNum(1)>10*samNum(2) || samNum(2)>10*samNum(1))
        totForceCellPeri_CellArray = cellfun(@(x) cellfun(@mean,x),totalForce_CellPeri_Group,'unif',false);
    end
    h1=figure; 
    boxPlotCellArray(totForceCellPeri_CellArray,nameList,1,false,true)
    ylabel('Total force (nN)')
    title('Total force in a cell periphery')
    hgexport(h1,strcat(figPath,'/totForceCellPeri'),hgexport('factorystyle'),'Format','eps')
    hgsave(h1,strcat(figPath,'/totForceCellPeri'),'-v7.3')
    print(h1,strcat(figPath,'/totForceCellPeri.tif'),'-dtiff')

    tableTotalForce_CellPeri=table(totForceCellPeri_CellArray,'RowNames',nameList);
    writetable(tableTotalForce_CellPeri,strcat(dataPath,'/totalForce_CellPeri.csv'),'WriteRowNames',true)
end
%% Total force - CellInside
if isCellSeg
    totForceCellInside_CellArray = cellfun(@(x) cell2mat(x'),totalForce_CellInside_Group,'unif',false);
    samNum=cellfun(@numel,totForceCellInside_CellArray);
    if length(samNum)>1 && (samNum(1)>10*samNum(2) || samNum(2)>10*samNum(1))
        totForceCellInside_CellArray = cellfun(@(x) cellfun(@mean,x),totalForce_CellInside_Group,'unif',false);
    end
    h1=figure; 
    boxPlotCellArray(totForceCellInside_CellArray,nameList,1,false,true)
    ylabel('Total force (nN)')
    title('Total force in a cell interior')
    hgexport(h1,strcat(figPath,'/totForceCellInside'),hgexport('factorystyle'),'Format','eps')
    hgsave(h1,strcat(figPath,'/totForceCellInside'),'-v7.3')
    print(h1,strcat(figPath,'/totForceCellInside.tif'),'-dtiff')

    tableTotalForce_CellInside=table(totForceCellInside_CellArray,'RowNames',nameList);
    writetable(tableTotalForce_CellInside,strcat(dataPath,'/totalForce_CellInside.csv'),'WriteRowNames',true)
end
%% Total force - FOV
totForceFOV_CellArray = cellfun(@(x) cell2mat(x),totalForce_FOV_Group,'unif',false);
    samNum=cellfun(@numel,totForceFOV_CellArray);
    if length(samNum)>1 && (samNum(1)>10*samNum(2) || samNum(2)>10*samNum(1))
        totForceFOV_CellArray = cellfun(@(x) cellfun(@mean,x),totalForce_FOV_Group,'unif',false);
    end
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
if isCellSeg
    spreadArea_CellArray = cellfun(@(x) cell2mat(x),spreadArea_Group,'unif',false);
    samNum=cellfun(@numel,spreadArea_CellArray);
    if length(samNum)>1 && (samNum(1)>10*samNum(2) || samNum(2)>10*samNum(1))
        spreadArea_CellArray = cellfun(@(x) cellfun(@mean,x),spreadArea_Group,'unif',false);
    end
    h1=figure; 
    boxPlotCellArray(spreadArea_CellArray,nameList,1,false,true)
    ylabel('Spread area (um2)')
    title('Cell spread area')
    hgexport(h1,strcat(figPath,'/spreadArea'),hgexport('factorystyle'),'Format','eps')
    hgsave(h1,strcat(figPath,'/spreadArea'),'-v7.3')
    print(h1,strcat(figPath,'/spreadArea.tif'),'-dtiff')

    tableSpreadArea=table(spreadArea_CellArray,'RowNames',nameList);
    writetable(tableSpreadArea,strcat(dataPath,'/spreadArea.csv'),'WriteRowNames',true)
end
%% save entire workspace for later
save([dataPath filesep 'allData.mat'])
