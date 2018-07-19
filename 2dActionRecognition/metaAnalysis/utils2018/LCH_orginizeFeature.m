%% 
% Obsolete function! Now implemented in LCH_orginizeFeatures2018
function [] = LCH_orginizeFeature(baseDname,dateStr,featsStr)

close all; clc;
always = true;

inFeatsFname = [baseDname filesep dateStr filesep dateStr '_' featsStr '.mat']; 
outFeatsDname = [baseDname filesep dateStr filesep featsStr];
if ~exist(outFeatsDname,'dir')
    mkdir(outFeatsDname);
end

outFeatsFname = [outFeatsDname filesep featsStr '.mat'];

if ~exist(outFeatsFname,'file') || always
    load(inFeatsFname); % allCellsMovieData
    
    nCells = length(allCellsMovieData);
    cellFeats = cell(1,nCells);
    cellTypes = cell(1,nCells);
    metEffs = nan(1,nCells);
    dates = cell(1,nCells);
    sources = cell(1,nCells);
    
    for icell = 1 : nCells
        curCell = allCellsMovieData{icell};
        cellFeats{icell} = curCell.featureMap(featsStr);
        cellTypes{icell} = curCell.cellType;
        metEffs(icell) = curCell.metEff;
        dates{icell} = curCell.date;
        sources{icell} = LCH_getSource(cellTypes{icell});
    end
    
    save(outFeatsFname,'nCells','cellFeats','cellTypes','metEffs','dates','sources');
end
end



	




