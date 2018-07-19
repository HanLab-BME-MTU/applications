%% 
% Link cell annotation to features (can be adjusted to any features)
% Input: Andres annotations, features
% Output: Linked input-output
function [] = linkAnnotationToFeature(annotationFname)

close all; clc;
always = true;

if nargin < 2
    annotationFname = 'C:\Assaf\LCH_CellExplorer\dopeCheckPoint_Anno09Jun2017_1037.mat';
    featuresFname = 'C:\Assaf\LCH_CellExplorer\14-May-2017_LBP_dLBP_1.mat';
    annotation2FeaturesFname = 'C:\Assaf\LCH_CellExplorer\annotation2feats_14-May-2017_LBP_dLBP_1.mat';
end

if ~exist(annotation2FeaturesFname,'file') || always
    
    load(annotationFname); % cellAnnotations
    
    cellData = cellAnnotations.cellData;
    keys = {cellData.key};
    annos = {cellData.annotations};
    repeats = cell2mat({cellData.repeatFlag});
    cellexpr = {cellData.cellexpr};
    
    nAnnotations = length(keys);
    
    annotationInfo = {};
    nAnnotCells = 0;
    for icell = 1 : nAnnotations
        if ~repeats(icell)
            nAnnotCells = nAnnotCells + 1;
            annotationInfo{nAnnotCells} = parseCellExplorerKey(keys{icell},cellexpr{icell},annos{icell});
        end
    end
    
    %% link cell annotation --> feature
    load(featuresFname); % allCellsMovieData
    
    cellFeats = cell(1,nAnnotCells);
    cellTypes = cell(1,nAnnotCells);
    metEffs = nan(1,nAnnotCells);
    dates = cell(1,nAnnotCells);    
    sources = cell(1,nAnnotCells);
    cellAnnotations = cell(1,nAnnotCells);
    
    for icell = 1 : nAnnotCells
        curCellAnnot = annotationInfo{icell};
        [cellFeats{icell},cellTypes{icell},metEffs(icell),dates{icell},sources{icell},cellAnnotations{icell}] = ...
            linkCellAnnot2Feats(curCellAnnot,allCellsMovieData);         
        cellTypes{icell} = curCellAnnot.cellType;        
    end
    
    save(annotation2FeaturesFname,'cellFeats','cellTypes','metEffs','dates','nAnnotCells','annotationInfo','sources','cellAnnotations');
else
    load(annotation2FeaturesFname);
end
%%

end

%%
function [feats,cellType,metEff,dateStr,sourceStr,cellAnnot] = linkCellAnnot2Feats(curCellAnnot,allCellsMovieData)
ncells = length(allCellsMovieData);

x = curCellAnnot.x;
y = curCellAnnot.y;
exp = curCellAnnot.exp;
location = curCellAnnot.location(2:end);
stime = curCellAnnot.stime;

for icell = 1 : ncells
    curCell = allCellsMovieData{icell};
    if ~strcmp(exp,curCell.expStr) || ~strcmp(location,curCell.locationStr)
        continue;
    end
    if stime ~= curCell.ts(1) % Caveat: changing tracking / feature extraction could change times..
        continue;
    end
    if abs(x - curCell.xs(1)) > 10 || abs(y - curCell.ys(1)) > 10
        continue;
    end
    feats = curCell.featureMap('LBP_dLBP_1');
    cellType = curCell.cellType;
    metEff = curCell.metEff;
    dateStr = curCell.date;
    sourceStr = getSource(cellType);
    cellAnnot = curCellAnnot.annotation;
    return;
end

warning('%s_s%s_t%d_x%d_y%d annotation->features failed',exp,location,stime,x,y);
feats = nan;
cellType = nan;
metEff = nan;
dateStr = nan;
sourceStr = nan;

end

%% 
function sourceStr = getSource(cellType)
if sum(strcmpi(cellType,{'m481','m214','m530','m610','um12','m514','m405','m528','m498','ut8','m634','m597'}))
    sourceStr = 'Tumors';
else
    if sum(strcmpi(cellType,{'atcc','m116'}))
        sourceStr = 'Melanocytes';
    else
        if sum(strcmpi(cellType,{'a375','mv3','wm3670','wm1361','wm1366','skmel2'}))
            sourceStr = 'CellLines';
        else
            warning('Failed finding source %s', cellType);
        end
    end
end
end


	




