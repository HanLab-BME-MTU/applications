function allCellsMovieData = importByFeatureToCellExplorer2018(analysisDname,accFeatsFilename,allCellsDname)

addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition/'));

if nargin < 2        
    analysisDname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/All/';
                    
    featsFname = 'LEVER_LBP_SHAPE_all.mat';
    accFeatsFilename = ['/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/metaAnalysis/LEVER_LBP_SHAPE/' featsFname];
    dateStr = '1-Jan-2018';
    
    allCellsDname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/CellExplorerData/';        
    
end

load(accFeatsFilename); % allCells
nCellTypeDate = length(allCells.strs);

featsStrID = allCells.featsStrID;

% allCellsFname = [allCellsDname filesep date '_' featsStrID]; % should be
% dateStr
allCellsFname = [allCellsDname filesep dateStr '_' featsStrID];

allCellsMovieData = {};
ncells = 0;
for iCellTypeDate = 1 : nCellTypeDate 
    expStr = allCells.expStr{iCellTypeDate};
    expDate = allCells.date{iCellTypeDate};
    cellType = allCells.cellType{iCellTypeDate};
    metEff = allCells.metEff{iCellTypeDate};        
    
    locationData = allCells.locations{iCellTypeDate};
    nLocations = length(locationData.locationTXY);
    
    for ilocation = 1 : nLocations
        curLocationStr = locationData.locationStr{ilocation};
        curLocationFeats = locationData.locationFeats{ilocation};
        curLocationTXY = locationData.locationTXY{ilocation};        
                
        %% 
        nCellsLocation = length(curLocationTXY);
        
        for icell = 1 : nCellsLocation
            ncells = ncells + 1;
            curCellMovieData.expStr = expStr;
            curCellMovieData.date = expDate;
            curCellMovieData.cellType = cellType;
            curCellMovieData.metEff = metEff;
            
            curCellMovieData.ts = curLocationTXY{icell}.ts;
            curCellMovieData.xs = curLocationTXY{icell}.xs;
            curCellMovieData.ys = curLocationTXY{icell}.ys;
            
            curCellMovieData.locationStr = curLocationStr;                        
            
            % Unique cell ID (assumes the trajecotry length is constant)
            curCellMovieData.key = [...
                date '_' cellType ...
                '_s' curCellMovieData.locationStr ...
                '_t' num2str(curCellMovieData.ts(1)) ...
                '_x' num2str(round(curCellMovieData.xs(1))) ...
                '_y' num2str(round(curCellMovieData.ys(1)))];
            
            mdFname = [analysisDname curCellMovieData.expStr filesep...
                curCellMovieData.expStr '_s' curCellMovieData.locationStr filesep...
                curCellMovieData.expStr '_s' curCellMovieData.locationStr '.mat'];                        
            
            % For locations where the location do not have a trailing zero
            if ~exist(mdFname,'file')
                mdFname = [analysisDname curCellMovieData.expStr filesep...
                    curCellMovieData.expStr '_s' sprintf('%d',curCellMovieData.locationStr) filesep...
                    curCellMovieData.expStr '_s' sprintf('%d',curCellMovieData.locationStr) '.mat'];
            end
            
            % Movie data to access the cell movies
            curCellMovieData.MD = mdFname;
            
            % Empty annotation map
            curCellMovieData.annotations = containers.Map(...
                {'blebbing', 'moving', 'short extensions',...
                'dead','debris','out of focus','plastic',...
                'stretched','collagen deformation', 'long extensions',...
                'multiple cells', 'big', 'small', 'medium'},...
                {nan,nan,nan,nan,nan,...
                nan,nan,nan,nan,nan,...
                nan,nan,nan,nan});
            
            % features
            curFeats = curLocationFeats(:,icell);
            curCellMovieData.featureMap = containers.Map({featsStrID},{curFeats});
            
            % include cell
            allCellsMovieData{ncells} = curCellMovieData;
        end
    end    
end

save(allCellsFname,'allCellsMovieData');

end