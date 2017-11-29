function allCellsMovieData = importByFeatureToCellExplorer(analysisDname,accFeatsFilename,allCellsDname)

addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition/'));

if nargin < 2
    %     analysisDname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/ThirdGenData/';
    %     accFeatsFilename = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/metaAnalysis/LBP_dLBP_120_ThirdGen201701/LBP_dLBP_120_ThirdGen201701_1_allInfo.mat';
    
    analysisDname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/All/';
    %     featsFname = 'LBP_dLBP_120_All201705_1_allInfo.mat';
    %     featsFname = 'LBP_dLBP_120_All201705_2_allInfo.mat';
    %     featsFname = 'LBP_dLBP_120_All201705_3_allInfo.mat';
    %     featsFname = 'LBP_dLBP_120_All201705_4_allInfo.mat';
    %     accFeatsFilename = ['/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/metaAnalysis/LBP_dLBP_120_All201705/' featsFname];
    
    %     dateStr = '14-May-2017';
    
    %     featsFname = 'LBP_dLBP_120_Clones201708_1_allInfo.mat';
    %     featsFname = 'LBP_dLBP_120_Clones201708_2_allInfo.mat';
    %     featsFname = 'LBP_dLBP_120_Clones201708_3_allInfo.mat';
    %     featsFname = 'LBP_dLBP_120_Clones201708_4_allInfo.mat';
    %     accFeatsFilename = ['/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/metaAnalysis/LBP_dLBP_120_Clones201708/' featsFname];
    %     dateStr = '10-Aug-2017';
    
    %     featsFname = 'LBP_dLBP_120_All201708_1_allInfo.mat';
    %     featsFname = 'LBP_dLBP_120_All201708_2_allInfo.mat';
    %     featsFname = 'LBP_dLBP_120_All201708_3_allInfo.mat';
    %     featsFname = 'LBP_dLBP_120_All201708_4_allInfo.mat';
    %     accFeatsFilename = ['/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/metaAnalysis/LBP_dLBP_120_All201708/' featsFname];
    %     dateStr = '31-Aug-2017';
               
    featsFname = 'LBP_dLBP_60_All201711_1_allInfo.mat';
    accFeatsFilename = ['/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/metaAnalysis/LBP_dLBP_60_All201711/' featsFname];
    dateStr = '12-Nov-2017';
    
    allCellsDname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/CellExplorerData/';        
    
end

load(accFeatsFilename); % allInfo
nCellTypeDate = length(allInfo.strs);

featsStrID = allInfo.featsStrID;

% allCellsFname = [allCellsDname filesep date '_' featsStrID]; % should be
% dateStr
allCellsFname = [allCellsDname filesep dateStr '_' featsStrID];

allCellsMovieData = {};
ncells = 0;
for iCellTypeDate = 1 : nCellTypeDate 
    expStr = allInfo.expStr{iCellTypeDate};
    expDate = allInfo.date{iCellTypeDate};
    cellType = allInfo.cellType{iCellTypeDate};
    metEff = allInfo.metEff{iCellTypeDate};        
    
    locationData = allInfo.fov.locations{iCellTypeDate};
    nLocations = length(locationData.locationTXY);
    
    for ilocation = 1 : nLocations
        curLocationStr = locationData.locationStr{ilocation};
        curLocationFeats = locationData.locationFeats{ilocation};
        curLocationTXY = locationData.locationTXY{ilocation};        
                
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