% 
% sample API calls with leverjs.net. 
%   will read segmentations and generate masks for each image frame in each
%   dataset. remove the :1 and uncomment the length(...) in the for
%   statements to run the whole tamale
%
%       (c) all rights andy cohen. 6/8/2017. acohen@coe.drexel.edu
%       this is code/api still in development. it is not licensed for
%       redistribution (yet). stay tuned...
%
% code written against MATLAB 2017A. jsondecode needs at least 2016b.

% API Summary
% URL_ROOT/LEVER                - list all datasets at URL_ROOT
% URL_ROOT/:dataset/constants   - get image metatdata for :dataset
% URL_ROOT/:dataset/image/t     - get image from time :t for :dataset
% URL_ROOT/:dataset/cells/t     - get segmentations from time :t for :dataset

clc;

URL_ROOT='https://leverjs.net/Danuser_texture/';
% URL_ROOT='https://leverjs.net/Danuser/';
% URL_ROOT='https://leverjs.net/Danuser/May2017/';

% str=webread([URL_ROOT 'LEVER']);
szDatasetNames=webread([URL_ROOT 'LEVER']);
% szDatasetNames=jsondecode(str);

%% Windows
% outLeverDname = 'C:\Assaf\Research\LCH\data\LEVER\masks\';
% metaDataDname = 'C:\Users\assafza\Google Drive\Research\PostDoc\ResearchInProgress\QuantitativeLiveCellHistology\';

%% BioHPC
outLeverDname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/LEVER/masks/';
metaDataDname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/MetaData/';

expsMetaData = [metaDataDname 'ExperimentsAll20171016a.mat'];


load(expsMetaData);

% for iDataset=1:length(szDatasetNames) 
for iDataset=1:length(szDatasetNames) 
        
    try
        CONSTANTS=webread([URL_ROOT szDatasetNames{iDataset} '/constants']);
    catch e
        try
            CONSTANTS=webread([URL_ROOT szDatasetNames{iDataset} '/constants']);
        catch ee
            warning('failed twice connecting to URL %s',[URL_ROOT szDatasetNames{iDataset} '/constants']);
            continue;
        end
    end
    % CONSTANTS.imageData has resolution, metadata, etc.
    % CONSTANTS.ui* is internal
    dname = szDatasetNames{iDataset};
    strs = strsplit(dname,'_');
    
    fprintf(sprintf('\n\n%d\n%s\n',iDataset, dname));
    
    dateStr = strs{1};
    foundExp = length(find(strcmp(dateStr,metaData.experiments.date)));
    if foundExp ~= 1
        warning('%s found %d times in metaData!\n',dateStr,foundExp);
        continue;  
    end
    
    %% Adjust filename (# 0 digits)
    curDname = getNDigits(dateStr,metaData,dname);   
    
    curDnameOut = [outLeverDname filesep curDname filesep];
    
    if ~exist(curDnameOut,'dir')
        mkdir(curDnameOut);
    end
    
    for t=1:CONSTANTS.imageData.NumberOfFrames
        %     for t = 110 : 140
        
        fprintf(sprintf('%d ',t));
        
        curFnameOut = sprintf('%s%d.mat',curDnameOut,t);
        
        if exist(curFnameOut,'file')
            continue;
        end
        
        % get the image for time t
        im=webread([URL_ROOT szDatasetNames{iDataset} '/image/' num2str(t)]);
        %         imagesc(im);colormap(gray);hold on
        
        % get all the cells at time t
        tCells=webread([URL_ROOT szDatasetNames{iDataset} '/cells/' num2str(t)]);
        tCells=jsondecode(tCells);
        % tCells.surface is the outline of the cell
        % tCells.cellID is the segmentation id
        % tCells.trackID is the track ID for that segmentation. 
        
        % nCellMask will be the filled regions inside each cell
        % each region will be numbered by cellID
        ROI_LEVER=false(size(im));
        
        % create masks for each cell
        
        for iCell=1:length(tCells) 
            try
                bw=roipoly(im,tCells(iCell).surface(:,1),tCells(iCell).surface(:,2));
                ROI_LEVER = ROI_LEVER | bw;
            catch e
                load(sprintf('%s%d.mat',curDnameOut,t-1)); % use mask of previous frame
            end            
        end
        save(curFnameOut,'ROI_LEVER');%d
    end
end

%% getNDigits - based on # of locations imaged!
function dnameOut = getNDigits(dateStr,metaData,dname)
ind = find(strcmp(dateStr,metaData.experiments.date));
assert(length(ind) == 1);
n = metaData.experiments.ns{1};
if n > 9 
    nDigits = 2;
end

dnameOut = dname;

if nDigits == 2 && strcmp('s',dname(end-1))
    dnameOut = [dname(1:end-1) '0' dname(end)];
end

end

        