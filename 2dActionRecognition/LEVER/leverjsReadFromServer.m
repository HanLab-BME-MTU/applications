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


URL_ROOT='https://leverjs.net/Danuser/';
str=webread([URL_ROOT 'LEVER']);
szDatasetNames=jsondecode(str);

for iDataset=1:1 % length(szDatasetNames) 
    
    CONSTANTS=webread([URL_ROOT szDatasetNames{iDataset} '/constants']);
    % CONSTANTS.imageData has resolution, metadata, etc.
    % CONSTANTS.ui* is internal
    
    for t=1:CONSTANTS.imageData.NumberOfFrames
        
        % get the image for time t
        im=webread([URL_ROOT szDatasetNames{iDataset} '/image/' num2str(t)]);
        imagesc(im);colormap(gray);hold on
        
        % get all the cells at time t
        tCells=webread([URL_ROOT szDatasetNames{iDataset} '/cells/' num2str(t)]);
        tCells=jsondecode(tCells);
        % tCells.surface is the outline of the cell
        % tCells.cellID is the segmentation id
        % tCells.trackID is the track ID for that segmentation. 
        
        % nCellMask will be the filled regions inside each cell
        % each region will be numbered by cellID
        nCellMask=0*im;
        
        % draw each cell
        
        for iCell=1:length(tCells)
            % jsColor is the javascript (hex, 8 bit) RGB
            
            jsColor=tCells(iCell).trackInfo.color.background;
            jsColor=jsColor(2:end); % drop the leading #
            
            colorTrack=[]; % the matlab 1x3 double [r,g,b]
            colorTrack(1)=hex2dec(jsColor(1:2))/255;
            colorTrack(2)=hex2dec(jsColor(3:4))/255;
            colorTrack(3)=hex2dec(jsColor(5:6))/255;
            
            
            plot(tCells(iCell).surface(:,1),tCells(iCell).surface(:,2),'color',colorTrack);
            bw=roipoly(im,tCells(iCell).surface(:,1),tCells(iCell).surface(:,2));
            nCellMask=nCellMask + tCells(iCell).trackID*uint8(bw);
        end
        drawnow
        % now -- nCellMask is the filled interior pixels of each segmented
        % cell for the current image frame
        % each connected component is numbered by trackID
        
    end
end

        