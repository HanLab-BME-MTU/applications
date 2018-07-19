function [imageData, metadata] = uiGetBioImageData(dataFilePath, varargin)

    imageData = [];
    metadata = [];

    p = inputParser;
    p.CaseSensitive = false;
    p.addRequired('dataFilePath', @ischar);
    p.addParamValue('seriesId', [], @(x) (isnumeric(x) && isscalar(x)));
    p.addParamValue('channelDescriptionsAndNamesList', [], @iscell);
    p.parse(dataFilePath, varargin{:});
    
    PARAMETERS = p.Results;
    
    channelDescriptionsAndNamesList = PARAMETERS.channelDescriptionsAndNamesList;

    % initialize bio image reader
    hStatusDialog = waitbar(0, 'Loading Image Data');
    try 
        br = BioImageReader(dataFilePath);    
    catch err
        fprintf( 'ERROR: could not load image data data from file %s', dataFilePath );
        err
        errordlg(err.message, 'Error loading image data from file');
        closeStatusDialog(hStatusDialog);        
        return;
    end

    % ask user to select a series
    if isempty(PARAMETERS.seriesId)
        if br.getSeriesCount() > 1
            [usel, bsel] = settingsdlg('Description', 'Select series', ...
                                       'seriesId', num2cell(1:br.getSeriesCount()));
            if ~strcmpi(bsel, 'ok')
                closeStatusDialog(hStatusDialog);   
                return;
            end
            seriesId = usel.seriesId;
        else
            seriesId = 1;
        end
    else
        seriesId = PARAMETERS.seriesId;
    end

    % get series metadata
    metadata = br.getMetadata(seriesId);

    % make sure that the dataset is 3D
    if numel(metadata.imageSize) ~= 3 || metadata.imageSize(3) < 2
        errordlg( 'Dataset is not 3D' );
        closeStatusDialog(hStatusDialog);        
        return;
    end

    % get image data based on user input
    hDisp = imseriesshow_multichannel( br.getImageData(seriesId, 'timepointIds', 1), ...
                                       'spacing', metadata.pixelSize );  
    set(hDisp, 'Name', 'Image Data Preview');
    title('Image Data at Timepoint 1');

    chArgs = {'separator', 'Set pixel size', ...
              {'Pixel Size in X', 'pixelSize_X'}, metadata.pixelSize(1), ...
              {'Pixel Size in Y', 'pixelSize_Y'}, metadata.pixelSize(2), ...
              {'Pixel Size in Z', 'pixelSize_Z'}, metadata.pixelSize(3), ...
              'separator', 'Select stack of interest', ...
              {'Timepoint', 'timepointId'}, num2cell(1:metadata.numTimePoints)};
          
    if ~isempty(channelDescriptionsAndNamesList)
        for i = 1:numel(channelDescriptionsAndNamesList)
            chArgs{end+1} = channelDescriptionsAndNamesList{i};
            chArgs{end+1} = num2cell(1:metadata.numChannels);
        end
    end

    [usel, bsel] = settingsdlg('Title', 'Set Image properties', ...
                               'Description', sprintf('Image Format - %s\nImage Size - [%s]\nNumber of Channels - %d\nChannel Names - [ %s ]\nNumber of Timepoints - %d\nPixel Type - %s', ...
                                                       metadata.format, ... 
                                                       sprintf(' %d ', metadata.imageSize), ...
                                                       metadata.numChannels, sprintf(' %s ', metadata.channelNames{:}), ...
                                                       metadata.numTimePoints, metadata.pixelType), ...
                               chArgs{:});

    if ~strcmpi(bsel, 'ok')
        close( hDisp );
        closeStatusDialog(hStatusDialog);        
        return;
    end

    
    metadata.timepointId = usel.timepointId;
    
    if ~isempty(channelDescriptionsAndNamesList)
        for i = 1:numel(channelDescriptionsAndNamesList)
            channelIdName = channelDescriptionsAndNamesList{i}{2};
            metadata.(channelIdName) = usel.(channelIdName);
        end
    end
    
    metadata.pixelSize = [usel.pixelSize_X, usel.pixelSize_Y, usel.pixelSize_Z];

    imageData = br.getImageData(metadata.seriesId, ...
                                'timepointIds', metadata.timepointId);

    close( hDisp );
    closeStatusDialog(hStatusDialog);        

end