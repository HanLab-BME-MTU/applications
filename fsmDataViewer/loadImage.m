function I = loadImage(settings, iFrame)

% Loading the sequence depends on what other kind of data need to be
% displayed. There are 6 cases:
%
% 1. mask only: load the sequence according to settings.maskFileList
%
% 2. channels only: load the sequence according to
%    settings.channels{*}.fileNames
%
% 3. mask + channels: load the sequence according to
%    settings.channels{*}.fileNames 
%
% 4. mask + layers: load the sequence according to
%    settings.layers{settings.inumLayerFiles}.fileNames
%
% 5. layers + channels: load the sequence according to
%    settings.layers{settings.inumLayerFiles}.fileNames
%
% 6. mask + layers + channels: load the sequence according to
%    settings.layers{settings.inumLayerFiles}.fileNames

numChannels = settings.numChannels;
numChannelFiles = settings.numChannelFiles;
numMaskFiles = settings.numMaskFiles;
numLayerFiles = settings.numLayerFiles;

channelTypeNames = {'Raw Images'}; 
channelLoaders = {@imread};

% Get sequence dimension from the first image.

if ~numLayerFiles
    if numMaskFiles
        BW = imread([settings.maskPath filesep...
            settings.maskFileNames{iFrame}]);
    end
    if numChannelFiles        
        I = [];
        for iChannel = 1:numChannels
            channelTypeName = settings.channels{iChannel}.type;
            channelType = strmatch(channelTypeName, channelTypeNames);
            J = channelLoaders{channelType}([settings.channels{iChannel}.path...
                filesep settings.channels{1}.fileNames{iFrame}]);
            I = cat(3, I, J);
        end
        % Fill with empty channel if numChannels == 2
        if numChannels == 2
            [m n] = size(I);
            I = cat(3, I, zeros(m, n, class(I)));
        end
        
        % Make channel permutation according to channel colors
        % TODO
    end
    
    if numMaskFiles % else (case 2)
        if numChannelFiles
            I(BW == 0, :) = 0; % (case 3)
        else
            I = uint8(BW); % (case 1)
        end
        clear BW;
    end
else
    % Get the number of the first layer file name.
    [dummy, body, no] = getFilenameBody(...
        settings.channels{settings.iNumLayerFiles}.fileNames{iFrame});
    
    no = str2double(no);
    
    if numMaskFiles
        [file, found] = findNumberedFileInList(settings.maskFileNames, no);
        assert(found); % Test done in getSettings.m. Must be valid.
        BW = imread([settings.maskPath filesep file]);
    end
    if numChannelFiles
        [file, found] = findNumberedFileInList(settings.channels{1}.fileNames, no);
        assert(found); % Test done in getSettings.m. Must be valid.
        I = imread([settings.maskPath filesep file]);
    end
    
    if numMaskFiles % else (case 5)
        if numChannelFiles
            I(BW == 0) = 0; % (case 6)
        else
            I = uint8(BW); % (case 4)
        end
        clear BW;
    end
end

end