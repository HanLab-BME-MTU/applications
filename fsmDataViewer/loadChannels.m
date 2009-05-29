function I = loadChannels(settings, iFrame)

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

% Get the channel plugins list
channelPlugins = getPlugins();

if ~numLayerFiles
    if numMaskFiles
        BW = imread([settings.maskPath filesep...
            settings.maskFileNames{iFrame}]);
    end
    if numChannelFiles        
        I = [];
        colors = zeros(3, 1);
        for iChannel = 1:numChannels
            channelTypeID = settings.channels{iChannel}.type;
            channelColor = settings.channels{iChannel}.color;
            
            fileName = [settings.channels{iChannel}.path ...
                filesep settings.channels{1}.fileNames{iFrame}];
            
            J = channelPlugins(channelTypeID).load(fileName);
            I = cat(3, I, J);
            
            colors(iChannel) = channelColor;
        end
       
        [m n c] = size(I);

        if c > 1
            % RGB images are 3 x [0..1] so each channel has to be
            % normalized.
            I = single(I);
            for iChannel = 1:c
                J = I(:, :, iChannel);
                J = J - min(J(:));
                I(:, :, iChannel) = J / max(J(:));
            end
            
            if c == 2
                % Add a blank channel.
                I = cat(3, I, zeros(m, n, class(I)));
                % Set the color of the blank channel to the free color (i.e.
                % different form colors(1) and colors(2)).
                colors(3) = 6 - (colors(1) + colors(2));
            end
            
            % Make channel permutation according to channel colors
            I = I(:, :, colors);            
        end
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
        settings.layers{settings.iNumLayerFiles}.fileNames{iFrame});
    
    no = str2double(no);
    
    if numMaskFiles
        [fileName, found] = findNumberedFileInList(...
            settings.maskFileNames, no);
        assert(found); % Test done in getSettings.m. Must be valid.
        BW = imread([settings.maskPath filesep file]);
    end
    if numChannelFiles
        % TODO: this part is pretty much a copy-past of lines 37:73. Make
        % something to avoid that.
        I = [];
        colors = zeros(3, 1);
        for iChannel = 1:numChannels
            channelTypeID = settings.channels{iChannel}.type;
            channelColor = settings.channels{iChannel}.color;

            fileName = findNumberedFileInList(...
                settings.channels{iChannel}.fileNames, no);
            
            J = channelPlugins(channelTypeID).load([settings.channels{iChannel}.path ...
                filesep fileName]);
            
            I = cat(3, I, J);
            
            colors(iChannel) = channelColor;
        end
       
        [m n c] = size(I);

        if c > 1
            % RGB images are 3 x [0..1] so each channel has to be
            % normalized.
            I = double(I);
            for iChannel = 1:c
                J = I(:, :, iChannel);
                J = J - min(J(:));
                I(:, :, iChannel) = J / max(J(:));
            end
            
            if c == 2
                % Add a blank channel.
                I = cat(3, I, zeros(m, n, class(I)));
                % Set the color of the blank channel to the free color (i.e.
                % different form colors(1) and colors(2)).
                colors(3) = 6 - (colors(1) + colors(2));
            end
            
            % Make channel permutation according to channel colors
            I = I(:, :, colors);            
        end
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