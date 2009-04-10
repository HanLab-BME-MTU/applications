function [B cmap] = getBackground(settings, iFrame)

backgroundIndex = settings.backgroundIndex;
%numFrames = settings.numFrames;
firstNonEmptyChannel = settings.firstNonEmptyChannel;
numBackgroundFiles = settings.numBackgroundFiles;
%numLayerFiles = settings.numLayerFiles;
numMaskFiles = settings.numMaskFiles;

% Read channels
cmap = 'default';

I = [];
if firstNonEmptyChannel
    switch backgroundIndex
        case 1, % from raw data
            path = settings.channels{firstNonEmptyChannel}.path;
            fileName = settings.channels{firstNonEmptyChannel}.fileNames{iFrame};
            I = double(imread([path filesep fileName]));
            B = cat(3, zeros(size(I, 1), size(I, 2), firstNonEmptyChannel - 1),...
                I / max(I(:)),...
                zeros(size(I, 1), size(I, 2), 3 - firstNonEmptyChannel));
            for iChannel = firstNonEmptyChannel+1:3
                path = settings.channels{iChannel}.path;
                if ~isempty(path)
                    fileName = settings.channels{iChannel}.fileNames{iFrame};
                    I = double(imread([path filesep fileName]));
                    B(:, :, iChannel) = I / max(I(:));
                end
            end

        case 2, % from speed map
            path = settings.channels{firstNonEmptyChannel}.path;
            fileName = settings.channels{firstNonEmptyChannel}.fileNames{iFrame};
            load([path, filesep fileName]);
            if ~exist('speedMap', 'var')
                error(['Unable to load the ''speedMap'' variable for ' path filesep fileName '.']);
            end
            B = speedMap;
            cmap = 'jet';

        case 3, % from kinetics map (REDO IT)
            path = settings.channels{firstNonEmptyChannel}.path;
            fileName = settings.channels{firstNonEmptyChannel}.fileNames{iFrame};
            load([path, filesep fileName]);
            if exist('polyMap', 'var')
                I = polyMap;
                I = I / max(I(:));
                clear polyMap;
            else
                if exist('depolyMap', 'var')
                    I = depolyMap;
                    I = I / min(I(:));
                    clear depolyMap;
                else
                    error(['Unable to load the poly/depoly map variable for ' path filesep fileName '.']);
                end
            end
            B = cat(3, zeros(size(I, 1), size(I, 2), firstNonEmptyChannel - 1),...
                I,...
                zeros(size(I, 1), size(I, 2), 3 - firstNonEmptyChannel));
            for iChannel = firstNonEmptyChannel+1:3
                path = settings.channels{iChannel}.path;
                if ~isempty(path)
                    fileName = settings.channels{iChannel}.fileNames{iFrame};
                    load([path, filesep fileName]);
                    if exist('polyMap', 'var')
                        I = polyMap;
                        I = I / max(I(:));
                        clear polyMap;
                    else
                        if exist('depolyMap', 'var')
                            I = depolyMap;
                            I = I / min(I(:));
                            clear depolyMap;
                        else
                            error(['Unable to load the poly/depoly map variable for ' path filesep fileName '.']);
                        end
                    end
                    B(:, :, iChannel) = I;
                end
            end
    end
end

% Add the mask
if numMaskFiles
    path = settings.mask.path;

    if numBackgroundFiles
        [dummy1, dummy2, no] = getFilenameBody([settings.channels{firstNonEmptyChannel}.path ...
            filesep ...
            settings.channels{firstNonEmptyChannel}.fileNames{iFrame}]);

        fileName = findNumberedFileInList(settings.mask.fileNames, str2double(no));

        M = imread([path filesep fileName]);

        for iChannel = 1:size(B, 3)
            B(:, :, iChannel) = B(:, :, iChannel) .* M;
        end
    else
        fileName = settings.mask.fileNames{iFrame};

        M = imread([path filesep fileName]);

        B = M;
    end
end
end