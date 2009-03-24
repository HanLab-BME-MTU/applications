function displayFrame(h, iFrame)

settings = get(h, 'UserData');

backgroundIndex = settings.backgroundIndex;
numFrames = settings.numFrames;
firstNonEmptyChannel = settings.firstNonEmptyChannel;
numBackgroundFiles = settings.numBackgroundFiles;
%numLayerFiles = settings.numLayerFiles;
numMaskFiles = settings.numMaskFiles;

% Set the title of the window
set(h, 'Name', ['fsmDataViewer: frame (' num2str(iFrame) '/' num2str(numFrames) ')' ]);

% Read channels

cmap = 'default';

C = [];
if firstNonEmptyChannel
    switch backgroundIndex
        case 1, % from raw data
            path = settings.channels{firstNonEmptyChannel}.path;
            fileName = settings.channels{firstNonEmptyChannel}.fileNames{iFrame};
            I = double(imread([path filesep fileName]));            
            C = cat(3, zeros(size(I, 1), size(I, 2), firstNonEmptyChannel - 1),...
                I / max(I(:)),...
                zeros(size(I, 1), size(I, 2), 3 - firstNonEmptyChannel));
            for iChannel = firstNonEmptyChannel+1:3
                path = settings.channels{iChannel}.path;
                if ~isempty(path)
                    fileName = settings.channels{iChannel}.fileNames{iFrame};
                    I = double(imread([path filesep fileName]));
                    C(:, :, iChannel) = I / max(I(:));
                end
            end
            
        case 2, % from speed map
            path = settings.channels{firstNonEmptyChannel}.path;
            fileName = settings.channels{firstNonEmptyChannel}.fileNames{iFrame};
            load([path, filesep fileName]);
            if ~exist('speedMap', 'var')
                error(['Unable to load the ''speedMap'' variable for ' path filesep fileName '.']);
            end
            C = speedMap;
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
            C = cat(3, zeros(size(I, 1), size(I, 2), firstNonEmptyChannel - 1),...
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
                    C(:, :, iChannel) = I;
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
        
        for iChannel = 1:size(C, 3)
            C(:, :, iChannel) = C(:, :, iChannel) .* M;
        end        
    else
        fileName = settings.mask.fileNames{iFrame};
        
        M = imread([path filesep fileName]);
        
        C = M;
    end
end

% Add layers (TODO)

% Display channels
imshow(C, []);
colormap(cmap);

end