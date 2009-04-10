function [B cmap] = getBackground(settings, iFrame)

backgroundIndex = settings.backgroundIndex;
firstNonEmptyChannel = settings.firstNonEmptyChannel;
numMaskFiles = settings.numMaskFiles;

% Set default colormap to gray.
cmap = 'gray';

% Get the number of non-empty channels.
numNonEmptyChannels = 0;
for iChannel = 1:3
    path = settings.channels{iChannel}.path;
    if ~isempty(path)
        numNonEmptyChannels = numNonEmptyChannels + 1;
    end
end

if numNonEmptyChannels
    switch backgroundIndex
        case 1, % From raw data
            
            % Open the image corresponding to the first channel.
            path = settings.channels{firstNonEmptyChannel}.path;
            fileName = settings.channels{firstNonEmptyChannel}.fileNames{iFrame};
            I = double(imread([path filesep fileName]));
            [m n] = size(I);
            
            if numNonEmptyChannels ~= 1
                % Create an RGB image
                B = cat(3, zeros(m, n, firstNonEmptyChannel - 1),...
                    I / max(I(:)),...
                    zeros(m, n, 3 - firstNonEmptyChannel));
                
                for iChannel = firstNonEmptyChannel+1:3
                    path = settings.channels{iChannel}.path;
                    if ~isempty(path)
                        fileName = settings.channels{iChannel}.fileNames{iFrame};
                        I = double(imread([path filesep fileName]));
                        B(:, :, iChannel) = I / max(I(:));
                    end
                end
            else
                % Otherwise, a gray image
                B = I;
            end

        case 2, % From speed map
            
            % Open the speedmap file.
            path = settings.channels{firstNonEmptyChannel}.path;
            fileName = settings.channels{firstNonEmptyChannel}.fileNames{iFrame};
            load([path, filesep fileName]);
            
            % Check the validity of variables.
            if ~exist('speedMap', 'var')
                error(['Unable to load the ''speedMap'' variable for '...
                    path filesep fileName '.']);
            end
            
            B = speedMap;
            
            % Set the proper colormap.
            cmap = 'jet';
            
            % Clear loaded variables.
            clear speedMap;

        case 3, % from kinetics map
            
            % Open the speedmap file for the first channel.
            path = settings.channels{firstNonEmptyChannel}.path;
            fileName = settings.channels{firstNonEmptyChannel}.fileNames{iFrame};
            load([path, filesep fileName]);
            
            % Check the validity of variables.
            if exist('polyMap', 'var')
                I = polyMap;                
                I = I / max(I(:));
                clear polyMap;
            elseif exist('depolyMap', 'var')
                I = depolyMap;
                I = I / min(I(:));
                clear depolyMap;
            else
                error(['Unable to load the poly/depoly map variable for '...
                    path filesep fileName '.']);
            end

            [m n] = size(I);
            
            if numNonEmptyChannels ~= 1
                % Create an RGB image
                B = cat(3, zeros(m, n, firstNonEmptyChannel - 1),...
                    I,...
                    zeros(m, n, 3 - firstNonEmptyChannel));
                
                for iChannel = firstNonEmptyChannel+1:3
                    path = settings.channels{iChannel}.path;
                    if ~isempty(path)
                        fileName = settings.channels{iChannel}.fileNames{iFrame};
                        load([path, filesep fileName]);
                        
                        % Check the validity of variables.
                        if exist('polyMap', 'var')
                            I = polyMap;
                            I = I / max(I(:));
                            clear polyMap;
                        elseif exist('depolyMap', 'var')
                            I = depolyMap;
                            I = I / min(I(:));
                            clear depolyMap;
                        else
                            error(['Unable to load the poly/depoly map variable for '...
                                path filesep fileName '.']);
                        end
                        B(:, :, iChannel) = I;
                    end
                end
            else
                % Otherwise, a gray image
                B = I;
            end
    end
end

% Add the masks
if numMaskFiles
    maskPath = settings.mask.path;

    if numNonEmptyChannels
        path = settings.channels{firstNonEmptyChannel}.path;
        fileName = settings.channels{firstNonEmptyChannel}.fileNames{iFrame};
        [dummy1, dummy2, no] = getFilenameBody([path filesep fileName]);

        fileName = findNumberedFileInList(settings.mask.fileNames, str2double(no));

        M = imread([maskPath filesep fileName]);

        for iChannel = 1:size(B, 3)
            B(:, :, iChannel) = B(:, :, iChannel) .* M;
        end
    else
        fileName = settings.mask.fileNames{iFrame};

        M = imread([maskPath filesep fileName]);

        B = M;
    end
end

end

