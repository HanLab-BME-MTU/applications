function displayFrame(h, iFrame)

settings = get(h, 'UserData');

backgroundIndex = settings.backgroundIndex;
numImages = settings.numImages;
firstNonZeroChannel = settings.firstNonZeroChannel;

% Set the title of the window
set(h, 'Name', ['fsmDataViewer: frame (' num2str(iFrame) '/' num2str(numImages) ')' ]);

% Read channels

cmap = 'default';

C = [];
if firstNonZeroChannel
    switch backgroundIndex
        case 1, % from raw data
            path = settings.channels{firstNonZeroChannel}.path;
            fileName = settings.channels{firstNonZeroChannel}.fileList(iFrame).name;
            I = double(imread([path filesep fileName]));
            C = cat(3, zeros(size(I, 1), size(I, 2), firstNonZeroChannel - 1),...
                I / max(I(:)),...
                zeros(size(I, 1), size(I, 2), 3 - firstNonZeroChannel));
            for iChannel = firstNonZeroChannel+1:3
                path = settings.channels{iChannel}.path;
                if ~isempty(path)
                    fileName = settings.channels{iChannel}.fileList(iFrame).name;
                    I = double(imread([path filesep fileName]));
                    C(:, :, iChannel) = I / max(I(:));
                end
            end
            
        case 2, % from speed map
            path = settings.channels{firstNonZeroChannel}.path;
            fileName = settings.channels{firstNonZeroChannel}.fileList(iFrame).name;
            load([path, filesep fileName]);
            if ~exist('speedMap', 'var')
                error(['Unable to load the ''speedMap'' variable from disk at frame ' iFrame]);
            end
            C = speedMap;
            cmap = 'jet';
       
        case 3, % from kinetics map
    end
end

% Add the mask

M = [];
if numel(settings.mask.fileList)
    path = settings.mask.path;
    fileName = settings.mask.fileList(iFrame).name;
    M = imread([path filesep fileName]);
end

if ~isempty(C)
    if ~isempty(M)
        for iChannel = 1:size(C, 3)
            C(:, :, iChannel) = C(:, :, iChannel) .* M(:, :);
        end
    end
else
    C = M;
end

% Display channels

imshow(C, []);
colormap(cmap);

end