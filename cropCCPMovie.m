% Function for cropping movies stored as individual TIFF files;
% uses the data structure returned by loadConditionData()

% Francois Aguet, 06/27/2011

function cropCCPMovie(data)

nCh = length(data.channels);

pDirs = cellfun(@(x) getParentDir(x), data.channels, 'UniformOutput', false);
if ~all(strcmp(pDirs{1}, pDirs)) || nCh==1
    error('This function currently only works for multi-channel data with channel folders.');
end

plotFrame(data, [], 1, 1:min(nCh,3), 'Mode', 'RGB');
set(gcf, 'NumberTitle', 'off', 'Name', 'Crop: first frame');
hr = imrect();
pos = round(wait(hr));
close(gcf);

plotFrame(data, [], data.movieLength, 1:min(nCh,3), 'Mode', 'RGB');
set(gcf, 'NumberTitle', 'off', 'Name', 'Crop: last frame');
hr = imrect(gca, pos);
pos = round(wait(hr));
close(gcf);
drawnow;

fprintf('Movie crop: writing files...');
for c = 1:nCh
    chDir = data.channels{c};
    chDir = [chDir(1:end-1) '_crop' filesep];
    [~,~] = mkdir(chDir);
    parfor f = 1:data.movieLength
        frame = imread(data.framePaths{c}{f});
        frame = frame(pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3));
        [~, fname, fext] = fileparts(data.framePaths{c}{f});
        imwrite(frame, [chDir fname '_crop' fext]);
    end
end
fprintf(' done.\n');