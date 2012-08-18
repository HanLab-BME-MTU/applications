% only supports single-channel data

function subsampleData(data, destPath)

nd = numel(data);

% create condition dir
[~,~] = mkdir([destPath getDirFromPath(getExpDir(data))]);

for i = 1:nd
    cpath = [destPath getDirFromPath(getExpDir(data)) filesep getCellDir(data(i)) '_sub' filesep];
    [~,~] = mkdir(cpath);
    
    fpaths = data(i).framePaths{1};
    fpaths = fpaths(1:2:end);
    
    for f = 1:numel(fpaths)
        copyfile(fpaths{f}, [cpath getCellDir(data(i)) '_' num2str(f, '%.3d') '.tif']);
    end
end

