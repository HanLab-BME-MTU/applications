function dataLayer = dispatchFilesToFrames(fileList, nFrames)

if numel(fileList) ~= nFrames
    error('Number of layer files does not match number of frames.');
end

dataLayer = cell(nFrames,1);

for iFrame = 1:nFrames
    tmp = load(fileList{iFrame});
    names = fieldnames(tmp);
    if numel(names) ~= 1
        error('The number of variables in the file should be 1.');
    end
    
    dataLayer{iFrame} = tmp.(names{1});
end
    