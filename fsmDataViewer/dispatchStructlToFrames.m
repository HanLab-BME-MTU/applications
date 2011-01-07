function dataLayer = dispatchStructlToFrames(fileList, nFrames)

if numel(fileList) ~= 1
    error('Only 1 file is expected');
end

tmp = load(fileList{1});

names = fieldnames(tmp);

if numel(names) ~= 1
    error('The number of variables in the file should be 1.');
end

S = tmp.(names{1});

if numel(S) ~= nFrames
    error('The struct array size does not match the number of frames.');
end

dataLayer = cell(nFrames,1);

for iFrame = 1:nFrames
    dataLayer{iFrame} = S(iFrame);
end
