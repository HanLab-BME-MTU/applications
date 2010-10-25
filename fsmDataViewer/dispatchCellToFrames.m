function dataLayer = dispatchCellToFrames(fileList, nFrames)

if numel(fileList) ~= 1
    error('Only 1 file is expected');
end

tmp = load(fileList{1});

names = fieldnames(tmp);

if numel(names) ~= 1
    error('The number of variables in the file should be 1.');
end

C = tmp.(names{1});

if numel(C) ~= nFrames
    error('The cell size does not match the number of frames.');
end

dataLayer = cell(nFrames,1);

for iFrame = 1:nFrames
    dataLayer{iFrame} = C{iFrame};
end