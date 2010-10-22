function dataLayer = dispatchMatrix3ToFrames(fileList, nFrames)

if numel(fileList) ~= 1
    error('Only 1 file is expected');
end

tmp = load(fileList{1});

names = fieldnames(tmp);

if numel(names) ~= 1
    error('The number of variables in the file should be 1.');
end

M = tmp.(names{1});

[nx,ny,nz] = size(M);

if nz ~= nFrames
    error('The size of the matrix does not match the number of frames.');
end

dataLayer = cell(nFrames,1);

for iFrame = 1:nFrames
    dataLayer{iFrame} = M(:,:,iFrame);
end