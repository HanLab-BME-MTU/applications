function dataLayer = dispatchThetaMap(fileList, nFrames)

if numel(fileList) ~= nFrames
    error('Number of layer files does not match number of frames.');
end

dataLayer = cell(nFrames,1);

for iFrame = 1:nFrames
    load(fileList{iFrame});
    if ~exist('thetaMap', 'var')
        error('Unable to find a variable thetaMap in the file.');
    end

    [ny nx] = size(thetaMap);
    [X Y] = meshgrid(1:nx,1:ny);
    dU = cos(thetaMap);
    dV = sin(thetaMap);
    
    dataLayer{iFrame} = horzcat(Y(:), X(:), dV(:), dU(:));
end
    