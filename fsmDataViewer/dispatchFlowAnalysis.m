function dataLayer = dispatchFlowAnalysis(fileList, nFrames)

[~,body] = getFilenameBody(fileList{1});

nAvgFrames = regexp(body,'frames=[0-9]+', 'match');
nAvgFrames = str2double(nAvgFrames{1}(8:end));
hside = floor(nAvgFrames / 2);

if numel(fileList) + 2 * hside ~= nFrames-1
    error('Not enough layer files.');
end

dataLayer = cell(nFrames,1);

for iFile = 1:numel(fileList)
    
    load(fileList{iFile});

    if ~exist('F', 'var')
        error('Unable to find field variable.');
    end

    dataLayer{iFile + hside} = F;
end
