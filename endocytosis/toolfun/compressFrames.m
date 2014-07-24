%compressFrames(data) re-writes the frames of the input data sets as LWZ-compressed TIFFs

% Francois Aguet, 07/03/2013

function compressFrames(data)

for i = 1:numel(data)
    for c = 1:numel(data(i).channels)
        parfor f = 1:data(i).movieLength
            frame = imread(data(i).framePaths{c}{f}); %#ok<PFBNS>
            imwrite(frame, data(i).framePaths{c}{f}, 'Compression', 'lzw');
        end
    end
end
