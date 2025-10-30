% averaging z-stack and exporting it to tiff from MD
function getMeanLayerFromStack(MD)
% Get channels
nChans = numel(MD.channels_);
nFrames = MD.nFrames_;
nZSlices = MD.zSize_;

for ii=nChans:-1:1
    % per channel
    curChan = MD.channels_(ii);
    % get average per frame
    for jj=1:nFrmaes
        
    end
end

end
