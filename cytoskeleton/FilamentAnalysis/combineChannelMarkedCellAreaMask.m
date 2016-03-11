function combineChannelCellMaskCell = combineChannelMarkedCellAreaMask(MD)
% function to get cell mask from marked single cells, combining all channels
% Liya Ding, Jan, 2015
%
% Input:
%   MD:     The movieList object loaded before running this function

% initialize output, as a cell of all empty []
combineChannelCellMaskCell  = cell(1, MD.nFrames_);
for iFrame = 1 : MD.nFrames_
    combineChannelCellMaskCell{iFrame}=[];
end

% get the per channel per frame result
cellMaskCell = markedCellAreaMask(MD);


for iChannel = 1 : numel(MD.channels_)
    for iFrame = 1 : MD.nFrames_
        if(~isempty(cellMaskCell{iChannel, iFrame}))
            
            %if never anything, put zero image
            if(isempty(combineChannelCellMaskCell{iFrame}))
                combineChannelCellMaskCell{iFrame} = zeros(MD.imSize_);
            end
            
            % combine
            combineChannelCellMaskCell{iFrame} = (combineChannelCellMaskCell{iFrame} + cellMaskCell{iChannel, iFrame}) >0;
            
        end
    end
end
