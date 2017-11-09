 function branch_analysis_combined_masks_movieList(ML)
% function to combine cell masks for filament segmentation
% Liya Ding,Jan, 2015
%
% Input:
%   ML:     The movieList object loaded before running this function

% the number of movies
movieNumber =  length(ML.movieDataFile_);
BA_output_ML_cell= cell(1,1);

for iM  = 1:movieNumber
    
    clearvars -except 'movieNumber' 'iM' 'ML' 
    
    close all;
    
    % load this movie
    load(ML.movieDataFile_{iM});
    % the number of channels
    display('======================================================================');
    display(['iM:', num2str(iM)]);
    combineChannelCellMaskCell = combineChannelMarkedCellAreaMask(MD);
    
    mkdir([MD.outputDirectory_,filesep,'combinedCellMasks']);
    for iFrame = 1 : MD.nFrames_
        if(~isempty(combineChannelCellMaskCell{iFrame}))
            imwrite(combineChannelCellMaskCell{iFrame},[MD.outputDirectory_,filesep,'combinedCellMasks',filesep,'combined_cell_masks_f',num2str(iFrame,'%03d'),'.tif']);
        end
    end    
end
