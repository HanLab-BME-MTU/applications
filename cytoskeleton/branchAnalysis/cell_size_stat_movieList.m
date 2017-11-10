function [cell_size_matrix, total_vim_matrix] = cell_size_stat_movieList(ML)
% function to do calculate the cell size for a whole movielist
% Liya Ding, April, 2014
%
% Input:
%   ML:     The movieList object loaded before running this function

% the number of movies
movieNumber =  length(ML.movieDataFile_);

cell_size_matrix = nan(movieNumber,10);
total_vim_matrix = nan(movieNumber,10);

for iM  = 1 :movieNumber
    % load this movie
    load(ML.movieDataFile_{iM});
    % the number of channels
    nChannel = length(MD.channels_);
    
     ROOT_DIR = MD.outputDirectory_;        
     
     if(exist([ROOT_DIR,'\FilamentAnalysisPackage\refined_masks\'],'dir'))
         PackageName = 'FilamentAnalysisPackage';
     else
         PackageName = 'SegmentationPackage';
     end
    
    for iChannel = 1 : nChannel        
        for iCell = 1 : 10
            % the folder name if there is marking
            truthPath = [ROOT_DIR,'\',PackageName,'\FixedChannel',num2str(iChannel),'Cell',num2str(iCell)];
            % check if this folder exist
            if(exist(truthPath,'dir'))
                % if it exist, try to do the branch analysis
%                 try
                    [cell_size_mean, total_vim_mean] = ...
                        cell_size_after_branch_analysis_marked_cell(MD, ...
                        iChannel, iCell);
                    % for the output, ignore the channel number since
                    % there is only one channel marked.
                    cell_size_matrix(iM, iCell) = cell_size_mean;
                    total_vim_matrix(iM, iCell) = total_vim_mean;
%                 end
            end
        end        
    end
end

 ML_ROOT_DIR = ML.outputDirectory_;
 save([ML_ROOT_DIR,'\cell_size_total_vim_ML.mat'],'cell_size_matrix','total_vim_matrix');

 