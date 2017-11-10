function BA_output_cell = branch_analysis_movieList(ML,figure_flag)
% function to do branch analysis for a whole movielist
% Liya Ding, March, 2014
%
% Input:
%   ML:     The movieList object loaded before running this function

% the number of movies
movieNumber =  length(ML.movieDataFile_);
BA_output_cell= cell(1,1);

for iM  = 1 :movieNumber
    % load this movie
    load(ML.movieDataFile_{iM});
    % the number of channels
    nChannel = length(MD.channels_);
    
     ROOT_DIR = MD.outputDirectory_;        
    
    for iChannel = 1 : nChannel        
        for iCell = 1 : 10
            
            display(['iM:', num2str(iM), ', iChannel:', num2str(iChannel),', iCell:', num2str(iCell)]);
            
            FilamentAnalysisPackage_complete_frames_file_name = [ROOT_DIR,'\','FilamentAnalysisPackage','\completedFramesChannel',num2str(iChannel),'Cell',num2str(iCell),'\completedFrames.mat'];
            SegmentationPackage_complete_frames_file_name = [ROOT_DIR,'\','SegmentationPackage','\completedFramesChannel',num2str(iChannel),'Cell',num2str(iCell),'\completedFrames.mat'];
            
            if(exist(FilamentAnalysisPackage_complete_frames_file_name,'file'))
                PackageName = 'FilamentAnalysisPackage';
            end
            
            if(exist(SegmentationPackage_complete_frames_file_name,'file'))
                PackageName = 'SegmentationPackage';
            end
            
            % the folder name if there is marking
            truthPath = [ROOT_DIR,'\',PackageName,'\FixedChannel',num2str(iChannel),'Cell',num2str(iCell)];
            % check if this folder exist
            if(exist(truthPath,'dir'))
                % if it exist, try to do the branch analysis
%                 try
                    BA_output = branch_analysis_marked_cell(MD, iChannel, iCell,figure_flag);
                    % for the output, ignore the channel number since
                    % there is only one channel marked.
                    BA_output_cell{iM, iCell} = BA_output;
%                 end
            end
        end        
    end
end

 ML_ROOT_DIR = ML.outputDirectory_;
 save([ML_ROOT_DIR,'\movieList_BA_output.mat'],'BA_output_cell');

 