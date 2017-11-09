function BA_output_cell = branch_analysis_movieData(MD,half_size,min_branch_size_Threshold,filament_stat_flag,figure_flag)
% function to do branch analysis for a whole movieData
% Liya Ding, March, 2014
%
% Input:
%   MD:     The movieList object loaded before running this function


if(nargin<2)
    half_size=150;
end

if(nargin<3)
    min_branch_size_Threshold=100;
end

if(nargin<4)
    filament_stat_flag=0;
end

if(nargin<5)
    figure_flag=0;
end
BA_output_cell = cell(1,1);

% the number of channels
nChannel = length(MD.channels_);

ROOT_DIR = MD.outputDirectory_;

for iChannel = 1 : nChannel
    for iCell = 1 : 20        
        display(['iChannel:', num2str(iChannel),', iCell:', num2str(iCell)]);
        
        FilamentAnalysisPackage_complete_frames_file_name = [ROOT_DIR,filesep,'FilamentAnalysisPackage',filesep,'completedFramesChannel',num2str(iChannel),'Cell',num2str(iCell),filesep,'completedFrames.mat'];
        SegmentationPackage_complete_frames_file_name = [ROOT_DIR,filesep,'SegmentationPackage',filesep,'completedFramesChannel',num2str(iChannel),'Cell',num2str(iCell),filesep,'completedFrames.mat'];
        PackageName=[];
        
        if(exist(FilamentAnalysisPackage_complete_frames_file_name,'file'))
            PackageName = 'FilamentAnalysisPackage';
        end
        
        if(exist(SegmentationPackage_complete_frames_file_name,'file'))
            PackageName = 'SegmentationPackage';
        end
        
        if(isempty(PackageName))
            continue;
        end
        
        % the folder name if there is marking
        truthPath = [ROOT_DIR,filesep,PackageName,filesep,'OnlyFixedChannel',num2str(iChannel),'Cell',num2str(iCell)];
        old_truthPath = [ROOT_DIR,filesep,PackageName,filesep,'FixedChannel',num2str(iChannel),'Cell',num2str(iCell)];
        % check if this folder exist
        if(exist(truthPath,'dir') || exist(old_truthPath,'dir'))
            % if it exist, try to do the branch analysis
%             try
            BA_output = branch_analysis_marked_cell(MD, iChannel, iCell,half_size,min_branch_size_Threshold,filament_stat_flag,figure_flag);
%             catch
%                 display_msg = 'This cell is not working, check back later!'
%             end
            % for the output, ignore the channel number since
            % there is only one channel marked.
            if(~isempty(BA_output))
                BA_output_cell{iChannel, iCell} = BA_output;            
            end
        end
    end
end


save([ROOT_DIR,filesep,'movieData_BA_output_balloon.mat'],'BA_output_cell');

