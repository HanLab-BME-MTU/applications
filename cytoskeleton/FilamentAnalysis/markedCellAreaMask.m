function cellMaskCell = markedCellAreaMask(MD)
% function to get cell mask from marked single cells
% Liya Ding, Jan, 2015

% Input:
%   MD:     The movieList object loaded before running this function

% the number of channels
nChannel = numel(MD.channels_);
ROOT_DIR = MD.outputDirectory_;

% initialize channel-frame masks, as []
cellMaskCell = cell(numel(MD.channels_), MD.nFrames_);

for iChannel = 1 : nChannel
    for iFrame = 1 : MD.nFrames_
        cellMaskCell{iChannel, iFrame} = [];
    end
end


for iChannel = 1 : nChannel
    % For each channel check the existance of marked cells
    
    for iCell = 1 : 10
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
            % if both not exist, exit
            % for all possible cell number, all exited, then everything is
            % still []
            continue;
        end
        
        % the folder name if there is marking
        truthPath = [ROOT_DIR,filesep,PackageName,filesep,'OnlyFixedChannel',num2str(iChannel),'Cell',num2str(iCell)];
        old_truthPath = [ROOT_DIR,filesep,PackageName,filesep,'FixedChannel',num2str(iChannel),'Cell',num2str(iCell)];
        
        % check if this folder exist
        if(exist(truthPath,'dir') || exist(old_truthPath,'dir'))
            
            load([ROOT_DIR,filesep,PackageName,filesep,'completedFramesChannel',num2str(iChannel),'Cell',num2str(iCell),filesep,'completedFrames.mat'],'isCompleted');
            truthPath = [ROOT_DIR,filesep,PackageName,filesep,'OnlyFixedChannel',num2str(iChannel),'Cell',num2str(iCell)];
            old_truthPath = [ROOT_DIR,filesep,PackageName,filesep,'FixedChannel',num2str(iChannel),'Cell',num2str(iCell)];
            
            if  (sum(isCompleted)==0)
                return;
            end
            
            fileExistance = 0*isCompleted;
            
            for iFrame = 1 : MD.nFrames_
                if(exist([old_truthPath,filesep,'mask_',num2str(iFrame),'.tif'],'file') ...
                        || exist([truthPath,filesep,'mask_',num2str(iFrame),'.tif'],'file'))
                    fileExistance(iFrame)=1;
                end
            end
            
            isCompletedFileExist = (isCompleted.*fileExistance)>0;
            
            for iFrame = 1 : MD.nFrames_
                
                if(isCompletedFileExist(iFrame)>0)   
                    %% Load the mask
                    if(exist([truthPath,filesep,'mask_',num2str(iFrame),'.tif'],'file'))
                        current_mask_iCh_iCell_iF = imread([truthPath,filesep,'mask_',num2str(iFrame),'.tif']);
                    else
                        if(exist([old_truthPath,filesep,'mask_',num2str(iFrame),'.tif'],'file'))
                            current_mask_iCh_iCell_iF = imread([old_truthPath,filesep,'mask_',num2str(iFrame),'.tif']);
                        else
                            disp('no mask found');
                            break;
                        end
                    end
                    current_mask_iCh_iCell_iF([ 1 end], :)=0;
                    current_mask_iCh_iCell_iF(:, [ 1 end])=0;
                    current_mask_iCh_iCell_iF = keep_largest_area(current_mask_iCh_iCell_iF);
                    
                    % if never has content, assign zero image
                    if(isempty(cellMaskCell{iChannel,iFrame}))
                        cellMaskCell{iChannel,iFrame} = zeros(MD.imSize_)>0;
                    end
                    
                    % combine mask from all cells     
                    % this addition is done for all iCells
                    cellMaskCell{iChannel,iFrame} = (cellMaskCell{iChannel,iFrame} + current_mask_iCh_iCell_iF) >0;
                end
            end
        end
        
    end
end
