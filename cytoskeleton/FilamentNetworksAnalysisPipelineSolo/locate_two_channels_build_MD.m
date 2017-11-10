function MD = locate_two_channels_build_MD(folder_name_channel_1,folder_name_channel_2,analysis_folder)
% locate the two channels and build the MD object and file
% Input: folder_name_channel_1: the folder of image of channel1
%        folder_name_channel_2: the folder of image of channel2
%        note: these two folders need to have same number of image with
%        same dimensions. (No tif recognized, need to be tiff)
%        analysisFolder: the folder for MD file and analysis output
% Output: the MD object, which is also save to disk, in the analysisFolder


%% locate the two channel data
channels_obj_one = Channel(folder_name_channel_1);
channels_obj_one.getImageFileNames();
channels_obj_one.sanityCheck();

channels_obj_two = Channel(folder_name_channel_2);
channels_obj_two.getImageFileNames();
channels_obj_two.sanityCheck();

channels_obj_all = [channels_obj_one channels_obj_two];


%% if no output dir, mk one
if(~exist(analysis_folder,'dir'))
    mkdir(analysis_folder);
end

%% build MD

MD = MovieData(channels_obj_all,analysis_folder);
MD.setPath(analysis_folder);
MD.setFilename('movieData.mat')
MD.sanityCheck();