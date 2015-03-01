function [Length_pool,NMS_pool,ST_pool,INT_pool] = load_whole_movie_stat_MD(movieData,iChannel)
% load previous made pools for whole movie
% Slimest code, take it for granted that the filament segmentation is run
% in typical way.

Length_pool = [];
NMS_pool = [];
ST_pool = [];
INT_pool = [];

% Output Directories
% default filament segmention process output dir
FilamentSegmentationProcessOutputDir = [movieData.outputDirectory_, filesep 'FilamentSegmentation'];

if (~exist(FilamentSegmentationProcessOutputDir,'dir'))
    display('No previous stat done for whole movie');
    return;
end

% default filament segmention channel output dir
FilamentSegmentationChannelOutputDir = [FilamentSegmentationProcessOutputDir,filesep,'Channel',num2str(iChannel)];

if(~exist(FilamentSegmentationChannelOutputDir,'dir'))
    display('No previous stat done for whole movie in this channel');
    return;
end

% load the pools
load([FilamentSegmentationChannelOutputDir, filesep,...
    'pool_whole_movie_stat_channel', num2str(iChannel),'.mat'],...
    'Length_pool','NMS_pool','ST_pool','INT_pool');



  
