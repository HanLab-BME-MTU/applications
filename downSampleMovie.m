function downSampleMovie(experiment,subSample)

% determinePitDensities calculates density of pits in whole movie
% 
% INPUT:    experiment    =   struct containing all data for one condition 
%                             (e.g. clathrin control);in particular must
%                             have source field with the address for each
%                             cell folder
%
%           subSample    = time in seconds that you would like new 
%                           framerate to be (e.g. 2 if you would like to
%                           subsample new movie to 2 second framerate)
%
% OUTPUT:   none
%
% Note: Only use this function after running runLocalization on original
% data; the image files are not copied over
%
% D Nunez 2011-02-22


if round(subSample/experiment(1).framerate) ~= subSample/experiment(1).framerate
    error('new sampling rate must be multiple of old sampling rate')
end

%for each movie
for iexp = 1:length(experiment)
    
    %find tif files
    tifFiles = dir([experiment(iexp).source filesep '*.tif']);
    
    %load Detection results
    load([experiment(iexp).source 'Detection' filesep 'detectionResults.mat']);
    
    frameInfo = frameInfo(1:subSample/experiment(iexp).framerate:length(frameInfo));
    
    %until I figure out a better way I will cd to the date directory and
    %make a new cell folder there that is labeled as subsampled
    cd(experiment(iexp).source)
    cd ..
    %make sub sampled cell directory
    mkdir([experiment(iexp).source(1:end-1) '_subSampled_' num2str(subSample) 's'])
    %make Detection directory under this one
    mkdir([experiment(iexp).source(1:end-1) '_subSampled_' num2str(subSample) 's' ...
        filesep 'Detection'])
    %save sub sampled detection results
    save([experiment(iexp).source(1:end-1) '_subSampled_' num2str(subSample) 's' ...
        filesep 'Detection' filesep 'detectionResults.mat'],'frameInfo');
    %move one tif file
    movefile([experiment(iexp).source filesep tifFiles(1).name],[experiment(iexp).source(1:end-1) '_subSampled_' num2str(subSample) 's' ]);
    
end %of for each movie

end %of function