function [filament_orientation_pool, network_feature_ML] = load_ML_network_for_VIF_anisotropy(ML,CellROI,radius,downsample_index)
% function to pooling filament orientation data together
% input:    ML:    the loaded movieList object.
%           CellROI:   the user input CellROI, if not defined [], will use the whole
%                  image area
%           radius: the definition of neighborhood
%           downsample_index: if 1, no downsampling; other wise 1 in
%                             downsample_index data will be taken.
% output:   filament_orientation_pool: the pool of filament orientation std
%           data for the whole movie list, for each channel
%           network_feature, a cell structure for each movies, then each channel, each frame.
%           Each struct with field of the 2 features as above

% Liya Ding 2015

% if no input, CellROI is full image
if(nargin<2)
    CellROI = [];
end

% if no input for radius, set it as default 20
if(nargin<3)
    radius = 20;
end

% if no downsample index, no downsampling
if(nargin<4)
    downsample_index = 1;
end

%%
% the number of movies
nMovie =  length(ML.movieDataFile_);

MD = MovieData.load(ML.movieDataFile_{1});
% all the movie should have the same number of channels, for same content
% each channel
nChannel = length(MD.channels_);
nFrame = MD.nFrames_;

% initialize output
network_feature_ML = cell(1,nMovie);
for iM = 1 : nMovie
    network_feature_ML{iM} =cell(nChannel, nFrame);    
end

filament_orientation_pool= cell(1,nChannel);
for iChannel = 1 : nChannel
    filament_orientation_pool{iChannel} =[];
end

% get analysis done
for iM  = 1 : nMovie
    % load this movie
    MD = MovieData.load(ML.movieDataFile_{iM});
    % the number of channels
    display('=============================================================');
    display(['Movie ', num2str(iM)]);
    display('Calculating:');
    network_feature = load_MD_network_for_VIF_anisotropy(MD,[],radius);
    network_feature_ML{iM} = network_feature;
end

% gather the pool for each channel
for iChannel = 1 : nChannel
    display('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&');
    display(['Channel:', num2str(iChannel)]);
    display('Pooling:');
    
    
    for iM  = 1 : nMovie
        display('=========================================================');
        display(['Movie ', num2str(iM)]);
        display('Pooling:');
        
        network_feature = network_feature_ML{iM};
        
        for iFrame = 1 : size(network_feature,2)
            
            if downsample_index > 1 % if need downsampled
                % get the feature data
                data = network_feature{iChannel,iFrame}.filament_orientation_STD(:);
                % get sampling index
                subsample_index = randsample(length(data), round(length(data)/downsample_index));
                
                %sample and add to pool
                filament_orientation_pool{iChannel} = [filament_orientation_pool{iChannel}; ...
                    data(subsample_index) ];
                
            else % no downsampling
                % add to pool
                if(~isempty(network_feature{iChannel,iFrame}))
                filament_orientation_pool{iChannel}  = [filament_orientation_pool{iChannel}; ...
                    network_feature{iChannel,iFrame}.filament_orientation_STD(:)];
                end
            end
        end
    end
    
end

save([ML.outputDirectory_,'wholemovie_eachch_filamentorient_STD.mat'],'filament_orientation_pool','network_feature_ML');



