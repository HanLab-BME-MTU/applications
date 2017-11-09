% load whole movie stat pools and gather them, calculate the stat for
% the whole movieList

% ML.sanityCheck();

% Assume the movie list consists movie with same setting in channels and everything else
MD =  MovieData.load(ML.movieDataFile_{1});
nChannel = numel(MD.channels_);
Whole_movie_stat_cell = cell(1,nChannel);


%% for length first

Group_length_pool = [];

for iChannel = 1 : nChannel
    
    %% Set everything to nan first
    Whole_movie_stat.mean_ST = nan;
    Whole_movie_stat.std_ST = nan;
    Whole_movie_stat.mode_ST = nan;
    Whole_movie_stat.otsu_ST = nan;
    Whole_movie_stat.otsu_mode_ST = nan;
    Whole_movie_stat.rosin_ST = nan;
    Whole_movie_stat.rosin_mode_ST = nan;
    
    
    Whole_movie_stat.mean_NMS = nan;
    Whole_movie_stat.std_NMS = nan;
    Whole_movie_stat.mode_NMS = nan;
    Whole_movie_stat.otsu_NMS = nan;
    Whole_movie_stat.otsu_mode_NMS = nan;
    Whole_movie_stat.rosin_NMS = nan;
    Whole_movie_stat.rosin_mode_NMS = nan;
    
    
    
    Whole_movie_stat.mean_INT = nan;
    Whole_movie_stat.std_INT = nan;
    Whole_movie_stat.mode_INT = nan;
    Whole_movie_stat.otsu_INT = nan;
    Whole_movie_stat.otsu_mode_INT = nan;
    
    Whole_movie_stat.mean_Length = nan;
    Whole_movie_stat.std_Length = nan;
    Whole_movie_stat.mode_Length = nan;
    Whole_movie_stat.rosin_Length = nan;
    Whole_movie_stat.rosin_mode_Length = nan;
    
    Whole_movie_stat_cell{iChannel} = [];
    
    %% Load movie by movie, just for length
    Group_length_pool=[];
    
    for iM = 1 : numel(ML.movieDataFile_);
        MD =  MovieData.load(ML.movieDataFile_{iM});
        [Length_pool,~,~ ,~] = load_whole_movie_stat_MD(MD,iChannel);
        if(~isempty(Length_pool))
            Group_length_pool = [Group_length_pool Length_pool];
        end
    end
    
    % if there is no information for this channel, leave it as empty
    if(isempty(Group_length_pool))
        continue;
    end
    
    Group_length_pool = Group_length_pool(Group_length_pool>1);
    [hist_n,bin] = hist(Group_length_pool,50);
    ind_mode = find(hist_n==max(hist_n));
    mode_Group_Length = bin(ind_mode(1));
    
    Whole_movie_stat.mean_Length = mean(Group_length_pool);
    Whole_movie_stat.std_Length = std(Group_length_pool);
    Whole_movie_stat.mode_Length = mode_Group_Length;
    Whole_movie_stat.rosin_Length = thresholdRosin(Group_length_pool);
    Whole_movie_stat.rosin_mode_Length = thresholdRosin(Group_length_pool(find(Group_length_pool>mode_Group_Length)));
    
    
    Group_length_pool=[];
    
    %% Load movie by movie, just for STNMS
    Group_NMS_pool=[];
    
    for iM = 1 : numel(ML.movieDataFile_);
        MD =  MovieData.load(ML.movieDataFile_{iM});
        [~,NMS_pool,~ ,~] = load_whole_movie_stat_MD(MD,iChannel);
        if(~isempty(NMS_pool))
            Group_NMS_pool = [Group_NMS_pool; NMS_pool];
        end
    end
    
    % if there is no information for this channel, leave it as empty
    if(isempty(Group_NMS_pool))
        continue;
    end
    
    [hist_n,bin] = hist(Group_NMS_pool,50);
    ind_mode = find(hist_n==max(hist_n));
    mode_Group_NMS = bin(ind_mode(1));
    
    Whole_movie_stat.mean_NMS = mean(Group_NMS_pool);
    Whole_movie_stat.std_NMS = std(Group_NMS_pool);
    Whole_movie_stat.mode_NMS = mode_Group_NMS;
    Whole_movie_stat.otsu_NMS = thresholdOtsu(Group_NMS_pool);
    Whole_movie_stat.otsu_mode_NMS = thresholdOtsu(Group_NMS_pool(find(Group_NMS_pool>mode_Group_NMS)));
    
    try
        Whole_movie_stat.rosin_NMS = thresholdRosin(Group_NMS_pool);
    catch
        Whole_movie_stat.rosin_NMS = Whole_movie_stat.otsu_NMS;
    end
    
    try
        Whole_movie_stat.rosin_mode_NMS = thresholdRosin(Group_NMS_pool(find(Group_NMS_pool>mode_Group_NMS)));
    catch
        Whole_movie_stat.rosin_mode_NMS = Whole_movie_stat.otsu_mode_NMS;
    end
    
    Group_NMS_pool=[];
    
    % ST
    % Load movie by movie, just for ST
    Group_ST_pool=[];
    
    for iM = 1 : numel(ML.movieDataFile_);
        MD =  MovieData.load(ML.movieDataFile_{iM});
        [~,~,ST_pool ,~] = load_whole_movie_stat_MD(MD,iChannel);
        if(~isempty(ST_pool))
            Group_ST_pool = [Group_ST_pool; ST_pool];
        end
    end
    
    % if there is no information for this channel, leave it as empty
    if(isempty(Group_ST_pool))
        continue;
    end
    
    [hist_n,bin] = hist(Group_ST_pool,50);
    ind_mode = find(hist_n==max(hist_n));
    mode_Group_ST = bin(ind_mode(1));
    
    Whole_movie_stat.mean_ST = mean(Group_ST_pool);
    Whole_movie_stat.std_ST = std(Group_ST_pool);
    Whole_movie_stat.mode_ST = mode_Group_ST;
    Whole_movie_stat.otsu_ST = thresholdOtsu(Group_ST_pool);
    Whole_movie_stat.otsu_mode_ST = thresholdOtsu(Group_ST_pool(find(Group_ST_pool>mode_Group_ST)));
    
    try
        Whole_movie_stat.rosin_ST = thresholdRosin(Group_ST_pool);
    catch
        Whole_movie_stat.rosin_ST = Whole_movie_stat.otsu_ST;
    end
    
    try
        Whole_movie_stat.rosin_mode_ST = thresholdRosin(Group_ST_pool(Group_ST_pool>mode_Group_ST));
    catch
        Whole_movie_stat.rosin_mode_ST = Whole_movie_stat.otsu_mode_ST;
    end
    
    Group_ST_pool=[];
    
    %% INT
    
    
    Group_INT_pool=[];
    
    for iM = 1 : numel(ML.movieDataFile_);
        MD =  MovieData.load(ML.movieDataFile_{iM});
        [~,~,~,INT_pool] = load_whole_movie_stat_MD(MD,iChannel);
        if(~isempty(INT_pool))
            Group_INT_pool = [Group_INT_pool INT_pool];
        end
    end
    
    % if there is no information for this channel, leave it as empty
    if(isempty(Group_INT_pool))
        continue;
    end
    
    Group_INT_pool =  double(Group_INT_pool(:));
    
    % find the mode of the intensity of the curves/lines
    [hist_n,bin] = hist((Group_INT_pool),50);
    ind_mode = find(hist_n==max(hist_n));
    mode_Group_INT = bin(ind_mode(1));
    
    Whole_movie_stat.mean_INT = mean(Group_INT_pool);
    Whole_movie_stat.std_INT = std(Group_INT_pool);
    Whole_movie_stat.mode_INT = mode_Group_INT;
    Whole_movie_stat.otsu_INT = thresholdOtsu(Group_INT_pool);
    Whole_movie_stat.otsu_mode_INT = thresholdOtsu(Group_INT_pool(find(Group_INT_pool>mode_Group_INT)));
    Group_INT_pool=[];
    %%
    Whole_movie_stat_cell{iChannel} = Whole_movie_stat;
end

save([ML.outputDirectory_,filesep,'whole_movie_list_stat.mat'],'Whole_movie_stat_cell');
