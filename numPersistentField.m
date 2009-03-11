function [numPers]=numPersistentField(experField, threshold)
% calculate the number of objects that have a lifetime 
% exceeding a specified threshold value, so-called 'persistent' objects
% INPUT: experFiedld      = experiment structure field
%        threshold  = threshold in seconds
% OUTPUT: numPers     = fraction of long-lived objects in every frame
%         
% NOTE: only movies with framerate 2 are used for this analysis
% NOTE: in a previous version of this function, the fraction of persistent 
% objects in each frame was calculated, this function looks for the
% fraction of the whole pool of objects - because of the higher frequency
% of short-lived objects, the contribution of persistent objects per whole
% will be lower!
%
% last modified : July 11, 2007 (Dinah)

od=cd;

numPers = [];

% read current framerate
framerate = experField.framerate;
    
% read number of frames in movie
numf = experField.movieLength;
    
% proceed if the movie is long enough for the specified threshold
% and has specified framerate
% AND to prevent multiple copies of resampled movies skewing the
% results, use only one of multiple resample copies
if ((framerate*numf)>min(threshold)) %& (framerate>=2)
    
        
        
    % read path
    path = experField.source;
    
    % load TrackInfo matrix
    loadvar = 0;
    if isfield(experField,'lftInfo')
        lftInfo = experField.lftInfo;
        if isempty(lftInfo)
            loadvar = 1;
        end
    else
        loadvar = 1;
    end
    
    if loadvar==1
        odir = cd;
        cd(path);
        if exist('LifetimeInfo')==7
            cd('LifetimeInfo')
            loadfile = load('lftInfo.mat');
            lftInfo = loadfile.lftInfo;
        end
        cd(odir);
    end
    
    
    if ~isempty(lftInfo)
        % convert the threshold (which is in seconds) into number of frames,
        % which depends on the framerate of this movies
        thresh_frames = round(threshold/framerate);
         
        if ( framerate*numf>threshold )
            numPers = numAboveThresh(lftInfo, thresh_frames);
        end
    end
    
           
       
end % of if movie is long enough
    
    
    
end % of function



function [numPers]=numAboveThresh(lftInfo, threshold)

% calculate the number of objects in the field that have a lifetime
% exceeding the threshold value, regardless of status


% first calculate lifetime and status matrices
Mlft = lftInfo.Mat_lifetime;
Mstat = lftInfo.Mat_status;

Mlft = full(Mlft);

[sx, sy] = size(Mlft);

% initialize results vector
numCat = zeros(4,length(threshold));

% loop over all specified threshold values
for t = 1 : length(threshold)
    th = threshold(t);
    
    
    % loop over all objects
    for b=1:sx
        % objects are divided into several categories:
        % 1. persistent (lifetime > th, regardless of status)
        % 2. analysis pool (full lifetime is known, i.e. status =1, and
        %    lifetime is longer than 4 frames and less than th)
        % 3. artifactual (lifetime is >=4 frames, therefore not used for
        %    the analysis of the different populations)
        % 4. unknown (full lifetime isn't known AND >th)
        %
        % The final fraction of persistent pits is the contributions of
        % cat1/(cat1+cat2)
        
        % current object lifetime 
        currLifetime = max(Mlft(b,:));
        % current object status 
        currStatus = min(nonzeros(Mstat(b,:)));
        
        % if lifetime > th, then objevcts is persistent
        if (currLifetime>th)
            cat = 1;
        else
        % else if the trajectory is partial, the lifetime is unknown
            if (currStatus>1)
                cat = 4;
        % else the trajectory is known, and the category depends only on
        % whether the lifetime is above or below the artifact threshold of
        % 4 (frames)
            else
                if (currLifetime>4)
                    cat = 2;
                else
                    cat = 3;
                end
            end
        end
        
        % add to the appropriate category
        numCat(cat,t) = numCat(cat,t)+1;
        
        
    end % of for b
    
    % calculate fraction of category 1 out of the 'usable' trajectories
    numPers(t) = numCat(1,t);
    
end % of for t

end % of function

