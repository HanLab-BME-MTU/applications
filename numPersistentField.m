function [numPers]=numPersistentField(exp, threshold)
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
% Dinah Loerke, July 11, 2007
% Last modified: Francois Aguet, 01/21/2010


% proceed if the movie is long enough for the specified threshold
% and has specified framerate
% AND to prevent multiple copies of resampled movies skewing the
% results, use only one of multiple resample copies
if ((exp.framerate*exp.movieLength)>min(threshold)) %& (framerate>=2)
        
    if ~isfield(exp,'lftInfo') || isempty(exp.lftInfo)
    
        lftPath = [exp.source filesep 'LifetimeInfo' filesep 'lftInfo.mat'];
        if (exist(lftPath, 'file')==2)
            load(lftPath);
        end;
    end;
            
    if ~isempty(lftInfo)
        % convert the threshold (which is in seconds) into number of frames,
        % which depends on the framerate of this movies
        thresh_frames = round(threshold/exp.framerate);
        
        if ( exp.framerate*exp.movieLength>threshold )
            numPers = numAboveThresh(lftInfo, thresh_frames);
        end
    end
end
end


function [numPers]=numAboveThresh(lftInfo, threshold)

numPers = zeros(1,length(threshold));

% calculate the number of objects in the field that have a lifetime
% exceeding the threshold value, regardless of status

Mlft = full(lftInfo.Mat_lifetime);
Mstat = full(lftInfo.Mat_status);

% initialize results vector
numCat = zeros(4,length(threshold));

% loop over all specified threshold values
for t = 1:length(threshold)
    
    for b = 1:size(Mlft,1)
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
        if (currLifetime > threshold(t))
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
    end
    
    % calculate fraction of category 1 out of the 'usable' trajectories
    numPers(t) = numCat(1,t);
end
end