function [allResults]=pitClusterRipleyTWindow_project(data, restrict, udist, windowsize, stepiter, scale)
% calculate the Ripley clustering of a lifetime sub-population defined by 
% the restriction vector projected across a given time window in the movie
% (used for CCP appearances and disappearances); the clustering is averaged 
% for a sliding time window across a given movie, and then averaged over 
% all available movies
% SYNOPSIS: [allResults]=pitClusterRipleyTWindow_project(data, restrict, udist, scale, windowsize, stepiter)
%
% INPUT:   
% data      = data structure pointing to all the image/detection/tracking
%           data for a given condition; for e.g. 10 cell movies, the 
%           structure should have 10 entries; the data structure can be 
%           created with the function loadConditionData, and it needs to 
%           contain at least the field 
%           .source, which is the (path) location of the lifetime information folder
%           .framerate, which is the movie framerate (necessary for the 
%           lifetime restriction
%           .movieLength, i.e. the number of frames in the movie
%
% restrict  = restriction vector can have variable length;
%           minimum length is five, where the entries are
%           [stat da minfr minlft maxlft]
%           optionally, the length can be extended to nine,
%           where the additional entries are 
%           [... minint maxint minmot maxmot], see below  
% 
% udist     = distance vector for Ripley function, e.g. [1:20]
%
% windowsize  = size of tau time window(s) (in frames); this is the time 
%           window across which the points are projected. This input can be
%           a single value or a vector, e.g. 40, or [50:50:300]. 
%           NOTE that the maximum possible value is the movie length; the
%           minimum possible value is 1 frame, but the minimum useful value
%           is determined by the frequency of the clustering process
% stepiter   = step size for sliding movement of the projection window (in
%           frames) e.g. 20 
%           Example: If one particular movie length is 300 fr, windowsize 
%           50 fr, and the stepsize 20 fr, then the clustering is averaged
%           for this movie for the 13 possible positions of the projection
%           window ranging from [1:50], [21:70], [41:90] up to [241:290]
%           The step size determines the degree of inter-movie averaging; 
%           the maximum useful step size is stepiter=windowsize, the 
%           minimum useful step size is stepiter=1 (although this is of 
%           course computationally expensive)
%
% scale     (OPTIONAL) distance scaling, e.g. 0.067 (meaning 1 pixel
%           corresponds to 0.067 um); DEFAULT value = 1
%
% OUTPUT:   
% allResults      = structure containing the fields:
%   .restriction    restrict vector for later reference
%   .dist           distance vector for later reference
%   .LR             all LR functions
%   .LRmean         LR averaged over all movies
%   .LRstd          std over all movies
%   .num            number of objects
%
%
% written by Dinah Loerke
% last modified June 29th, 2009


%% EXPLANATION of restriction values:
% rest = [stat da minfr minlft maxlft minint maxint minmot maxmot]
%                     
% stat  =   object status, possible values are 1, 2, 3
%           1: object appears and disappears in the movie
%           2: object is present for entire length of the movie
%           3: object is present in either first or last frame of the movie
% da    =   appearance/disappearance status, possible values are 1,0,-1
%           1: marks first frame
%           0: marks middle frame
%           -1: marks last frame
% minfr =   minimum lifetime/length in FRAMES - e.g. 4, to exclude tracking
%           artifacts of false detection positives
% minlft =  minimum lifetime in SECONDS - e.g. 60 to select for productive
%           clathrin-coated pits
% maxlft =  maximum lifetime in SECONDS - e.g. 25 to select for abortive
%           clathrin-coated pits
%
% OPTIONAL:
%
% minint =  minimum normalized intensity (ranging from 0 to 1)
% maxint =  maximum normalized intensity (ranging from 0 to 1)
% minmot =  minimum normalized motility (ranging from 0 to 1)
% maxmot =  maximum normalized motility (ranging from 0 to 1)
%
% These latter criteria allow you to select e.g. the brightest 10% of the
% population, or the faster 50%.
% 
% EXAMPLE:  to select the positions where productive pits appear
%           rest1 = [1 1 4 60 300]
%           to select the positions where abortive pits are located in each
%           frame
%           rest2 = [1 0 4 8 25]
%
% EXAMPLE of function call:
% [res] = pitClusterRipleyTWindow_project(dataControl, [1 1 4 60 300], [1:30], [50:50:300], 20, 0.067);



%%
% number of restriction 'sub-populations'
[ccx,ccy] = size(restrict);
% total number of self- and cross-correlations
numcc = sum([1:ccx]);

% allLR
allLR = zeros(length(udist),numcc,length(data));

orDir = cd;

% loop over all fields in the data
for i=1:length(data)
    
    fprintf('movie #%02d',i);
    
    path = data(i).source;
    cd(path);
    
    
    % load lifetime data
    lftInfo = [];
    % check for Lifetime Info data
    if exist('LifetimeInfo')==7
        cd('LifetimeInfo')
        if exist('lftInfo.mat')==2
            loadfile = load('lftInfo.mat');
            lftInfo = loadfile.lftInfo;
        end
        cd(path);
    end

    % if no leftInfo data is found, skip this movie
    if isempty(lftInfo)
        allLR(:,:,i) = nan;
        continue
    end
    
    % status matrix
    mat_stat    = lftInfo.Mat_status;
    % lifetime matrix
    mat_lft     = lftInfo.Mat_lifetime;
    % x-coordinate matrix
    mat_x       = lftInfo.Mat_xcoord;
    % y-coordinate matrix
    mat_y       = lftInfo.Mat_ycoord;
    % disapp status matrix
    mat_da      = lftInfo.Mat_disapp;
    % framerate
    fr          = data(i).framerate;
    % image size
    imsiz       = data(i).imagesize;
    
    msx = imsiz(1);
    msy = imsiz(2);
    imsizS = [imsiz(2) imsiz(1)];
       
    
    % construct convex hull out of complete point distribution
    % combined mpm
    selx = full(mat_x); selx = selx(isfinite(selx)); selx = nonzeros(selx(:));
    sely = full(mat_y); sely = sely(isfinite(sely)); sely = nonzeros(sely(:));
    combMPM = [ selx sely ];
    K = convhull(combMPM(:,1),combMPM(:,2));
    % edge points of the convex hull
    cpointsx = combMPM(K,1);
    cpointsy = combMPM(K,2);
    % create mask
    areamask = poly2mask(cpointsx,cpointsy,msx,msy);
    
    
    % CREATE CORRECTION FACTOR MATRIX FOR THIS MOVIE using all objects
    corfacmatM = makeCorrFactorMatrix(imsizS, udist, 10, areamask'); 
    normArea = sum(areamask(:));
    
    cc_counter = 1;
    
   
        
    % for the current value of restvector, collect those locations of 
    % points that fulfill the specified restrictions, regarding
    % appearance/disappearance status and minimum/maximum lifetime

    %desired minimum lifetime in seconds
    dstat   = restrict(1);
    dapp    = restrict(2);
    dminfr  = restrict(3);
    minlft  = restrict(4);
    maxlft  = restrict(5);
    minlft_fr = round(minlft/fr);
    maxlft_fr = round(maxlft/fr);

    [kx,ky]=size(mat_da);

    % status vectors for the different restriction conditions
    if isfinite(dstat)
        svec_stat = (mat_stat==dstat);
    else
        svec_stat = (mat_stat==mat_stat);
    end

    svec_da     = (mat_da==dapp);
    svec_lft1   = (mat_lft>dminfr);
    svec_lft2   = (mat_lft>minlft_fr);
    svec_lft3   = (mat_lft<maxlft_fr);
    
    % positions where the detected objects fulfill the requirements
    findpos = find( svec_stat & svec_da &  svec_lft1 & svec_lft2 & svec_lft3);
    % convert positions back into row,column format
    [findpos_row,findpos_column] = ind2sub(size(svec_da),findpos);

    % loop across all specified values for the window size
    for w=1:length(windowsize)
        
        % current windowsize
        cwin = windowsize(w);
        fprintf(' w=%04d',cwin);
        
        % tvec is a vector containing the start positions of the projection
        % window determined by the step size and the window size (the last
        % projection window still has to fit inside the movie length;
        % partial windows are not used)
        tvec = [1:stepiter:(ky-cwin+1)];
        
        % clear/initialize results matrices
        krMat_t = [];
        lrMat_t = [];
        prMat_t = [];
        num_t   = [];
        
        % loop across all start positions of the projection window        
        for t=1:length(tvec)
            
            fprintf(' t=%02d',t);
            
            % current start position
            tstart = tvec(t);
            % current end position
            tend = tstart+cwin-1;
            
            % findpos contains the positions where the detected objects
            % fulfill the requirements specified by the restriction vector;
            % these positions are now further restricted to include only
            % the objects occuring in the frames (columns) that belong to
            % this projection window
            fpos_t = find( (findpos_column>=tstart) & (findpos_column<=tend));
            % number of objects in this sub-window
            num_t(t) = length(fpos_t);
            
            % x- and y-positions of these objects
            findx_t = full(mat_x(findpos(fpos_t)));
            findy_t = full(mat_y(findpos(fpos_t)));
            currMPM_t = [findx_t findy_t];
            
            % plot(currMPM_t(:,1),currMPM_t(:,2),'b.');  pause(0.1);
            
            % calculate clustering
            [kr_t,lr_t,pr_t] = RipleysKfunction(currMPM_t, currMPM_t, imsizS, udist, corfacmatM, normArea);
            % write the results into results matrix
            krMat_t(t,:) = kr_t;
            lrMat_t(t,:) = lr_t;
            prMat_t(t,:) = pr_t;

            fprintf('\b\b\b\b\b\b');
        
        end


        % average across the different start positions of the projection
        % window (written into successive rows in the matrix)
        ckr = nanmean(krMat_t,1);
        clr = nanmean(lrMat_t,1);
        cpr = nanmean(prMat_t,1);
        cnum = nanmean(num_t);
        
        % write result into global matrix (that collects the results for
        % the different movies)
        KRmat_w(i,:,w) = ckr;
        LRmat_w(i,:,w) = clr;
        PRmat_w(i,:,w) = cpr;
        num_vec(i,w) = cnum;
        
        fprintf('\b\b\b\b\b\b\b');
        
    end % of w-loop
                 
    fprintf('\n');
    
end % of for-i loop


% average across the different movies
for n=1:size(KRmat_w,3)
    
    KRmat_ave(n,:) = nanmean(KRmat_w(:,:,n),1);
    LRmat_ave(n,:) = nanmean(LRmat_w(:,:,n),1);
    PRmat_ave(n,:) = nanmean(PRmat_w(:,:,n),1);
    
    KRmat_std(n,:) = nanstd(KRmat_w(:,:,n),1);
    LRmat_std(n,:) = nanstd(LRmat_w(:,:,n),1);
    PRmat_std(n,:) = nanstd(PRmat_w(:,:,n),1);
    
end



allResults.restriction  = restrict;
allResults.dist         = udist;
allResults.KRmat        = KRmat_w;

allResults.KRmean   = KRmat_ave;
allResults.KRstd    = KRmat_std;
allResults.LRmean   = LRmat_ave;
allResults.LRstd    = LRmat_std;
allResults.NDmean   = PRmat_ave;
allResults.NDstd    = PRmat_std;
allResults.num      = num_vec;

sc = 1;
if nargin>5
    if ~isempty(scale)
        sc = scale;
    end
end


%% display results - add your own favorite display code as necessary here
figure; hold on;
plot(sc*udist,PRmat_ave)

% return to original directory
cd(orDir);


end % of function
       