function [allResults]=pitClusterRipleyCross_project(data, restrict, udist, scale)
% calculate the Ripley clustering and cross-clustering of two populations 
% defined by the restriction vector projected over all frames of the
% movie, then averaging over all movies
% SYNOPSIS: [allResults]=pitClusterRipleyCross_project(data, udist, restrict)
%
% INPUT:   data     =   data structure pointing to all the
%                       image/detection/tracking data for a given condition;
%                       for e.g. 10 cell movies, the structure should have 10
%                       entries; the data structure can be created with the 
%                       function loadConditionData, and it needs to contain
%                       at least the field 
%                       .source, which is the (path) location of the 
%                       lifetime information folder
%                       .framerate, which is the movie framerate, which is
%                       necessary for the lifetime restriction
%                       .movieLength, i.e. the number of frames in the
%                       movie                       
%           restrict  = restriction vector can have variable length;
%                       minimum length is five, where the entries are
%                       [stat da minfr minlft maxlft]
%                       optionally, the length can be extended to nine,
%                       where the additional entries are 
%                       [... minint maxint minmot maxmot], see below  
%                       NOTE: For calculating the cross-correlation of 
%                       two populations, the vector needs to have two rows,
%                       where each row constitutes the restriction for one
%                       population                       
%           udist   =   distance vector for Ripley function, e.g. [1:20]
%           scale   =   distance scaling (optional), e.g. 0.067 (um/pix)
%
% OUTPUT:   allResults      = structure containing the fields:
%           .restriction    restrict vector for later reference
%           .dist          distance vector for later reference
%           .LR             all LR functions
%           .LRmean         LR averaged over all movies
%           .LRstd          std over all movies
%           .LRsem          standard error of the mean over all movies
%           .dens            density of points of the two populations in
%                           each movie
%
% written by Dinah Loerke
% last modified January 28th, 2008
% last modified Feb 13th, 2008
%


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
% [res] = pitClusterRipleyCross_project(dataClathrinLatr,...
%         [[1 1 4 60 300];[1 1 4 8 25]], [1:20], 0.67);



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
    
    % loop over all rows of the restriction values, which designate
    % different conditions
    
    for r=1:ccx
        
        % for the current value of restvector, collect those locations of 
        % points that fulfill the specified restrictions, regarding
        % appearance/disappearance status and minimum/maximum lifetime

        %desired minimum lifetime in seconds
        dstat   = restrict(r,1);
        dapp    = restrict(r,2);
        dminfr  = restrict(r,3);
        minlft  = restrict(r,4);
        maxlft  = restrict(r,5);
        minlft_fr = round(minlft/fr);
        maxlft_fr = round(maxlft/fr);

        [kx,ky]=size(mat_da);
    
        findpos = find( (mat_stat==dstat) & (mat_da==dapp) &...
            (mat_lft>dminfr) & (mat_lft>minlft_fr) & (mat_lft<maxlft_fr));
        findx = full(mat_x(findpos));
        findy = full(mat_y(findpos));
        
        currMPM = [findx findy];
        restrictedMPMs(r).mpm = currMPM;
        
        % density of objects
        allDen(r,i) = 1000 * length(findpos)/sum(areamask(:));
        
        %calculate average lifetime of restricted objects for each
        %restriction and each movie
        LifetimesPerMovie = max(mat_lft,[],2);  %find lifetimes for each trajectory
        findRestrictedLifetimes = find((LifetimesPerMovie>minlft_fr) & (LifetimesPerMovie<maxlft_fr)); %find where lifetimes fall within specified restrictions
        avgLifetimeForMovie(r,i) = nanmean(LifetimesPerMovie(findRestrictedLifetimes)); %average all lifetimes that fall within specified restrictions
        
         % print the average number of points of each population in each 
         % frame - if this number is extremely low, this acts as a 
         % cautionary in interpreting the data
        print_np = length(findpos);
        fprintf([' nump',num2str(r)]);
        fprintf('=%05d',print_np);
            
    end
    
   
    
    
    % now loop over all combinations of the different restriction
    % conditions
    for k1=1:ccx
        
        currMPM1 = restrictedMPMs(k1).mpm;
                        
        % determine second mpm file
        for k2 = k1:ccx
            
            currMPM2 = restrictedMPMs(k2).mpm;
            
                       
            % determine mask for correction
            if (length(currMPM1)>20) & (length(currMPM2)>20)
                              
                % determine ripley clustering function for
                % cross-correlation of the two specified mpms (which may be
                % the same mpm in some cases)       

                % cross-correlation mpm1-mpm2
                %[krCross,lrCross]=RipleyKfunc_standalone(currMPM1, currMPM2, imsizS, udist, corfacmat);
                [krCross,lrCross]=RipleysKfunction(currMPM1, currMPM2, imsizS, udist, corfacmatM, normArea);
                
                % store for averaging
                allLR(:,cc_counter,i) = lrCross;

                                
                % cross_correlation combination
                ccComb(cc_counter,:) = [k1,k2];

                
                % display
                % hold on; plot(udist,lrCross,'b.-'); 
                % hold on; plot(udist,lrCross2,'r.-');
                % pause(0.1);
                
            else
                allLR(:,cc_counter,i) = nan;
                allDen(i) = nan;
            end % of if
            
            cc_counter = cc_counter+1;
            
                      
            
        end % of for k2
        
    end % of for k1
            
    fprintf('\n');
    
end




lrAV = nanmean(allLR,3);
lrSTD = nanstd(allLR,1,3);
lrSEM = lrSTD/sqrt(length(data));

allResults.restriction = restrict;
allResults.dist = udist;
allResults.LR = allLR;
allResults.LRmean = lrAV;
allResults.LRstd = lrSTD;
allResults.LRsem = lrSEM;
allResults.dens = allDen;
allResults.comb = ccComb;

 
areaplot1 = [ ((lrAV(:,1)-lrSEM(:,1))') ; 2*lrSEM(:,1)' ];
if ccx>1
    areaplot2 = [ ((lrAV(:,2)-lrSEM(:,2))') ; 2*lrSEM(:,2)' ];
    areaplot3 = [ ((lrAV(:,3)-lrSEM(:,3))') ; 2*lrSEM(:,3)' ];
end

pdist = udist;
if nargin>3
    pdist = scale*udist;
end

figure
hold on;
area(pdist,areaplot1');
plot(pdist,lrAV(:,1),'g-');
if ccx>1
    area(pdist,areaplot2');
    plot(pdist,lrAV(:,2),'r-');
    area(pdist,areaplot3');
    plot(pdist,lrAV(:,3),'c-');
end

%%=========================================================================
%           convert allLR to local densities
%==========================================================================


allDen = allLR;
[dlx,dly,dlz] = size(allLR);

carea = udist.^2;
careadiff = carea; careadiff(2:length(careadiff)) = diff(carea);
amat = repmat(careadiff',1,dly);

dmat = repmat(udist',1,dly);

for i=1:length(data)
    currLR = allLR(:,:,i);
    currKR = (currLR+dmat).^2;
    currKRdiff = currKR;
    currKRdiff(2:length(udist),:) = diff(currKR,1);
    currDen = currKRdiff./amat;
    
    allDen(:,:,i) = currDen;

end


ndenAV = nanmean(allDen,3);
ndenSTD = nanstd(allDen,1,3);
ndenSEM = ndenSTD/sqrt(length(data));

allResults.ND = allDen;
allResults.NDmean = ndenAV;
allResults.NDstd = ndenSTD;
allResults.NDsem = ndenSEM;

areaplot1d = [ ((ndenAV(:,1)-ndenSEM(:,1))') ; 2*ndenSEM(:,1)' ];
if ccx>1
    areaplot2d = [ ((ndenAV(:,2)-ndenSEM(:,2))') ; 2*ndenSEM(:,2)' ];
    areaplot3d = [ ((ndenAV(:,3)-ndenSEM(:,3))') ; 2*ndenSEM(:,3)' ];
end

figure
hold on;
area(pdist,areaplot1d');
plot(pdist,ndenAV(:,1),'g-');
if ccx>1
    area(pdist,areaplot2d');
    plot(pdist,ndenAV(:,2),'r-');
    area(pdist,areaplot3d');
    plot(pdist,ndenAV(:,3),'c-');
end


end % of function
       