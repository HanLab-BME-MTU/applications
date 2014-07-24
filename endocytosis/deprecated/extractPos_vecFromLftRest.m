function [posvec] = extractPos_vecFromLftRest(lftinfopath, restvector, framerate)
% this function extracts an MPM of points of interest from trackInfo using
% specified restrictions, and time delays in the third dimension
% 
% SYNOPSIS [MPMglobal] = CorrelateData2Pos_extractMPM(trackinfo,restvector, tvec, ttype);
%
% INPUT     lftinfopath       
%
%           restvector      = restriction vector, defined by
%               dstat       = restvector(1) = status [1,2,3];
%               dapp        = restvector(2) = disappearance [-1,0,1];
%               dminfr      = restvector(3) = min lifetime in frames;
%               minlft      = restvector(4) = min lifetime in seconds;;
%               maxlft      = restvector(5) = min lifetime in seconds;;
%               minlft_fr   = round(minlft/fr);
%               maxlft_fr   = round(maxlft/fr);
%               
%               OPTIONAL:
%               intrank_min = restvector(6)  = minimum intensity rank
%                             percentile, e.g. 0.80
%               intrank_max = restvector(7)  = maxmimum intensity rank
%                             percentile, e.g. 0.95
%
%           fr              = framerate; default 1
%
%
% OUTPUT:   posvec       = result MPM
%
% Dinah Loerke   01/30/2009
% Last modified: Francois Aguet 02/10/2010

if nargin<3 || isempty(framerate)
    framerate = 1;
end

% calculate lft matrices
load([lftinfopath filesep 'lftInfo.mat']);
mat_lft   = full(lftInfo.Mat_lifetime);
mat_stat  = full(lftInfo.Mat_status);

vec_lft   = max(mat_lft,[],2);
mat_stat(mat_stat==0) = nan;
vec_stat = nanmin(mat_stat,[],2);

% if intensity information is used, upload trackInfo and extract intensity data here
if length(restvector)>5
    tmat = load([strrep(lftinfopath, 'LifetimeInfo', 'TrackInfoMatrices') filesep 'trackInfo.mat']);
    tmat = full(tmat.trackInfo);
    imat = tmat(:,4:8:size(tmat,2));
    % maximum intensity
    ivec = nanmax(imat,[],2);
    % intensity vector of all usable positions (>4 frames lft)
    ivec_all = ivec;
    ivec_all(vec_lft<4) = NaN;
    % calculate intensity percentile rank vector
    ivec_rank_all = intensityPercentageRank(ivec_all);
end


% =====================================================================
% DESIRED CONDITIONS
% for all frames in the movie, collect those locations of points that
% fulfill a number of requirements
% specifically status, minimum/maximum lifetime
[rx,ry] = size(restvector);
if min(rx,ry)>1
    nres = size(restvector,2);
else
    nres = 1;
end

for r=1:nres
    
    if nres>1
        crestvector = restvector(r,:);
    else
        crestvector = restvector;
    end
       
    % desired conditions
    dstat       = crestvector(r,1);
    dminfr      = crestvector(r,3);
    minlft      = crestvector(r,4);
    maxlft      = crestvector(r,5);
      
    minlft_fr   = round(minlft/framerate);
    maxlft_fr   = round(maxlft/framerate);
    
    if length(restvector)>5
        min_pint = crestvector(r,6);
        max_pint = crestvector(r,7);
    end
    

    % =====================================================================

    % find positions that fulfill required conditions, as logical matrices
    % correct status
    if isfinite(dstat)
        findpos_stat = ( vec_stat==dstat );
    else
        findpos_stat = isfinite(vec_stat);
    end

    % correct lifetime    
    findpos_lft = ( (vec_lft>dminfr) & (vec_lft>=minlft_fr) & (vec_lft<=maxlft_fr) );
         
    % combine positions
    findpos_all = find(findpos_stat & findpos_lft );
    
    if length(restvector)>5
        
        % OPTION 1 (default) - continued below: 
        % USE PERCENTAGE RANK RELATIVE TO ALL USABLE POSITIONS
              
%             % OPTION 2: (uncomment this mini-paragraph)
%             % USE PERCENTAGE RANK RELATIVE TO THE POSITIONS DEFINED
%             % BY THE LIFETIME CRITERIA ABOVE
%             % limit intensity vector to positions defined by lifetime criteria
%             ivec_curr = nan*ivec_all;
%             ivec_curr(findpos_all) = ivec_all(findpos_all);
%             % calculate intensity percentile rank vector (modified compared to
%             % version from line 78
%             ivec_rank_all = intensityPercentageRank(ivec_curr);
          
        
        % OPTION 1: USE PERCENTAGE RANK RELATIVE TO ALL USABLE POSITIONS
        findpos_intRank = ( (ivec_rank_all>min_pint) & (ivec_rank_all<max_pint) );
        % modify combine positions
        findpos_all = find(findpos_stat & findpos_lft & findpos_intRank);
             
    end
    
    if nres>1
        posvec(r).posvec = findpos_all;
    else
        posvec = findpos_all;
    end
end % of for r-loop
end % of function


%%=========================================================================
%
%                           subfunctions
%
%==========================================================================

function [vec_intPercRank] = intensityPercentageRank(vec_int)

% sorting matrix
% first colums: intensities
imat_sort(:,1) = vec_int;
%second columns: original positions
imat_sort(:,2) = 1:length(vec_int);
% intensity-sorted matrix
imat_sorted1 = sortrows(imat_sort,1);
% assign rank into additional column
imat_sorted1(:,3) = 1:length(vec_int);
% maximum defined position
imaxpos = find(isfinite(imat_sorted1(:,1)), 1, 'last');
% convert absolute rank into percentile rank
imat_sorted1(:,3) = imat_sorted1(:,3)/imaxpos;
% fill non-defined positions from intensity vector with nans
imat_sorted1(imaxpos:length(vec_int),3) = nan;
% sort matrix back to original position (in column 2)
imat_sorted2 = sortrows(imat_sorted1,2);
% extract re-sorted percentage rank
vec_intPercRank = imat_sorted2(:,3);

end


function rankIntensityVect = rankIntensity(intensityVect)

idx = 1:length(intensityVect);
sortMat(:,1) = intensityVect;
sortMat(:,2) = idx;
sortMat = sortrows(sortMat, 1);
% sort places 'NaN' last -> find first valid max
maxIdx = find(isfinite(sortMat(:,1)), 1, 'last');
sortMat(:,3) = idx / maxIdx; % absolute rank -> percentile rank
sortMat(maxIdx+1:end, 3) = NaN;
sortMat = sortrows(sortMat, 2);
rankIntensityVect = sortMat(:,3);
end