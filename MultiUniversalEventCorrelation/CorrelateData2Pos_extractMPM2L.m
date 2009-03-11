function [MPMglobal] = CorrelateData2Pos_extractMPM2L(lftinfopath, restvector, tvec, fr, ttype);
% this function extracts an MPM of points of interest from trackInfo using
% specified restrictions, and time delays in the third dimension
% 
% SYNOPSIS [MPMglobal] = CorrelateData2Pos_extractMPM(trackinfo,restvector, tvec, ttype);
%
% INPUT     trackinfo       = track info file
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
%           tvec            = time vector with time shifts
%                               e.g. [-15:1:15]
%
%           fr              = framerate; default 1
%
%           ttype           = variable to indicate the strategy for
%                             choosing shifted positions (optional)
%                             0 (default): shifted positions are the same
%                             x,y-positions as the reference position
%                             1: shifted positions are the positions of the
%                             respective trajectory at this time point (if
%                             available)
%
% OUTPUT:   MPMglobal       = result MPM
%
% last modified: Dinah Loerke   08/28/2008


%% determine time mode
tmod = 0;
if nargin>2
    if ~isempty(tvec)
        tmod = 1;
    end
end

%% determine framerate
framerate = 1;
if nargin>3
    if ~isempty(fr)
        framerate = fr;
    end
end

%% determine position fill mode
fillmod = 0;
if nargin>4
    if ttype ==1
        fillmod = 1;
    end
end

%% calculate lft matrices
cd(lftinfopath);
loadfile = open('lftInfo.mat');
lmat = loadfile.lftInfo;

CMlft   = full(lmat.Mat_lifetime);
CMstat  = full(lmat.Mat_status);
CMx     = full(lmat.Mat_xcoord);
CMy     = full(lmat.Mat_ycoord);
CMda    = full(lmat.Mat_disapp);
CVlft   = max(CMlft,[],2);


clusterStatVec = nan * CVlft;
if exist('clusterProperties.mat')==2
    loadfile = open('clusterProperties.mat');
    clusterProps = loadfile.clusterProperties;
    clusterStatVec = clusterProps.clusterStatus;
end



%% extract current conditions
% EXISTING DATA CONDITIONS
% status matrix
mat_stat = full(CMstat);
% frame number
[nx,nf] = size(mat_stat);
% % status vector
% vec_stat = zeros(nx,1);
% for k=1:nx, vec_stat(k) = min(nonzeros(mat_stat(k,:))); end

% lifetime matrix
mat_lft = full(CMlft);
% % lifetime vector
% vec_lft = zeros(nx,1);
% for k=1:nx, vec_lft(k) = max(mat_lft(k,:)); end
    
% x-coordinate matrix
mat_x = full(CMx);
% y-coordinate matrix
mat_y = full(CMy);
% disapp status matrix
mat_da = full(CMda);


% frame number
[nx,nf] = size(mat_lft);
    
%     % intensity - rescale values to range between 0 and 1, where the value
%     % corresponds to the percentile of all available intensities
%     svec_int = [];
%     svec_int(:,1) = currLI.Vec_int; le = length(currLI.Vec_int);
%     svec_int(:,2) = 1:le;
%     svec_int = sortrows(svec_int,1);
%     svec_int(:,3) = [1:le]/le;
%     svec_int = sortrows(svec_int,2);
%     % mat_int = []; mat_int = repmat(svec_int(:,3),1,nf);
%     vec_int = svec_int(:,3);
%     
%     % motility - rescale values to range between 0 and 1
%     svec_mot = [];
%     svec_mot(:,1) = currLI.Vec_displ; le = length(currLI.Vec_displ);
%     svec_mot(:,2) = 1:le;
%     svec_mot = sortrows(svec_mot,1);
%     svec_mot(:,3) = [1:le]/le;
%     svec_mot = sortrows(svec_mot,2);
%     % mat_mot = []; mat_mot = repmat(svec_mot(:,3),1,nf);
%     vec_mot = svec_mot(:,3);
    


%% =====================================================================
% DESIRED CONDITIONS
% for all frames in the movie, collect those locations of points that
% fulfill a number of requirements
% e.g. status, minimum lifetime, maximum motility

if min(size(restvector))==1

    % desired minimum lifetime in seconds
    dstat       = restvector(1);
    dapp        = restvector(2);
    dminfr      = restvector(3);
    minlft      = restvector(4);
    maxlft      = restvector(5);
    
    minlft_fr   = round(minlft/framerate);
    maxlft_fr   = round(maxlft/framerate);

    % beginning and end frame
    amin        = 0; 
    amax        = nf;
    if length(restvector)>5
        if ~isnan(restvector(6))
            amin = restvector(6);
        end
        if ~isnan(restvector(7))
            amax = restvector(7); 
        end
    end
    
    % cluster status
    clus_stat = nan;
    if length(restvector)>7
        clus_stat = restvector(8); 
        if (isfinite(clus_stat)) & (sum(isfinite(clusterStatVec))==0)
            disp('warning: no specified cluster properties available');
        end
    end
    
    %% =====================================================================

    % find positions that fulfill required conditions, as logical matrices
    % correct status
    if isfinite(dstat)
        findpos_stat = ( mat_stat==dstat );
    else
        findpos_stat = (mat_x>0);
    end
    % correct disapperance status
    if isfinite(dapp)
        findpos_dapp   = ( mat_da==dapp );
    else
        findpos_dapp   = (mat_x>0);
    end
    % correct lifetime    
    findpos_lft = ( (mat_lft>dminfr) & (mat_lft>minlft_fr) & (mat_lft<maxlft_fr) );
    
    if isfinite(clus_stat)
        findpos_clusstat_vec = ( clusterStatVec==clus_stat );
        findpos_clusstat = repmat(findpos_clusstat_vec,[1 nf]);
    else
        findpos_clusstat = (mat_x>0);
    end
    
    % combine positions
    findpos_all = find(findpos_stat & findpos_dapp & findpos_lft & findpos_clusstat);
    % convert to subscript
    [findpos_p,findpos_f] = ind2sub(size(mat_da),findpos_all);

    % if necessary, kick out the results that are not in the desired frame
    % range
    fpos_framerange = find( (findpos_f>=amin) & (findpos_f<=amax) );
    findpos_f = findpos_f(fpos_framerange);
    findpos_p = findpos_p(fpos_framerange);
       
    findpos_all = sub2ind(size(mat_da),findpos_p,findpos_f);


    % read out x and y coordinates in these positions
    findx = nan*mat_x;
    findx(findpos_all) = mat_x(findpos_all);
    findy = nan*mat_y;
    findy(findpos_all) = mat_y(findpos_all);
    
    MPMreferenceZero = nan*mat_x;
    MPMreferenceZero(:,1:2:2*nf) = findx;
    MPMreferenceZero(:,2:2:2*nf) = findy;
    
else
    findx = restvector(:,1:2:size(restvector,2));
    findy = restvector(:,2:2:size(restvector,2));
    findpos_all = find((findx>0) & (findy>0));
    [findpos_p,findpos_f] = ind2sub(size(findx),findpos_all);
    MPMreferenceZero = restvector;
end

   


% if there's only one zero timeshift, the write to result and stop function
% here
if tmod==0
    MPMglobal = MPMreferenceZero;
    return
end

% initialize global results
MPMglobal = MPMreferenceZero;
MPMtvalues = MPMreferenceZero;


%% next step: fill in time shift values

fprintf('time shift #');

for t=1:length(tvec)
    
    fprintf('%02d',t);
    % time shift
    tau = round(tvec(t));
    %
    MPMshift = MPMreferenceZero;
    
    
    % for each non-zero time-shift, we fill out the position of each found
    % value in MPMreferenceZero, but shifted by tau frames
    % there are two possibilities for the actual [x,y] values that are
    % filled into this position:
    % 1. the same position as in MPMreferenceZero is filled in (default)
    % which means that we're tracking whatever parameter (intensity,
    % density, etc.) in the same event location, before and after the event
    % 2. if available, the position is filled in with the current values of
    % the trajectory, which means that we're following not the same spot,
    % but the same moving trajectory before and after
    
    if tau~=0
        
        % initialized matrices
        findxShift = nan*findx;
        findyShift = nan*findy;
        
        % shifted positions
        fillpos_p = findpos_p;
        fillpos_f = findpos_f+tau;
        
        badpos = find( (fillpos_f<1) | (fillpos_f>nf) );
        
        % shifted positions that lie outside the mpm size are set to nan
        fillpos_p(badpos) = nan;
        fillpos_f(badpos) = nan;
    
        % loop over all positions
        for i=1:length(fillpos_p)
            % if shifted position is inside size of mpm file, i.e. not-nan
            if ~isnan(fillpos_p(i))
                % fill with event position i
                findxShift(fillpos_p(i),fillpos_f(i)) = findx(findpos_p(i),findpos_f(i));
                findyShift(fillpos_p(i),fillpos_f(i)) = findy(findpos_p(i),findpos_f(i));
                
                % if fillmod is not zero, then overwrite with trajectory
                % position in this location, provided one is available
                if (fillmod>0) 
                    % status at this position
                    cstat = mat_stat(fillpos_p(i),fillpos_f(i));
                    % if there's a detected object, use this position
                    if (cstat>0) & (cstat<4) 
                        findxShift(fillpos_p(i),fillpos_f(i)) = mat_x(fillpos_p(i),fillpos_f(i));
                        findyShift(fillpos_p(i),fillpos_f(i)) = mat_y(fillpos_p(i),fillpos_f(i));
                    % else find the nearest previous/subsequent detected position and use it
                    else
                        % try previous positions first: vector of previous
                        % status values
                        svec = mat_stat(fillpos_p(i),1:fillpos_f(i));
                        % rightmost (highest) position where an object was
                        % present (no gap, notr empty)
                        mpos = max(find((svec<4)&(svec>0)));
                        if isempty(mpos)
                            svec = mat_stat(fillpos_p(i),:);
                            fpos = find( (svec<4)&(svec>0) );
                            mpos = fillpos_f(i) + min(fpos-fillpos_f(i));
                        end
                        findxShift(fillpos_p(i),fillpos_f(i)) = mat_x(fillpos_p(i),mpos);
                        findyShift(fillpos_p(i),fillpos_f(i)) = mat_y(fillpos_p(i),mpos);
                    end
                end % of if fillmod>0
            end % of if position not-nan
        end % of for i-loop
        
        MPMshift(:,1:2:2*nf) = findxShift;
        MPMshift(:,2:2:2*nf) = findyShift;
        
        fpos = isfinite(MPMshift);
        % defined positions for this time shift
        MPMglobal(fpos) = MPMshift(fpos);
        
        
    end % of if tau>0
    
    % for all layers (regardless of whether or not the shift is zero) set
    % the time shift value to t
    MPMtvalues(isfinite(MPMshift)) = t;    
        
    if t<length(tvec)
        fprintf('\b\b');
    end
    
end % of for t-loop
fprintf('\n');

MPMglobal(:,:,2) = MPMtvalues;

end % of function

