function [data]=determineLifetimeInfo(data, overwrite);
% determine lifetime info data for all entries in the data structure
%
% SYNOPSIS  [data]=determineLifetimeInfo(data);
%
% INPUT     data    :    experiment structure, which has to contain a field
%                       .source, which is the path to the data location; at
%                       this location, the function reads the trackInfo
%                       from a folder called TrackInfoMatrices
% OUTPUT    data           
% REMARKS 
%
% Dinah Loerke, last modified Feb 2008


% number of entries in data
lens = length(data);

ordir = cd;

owVar = 0;
if nargin>1
    if overwrite==1
        owVar = 1;
    end
end


for i=1:lens
    
    
    fprintf(' movie %02d',i);   
    
    % number of frames for this exp
    lenf = data(i).movieLength;
    
    % current path
    currPath = data(i).source;
    cd(currPath);
    
    lftexist = 0;
    
    if exist('LifetimeInfo')==7
        cd('LifetimeInfo')
        if exist('lftInfo.mat')==2
            lftexist = 1;
        end
        cd(currPath);
    end
    
    % if lifetime Info doesn't already exist or if overwrite is specified
    if (lftexist == 0) | (owVar==1)
        
        % check for TrackInfo
        if exist('TrackInfoMatrices')==7
            cd('TrackInfoMatrices')
            if exist('trackInfo.mat')==2
                % load trackInfo
                loadfile = load('trackInfo.mat');
                if isfield(loadfile,'trackInfo')
                    trackInfo = loadfile.trackInfo;
                elseif isfield(loadfile,'trackInfoMat')
                    trackInfo = loadfile.trackInfoMat;
                else
                    trackInfo = [];
                end
            elseif exist('tracksFinal.mat')==2
                % load trackInfo
                loadfile = load('tracksFinal.mat');
                if isfield(loadfile,'tracksFinal')
                    tracksFinal = loadfile.tracksFinal;
                    [trackInfo] = convStruct2MatNoMS(tracksFinal);
                else
                    trackInfo = [];
                end 
            end
            cd(currPath);
        end   
            
        if ~isempty(trackInfo)
            % calculate lftInfo
            [lftInfo] = calcLifetimeData(trackInfo);
            % make a directory where to save the data, unless it already exists
            if exist('LifetimeInfo')==7
                cd('LifetimeInfo');
            else
                mkdir('LifetimeInfo');
                cd('LifetimeInfo');
            end

            % save data
             save('lftInfo.mat','lftInfo');

        end % of if trackInfo isn't empty
    
    end % if lftInfo doesn't already exist
    
    cd(ordir);
    
    fprintf('\b\b\b\b\b\b\b\b\b');
    
end % of for i

fprintf('\n');

end % of function


%% SUBFUNCTION

function [lftInfo] = calcLifetimeData(trackinfo);
% find lifetimes and status of all tracks in the trackInfo matrix
% SYNOPSIS:
% [lftInfo] = calcLifetimeData(trackinfo);
%
% INPUT:    trackInfo   = trackinfo matrix from Khuloud's tracker
%
% OUTPUT:   lftInfo = structure with the following fields
%
%           Note: the ouput matrices have the same format as the input
%           matrix 'trackinfo'; if the input matrix has the size 
%           (numObj,8*numFrame), then the output matrices have the size 
%           (numObj,numFrame)
%
%           .lftMat      = lifetime matrix (lifetime is the total lifetime 
%                          of the trajectory, regardless of its status)
%           .statMat     = status matrix, the values signify:
%                          1 = valid trajectory (appears and disappears)
%                          2 = partial trajectory (cut off at beginning or
%                          end)
%                          3 = stationary trajectory (present in all frames)                           
%                          4 = "good" gap (shorter than adjacent trajectory pieces)
%                          5 = "bad" gap (longer than adjacent trajectory pieces)
%           .xMat, .yMat  = x,y positions of objects 
%                           NOTE: the positions of gaps are not empty, they
%                           are filled up with the last known position of 
%                           the object before gap closing
%           .disappMat    = matrix with appearance/disappearance events
%                           1=appearance; 0=no change; -1=disappearance
%           .lftVec 	  = vector containing the valid lifetimes (status 1 
%                           with 'good' gaps)
%
% created with MATLAB ver.: 7.1.0.246 (R14) Service Pack 3 on Windows_NT
%
% created by: dloerke
% last modified: 22-Feb-2008
%


%% =======================================================================
%%======        Step 0: initialize and set default values
%% =======================================================================

[stx,sty] = size(trackinfo);
tsiz = round(sty/8);

% statusMat matrix has number of columns corresponding to number of frames,
% and four layers: one for x coordinate, one for y-coordintae, one for
% status, and one for measured lifetime

lftMat  = zeros(stx,tsiz);
statMat = zeros(stx,tsiz);
xMat    = zeros(stx,tsiz);
yMat    = zeros(stx,tsiz);
disappMat = zeros(stx,tsiz);
lftVec  = nan*zeros(stx,1);

if (issparse(trackinfo))
    use_trackinfo = full(trackinfo);
else 
    use_trackinfo = trackinfo;
end



%% =======================================================================
%%======                    Step 1: find lifetimes
%% =======================================================================


% loop over rows (trajectories)

for i=1:stx
    
    % find defined x coordinates
    xcoords = use_trackinfo(i,1:8:sty);
    xcoords_valpos = find(xcoords>0);
    % find defined y coordinates
    ycoords = use_trackinfo(i,2:8:sty);
    ycoords_valpos = find(ycoords>0);
    
    if (length(xcoords_valpos)>0)
    
        pmin = min(xcoords_valpos);
        pmax = max(xcoords_valpos);
        traj_len = pmax-pmin+1;
        traj_np = length(xcoords_valpos);
    
        vec_x = xcoords; 
        vec_y =  ycoords; 
        vec_stat    =  zeros(tsiz,1); 
        vec_lft     =  zeros(tsiz,1);
        vec_disapp  = zeros(tsiz,1);
    
        %for the initial appearance, the value for vec_disapp is set to 1, for
        %the final disappearance, it's set to -1
        if pmin>1, vec_disapp(pmin) = 1; end
        if pmax<tsiz, vec_disapp(pmax) = -1; end 
    
        % if the track length is the same as the movie length, we have a
        % stationary point
        if traj_len == tsiz
            vec_stat(xcoords_valpos) = 3;
            vec_lft(xcoords_valpos) = tsiz;
        else
            % if the track starts after frame 1 and end before the last frame,
            % we have a valid full trajectory (i.e. one that appears and
            % disappears)
            if (pmin>1) & (pmax<tsiz)
                vec_stat(xcoords_valpos) = 1;
                vec_lft(xcoords_valpos) = traj_len;
                % else we have a partial trajectory (where the trajectory is cut 
                % off by the time series beginning or ending)
            else
                vec_stat(xcoords_valpos) = 2;
                vec_lft(xcoords_valpos) = traj_len;
            end
        end % of if trajectory is stationary, partial, or valid
    
        % regardless of the full trajectory status, now check for closed gaps
    
        % if gaps exist at all
        if (traj_len > traj_np)   
    
        % We want to be able to exclude those tracks which consist of a single
        % detection connected by a long closed gap to another single detection
        % instance, because these kinds of trajectories constitute random 
        % detection noise connected by the gap closer, rather than being 
        % sustained trajectories.
        % What kind of criterion do we apply?
        % In this simple version, a "bad" gap is a gap>2 flanked by single
        % detection events on either side - a "good" gap is everything else
    
            traj = xcoords;
            % gap status is 1 for existence of trajectory, 0 for empty
            gap_status = (xcoords>0);
    
           
        % We're looking for the position where the gap is located, which is 
        % determined from where diff(gap_status) is -1; however, nonzero
        % values for diff(gap_status) also occur for appearing and
        % disappearing of the trajectory. In the subsequent step, the
        % differential of the vector with the nonzero positions yields the
        % length of the gap and the adjacent trajectory pieces.
        % Since we're dealing with both partial and fully appearing/
        % disappearing trajectories, the vector is buffered with 1 and tsiz
        % to yield a trajectory which appears (if partial at frame one), 
        % has a gap, reappears, and then disappears (if partial at frame tsiz)
    
            diffgap = diff(gap_status);
            gap_posFirst = min(find(diffgap==-1));
            gap_posLast = max(find(diffgap==-1));
            gap_posAll = find(abs(diffgap)>0);
            % if the track doesn't appear, the values are buffered with 1
            if (min(gap_posAll)==gap_posFirst)
                gap_posAll = [1 gap_posAll];
            end
            % if the track doesn't disappear, the values are buffered with tsiz
            if (max(gap_posAll)>gap_posLast)
                gap_posAll = [gap_posAll tsiz];
            end

            gap_piecelengths = diff(gap_posAll);
            gp_len = length(gap_piecelengths);
        
            % loop over the (possibly multiple) gaps in the trajectory
            for p = 2:2:gp_len
                 % the gap is good if either the number of detected points is >2, or
                % if the gap between the points is not larger than 1 empty frame
                if ( traj_np > 2 ) | ( traj_len-traj_np < 2 )
                    vec_stat(gap_posAll(p)+1:gap_posAll(p+1)) = 4;
                    vec_lft(gap_posAll(p)+1:gap_posAll(p+1)) = gap_piecelengths(p);
                    % since x,y values in gaps are zero, fill up with the last
                    % known xy position
                    vec_x(gap_posAll(p)+1:gap_posAll(p+1)) = vec_x(gap_posAll(p));
                    vec_y(gap_posAll(p)+1:gap_posAll(p+1)) = vec_y(gap_posAll(p));
                % else the gap is "bad"
                else
                    vec_stat(gap_posAll(p)+1:gap_posAll(p+1)) = 5;
                    vec_lft(gap_posAll(p)+1:gap_posAll(p+1)) = gap_piecelengths(p);
                    % since x,y values in gaps are zero, fill up with the last
                    % known xy position
                    vec_x(gap_posAll(p)+1:gap_posAll(p+1)) = vec_x(gap_posAll(p));
                    vec_y(gap_posAll(p)+1:gap_posAll(p+1)) = vec_y(gap_posAll(p));
                end % of if
            end % of for p
        
        end % of if gaps exist
    
        xMat(i,:) = vec_x;
        yMat(i,:) = vec_y;
        statMat(i,:)    = vec_stat;
        lftMat(i,:)     = vec_lft;
        disappMat(i,:)  = vec_disapp;
    
            
    end % of if points exists in this row
    
    % if the status of the trajectory is ==1 and the value of the gaps is
    % ==4, then this trajectory is counted in the lifetimes vector
    if ( (min(nonzeros(vec_stat))==1) & (max(vec_stat)<5) )
        lftVec(i)=max(vec_lft);
    end
    
end % of for loop over rows



lftInfo.Mat_xcoord    = sparse(xMat);
lftInfo.Mat_ycoord    = sparse(yMat);
lftInfo.Mat_status = sparse(statMat);
lftInfo.Mat_lifetime  = sparse(lftMat);
lftInfo.Mat_disapp = sparse(disappMat);




end % of function

