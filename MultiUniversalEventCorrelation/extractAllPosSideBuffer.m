function [MPMallpos] = extractAllPosSideBuffer(lftinfopath, tbuffer);
% this function extracts an MPM of points of interest from trackInfo using
% a side buffer on both sides of the trajectory, i.e. with a specified
% numver of points before appearance and after disappearance of the object
% 
% SYNOPSIS  [MPMallpos] = extractAllPosSideBuffer(lftinfopath, tbuffer);
%
% INPUT     lftinfopath:    the path pointing to the location of the
%                           lftInfo data
%
%           tbuffer:        time buffer on both sides of the trajectory (in
%                           frames)
%
% OUTPUT:   MPMallpos:      result MPM
%                           first layer: x and y-coordinates
%                           second layer: time relative to appearance
%                           third layer: time relative to disappearance
%                       
%
% last modified: Dinah Loerke   01/30/2009




%% extract lft matrices
cd(lftinfopath);
loadfile = open('lftInfo.mat');
lmat = loadfile.lftInfo;

mat_lft   = full(lmat.Mat_lifetime);
mat_stat  = full(lmat.Mat_status);
mat_x     = full(lmat.Mat_xcoord);
mat_y     = full(lmat.Mat_ycoord);
mat_da    = full(lmat.Mat_disapp);
vec_lft   = max(mat_lft,[],2);

% frame number
[nx,nf] = size(mat_lft);

% initiate mpm sub-files containing x,y corrdinates and relative time shifts
mpm_x = mat_x;
mpm_y = mat_y;
mpm_t1 = nan*mpm_x;
mpm_t2 = nan*mpm_x;



%% fill in time shift values
% for each trajectory with lifetime > 4 frames, add the specified number of
% frames to the left and right; the position of the first or last visible
% frame is filled in

gtpos = find(vec_lft>4);

for i=1:length(gtpos)
    
    k = gtpos(i);
     
    % current traj
    ct_x = mat_x(k,:);
    ct_y = mat_y(k,:);
    ct_da = mat_da(k,:);
    
    % appearance and disappearance positions
    pos_a = find(ct_da==1);
    pos_da = find(ct_da==-1);
    
    % set appropriate buffered start and end positions (cannot fall outside
    % the range of the movie)
    if ~isempty(pos_a)
        startpos = max(1,pos_a-tbuffer);
    else
        startpos = 1;
    end
    
    if ~isempty(pos_da)
        endpos = min(nf,pos_da+tbuffer);
    else
        endpos = nf;
    end
    
    % if trajectory is seen to appear (i.e. has pre-appearance buffer),
    % define t1 as time from appearance
    if ~isempty(pos_a)        
        % fill in buffer zone positions
        mpm_x(k,startpos:pos_a-1) =  mpm_x(k,pos_a);
        mpm_y(k,startpos:pos_a-1) =  mpm_y(k,pos_a);
        % fill in pre-appearance time
        mpm_t1(k,startpos:pos_a-1) = -pos_a + [startpos:pos_a-1];
        % fill in trajectory times
        mpm_t1(k,pos_a:endpos) = [pos_a:endpos] - pos_a; 
    end
    
    % is trajectory is seen to disappear (i.e. has post-disappearance
    % buffer), define t2 as time from disappearance
    if ~isempty(pos_da)
        % fill in buffer zone
        mpm_x(k,pos_da+1:endpos) =  mpm_x(k,pos_da);
        mpm_y(k,pos_da+1:endpos) =  mpm_y(k,pos_da);
        % fill in post-disappearance time
        mpm_t2(k,pos_da+1:endpos) = [pos_da+1:endpos] - pos_da;
        % fill in trajectory times
        mpm_t2(k,startpos:pos_da) = -pos_da + [startpos:pos_da];
    end 
    
         
end % of for i-loop

%fprintf('\n');

MPMallpos = nan*zeros(nx,2*nf,3);

% merge information into one big MPM file
MPMallpos(:,1:2:2*nf,1) = mpm_x;
MPMallpos(:,2:2:2*nf,1) = mpm_y;
MPMallpos(:,1:2:2*nf,2) = mpm_t1;
MPMallpos(:,2:2:2*nf,2) = mpm_t1;
MPMallpos(:,1:2:2*nf,3) = mpm_t2;
MPMallpos(:,2:2:2*nf,3) = mpm_t2;
    
end % of function

