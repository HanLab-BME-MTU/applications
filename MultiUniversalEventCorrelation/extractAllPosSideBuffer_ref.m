function [MPMrefpos] = extractAllPosSideBuffer_ref(MPMallpos, im);
% this function produces an MPM of reference positions for a pre-determined
% set of points of interest, where the have different locations but the
% identical length and appearance/disapperance positions; this ensures that
% the reference positions adequately capture possible bleaching effects or
% global changes during the course of the movie
% 
% SYNOPSIS [MPMrefpos] = extractAllPosSideBuffer_ref(MPMallpos, im);
%
% INPUT     MPMallpos:      MPM with objects of interest
%
%           im:             can be either:
%                           mask of cell outline 
%                           or
%                           vector containing desired image size
%
% OUTPUT:   MPMrefpos:      result MPM
%                           first layer: x and y-coordinates
%                           second layer: time relative to appearance
%                           third layer: time relative to disappearance
%                       
%
% last modified: Dinah Loerke   031/17/2009


% read size of object MPM to define initialization parameters
[nx,ny,nz] = size(MPMallpos);
nf = ny/2;
[ix,iy] = size(im);

mpm_x = MPMallpos(:,1:2:ny,1);
mpm_y = MPMallpos(:,2:2:ny,1);
mpm_t1 = MPMallpos(:,1:2:ny,2);
mpm_t2 = MPMallpos(:,1:2:ny,3);

MPMrefpos = MPMallpos;

% number of reference positions per data point
nref = 8;

% if im is the image size, i.e. a vector [sx,sy], then choose shifted
% object positions as reference positions
if min(ix,iy)==1
    mode = 'shift';
    shiftpos(1,:) = [0 16]; 
    shiftpos(2,:) = [0 -16]; 
    shiftpos(3,:) = [16 0]; 
    shiftpos(4,:) = [-16 0]; 
    shiftpos(5,:) = [12 12];
    shiftpos(6,:) = [12 -12];
    shiftpos(7,:) = [-12 12];
    shiftpos(8,:) = [-12 -12];
    
    imsiz = im;

% if im is a mask, then choose 'scramble' mode and allow only positions
% that fall within the allowed area of the mask    
else
    mode = 'scramble';
    imsiz = size(im);
end

% maximum time values/shifts present in the time layers
stat_t1 = nanmax(abs(mpm_t1),[],2);
stat_t2 = nanmax(abs(mpm_t2),[],2);


%% fill in time shift values
% for each trajectory with lifetime > 4 frames, add the specified number of
% frames to the left and right; the position of the first or last visible
% frame is filled in

gtpos = find( (stat_t1>1) | (stat_t2>1) );


mpm_x_ref = nan*zeros(nref*nx,nf);
mpm_y_ref = nan*zeros(nref*nx,nf);
mpm_t1_ref = nan*zeros(nref*nx,nf);
mpm_t2_ref = nan*zeros(nref*nx,nf);


%% loop over usable trajectories

for i=1:length(gtpos)
    
    k = gtpos(i);
       
    % current traj
    ct_x = mpm_x(k,:); ct_x(ct_x==0)=nan;
    ct_y = mpm_y(k,:); ct_y(ct_y==0)=nan;
    ct_t1 = mpm_t1(k,:);
    ct_t2 = mpm_t2(k,:);
    
    if strcmp(mode,'shift')
        
        creftraj_x = nan*zeros(size(shiftpos,1),nf);
        creftraj_y = nan*zeros(size(shiftpos,1),nf);
        creftraj_t1 = nan*zeros(size(shiftpos,1),nf);
        creftraj_t2 = nan*zeros(size(shiftpos,1),nf);
        
        for s=1:size(shiftpos,1)
            
            shiftx = shiftpos(s,1);
            shifty = shiftpos(s,2);
            
            ct_x_shift = ct_x - shiftx;
            ct_y_shift = ct_y - shifty;
            
            if (nanmin(ct_x_shift)<1) | (nanmin(ct_y_shift)<1)
                continue
            elseif (max(ct_x_shift)>imsiz(2)) | (max(ct_y_shift)>imsiz(1))
                continue
            end
            
            creftraj_x(s,:) = ct_x_shift;
            creftraj_y(s,:) = ct_y_shift;
            creftraj_t1(s,:) = ct_t1;
            creftraj_t2(s,:) = ct_t2;
        end % of for s-loop
        
    elseif strcmp(mode,'scramble')
        
        ct = 1;
        wct = 1;
        while ct<=nref
            xcoor_rand = imsiz(1)*rand(1);
            ycoor_rand = imsiz(2)*rand(1);                 
            if im(ceil(xcoor_rand),ceil(ycoor_rand))>0
                xvec_rand(ct) = xcoor_rand;
                yvec_rand(ct) = ycoor_rand;
                ct = ct+1;
            end
            wct = wct+1;
            if wct>1000
                error('no random points are found');
            end
        end % of while
        
        for s=1:nref
            creftraj_x(s,:) = xvec_rand(s)+0*ct_x;
            creftraj_y(s,:) = yvec_rand(s)+0*ct_y;
            creftraj_t1(s,:) = ct_t1;
            creftraj_t2(s,:) = ct_t2;
        end
    end % of elseif
          
           
    % fill in total matrix
    mpm_x_ref( ((k-1)*nref)+1:(k*nref),: ) = creftraj_x;
    mpm_y_ref( ((k-1)*nref)+1:(k*nref),: ) = creftraj_y;
    mpm_t1_ref( ((k-1)*nref)+1:(k*nref),: ) = creftraj_t1;
    mpm_t2_ref( ((k-1)*nref)+1:(k*nref),: ) = creftraj_t2;
    
    
end % of for i-loop

%fprintf('\n');

MPMrefpos = nan*zeros(nref*nx,2*nf,3);

% merge information into one big MPM file
MPMrefpos(:,1:2:2*nf,1) = mpm_x_ref;
MPMrefpos(:,2:2:2*nf,1) = mpm_y_ref;
MPMrefpos(:,1:2:2*nf,2) = mpm_t1_ref;
MPMrefpos(:,2:2:2*nf,2) = mpm_t1_ref;
MPMrefpos(:,1:2:2*nf,3) = mpm_t2_ref;
MPMrefpos(:,2:2:2*nf,3) = mpm_t2_ref;
    
end % of function

