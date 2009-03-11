function [mpmRand] = makeRandomMPM(mpm, maskImage, randvar)
% make randomized MPM in the same format as input mpm
% 
% INPUT:    mpm    =   mpm file from which new positions are randomized
%           maskImage (optional) = if maskImage is entered, then randomized
%                   positions are only chosen if they fall into the 'good'
%                   area of the mask
%           randvar (optional) = randomizer variable; if randvar
%                   = 0 (default) positions are drawn from the reshuffled 
%                   of existing x and y positions separately
%                   = 1 positions are completely random in the interval
%                   [0 xmax] amd [0 ymax]
%           
%
% OUTPUT:   mpmRand    =   randomized mpm, in the same format as input mpm,
%                   and with the same number of points in each frame
%
%
%
% Dinah Loerke, modified 07/29/2008

mask = 0;
shuffle = 1;

if nargin>1
    [ix,iy] = size(maskImage);
    if min(ix,iy)>1
        mask = 1;
    end
    if nargin>2
        if randvar==1
            shuffle=0;
        end
    end
end


[sx,sy] = size(mpm);
nf = sy/2;

% initialize randomized mpm
mpmRand = mpm;

% extract x- and y-positions
xmat = mpm(:,1:2:sy);
ymat = mpm(:,2:2:sy);


% used positions in the matrix
umat = 0*xmat;
upos = find( xmat>0  &  ymat>0 );
umat(upos) = 1;

% all collected x- and y-positions
xvec = xmat(upos);
yvec = ymat(upos);
l = length(xvec);

max_x = max(xvec);
max_y = max(yvec);

% if no mask is specified
if mask==0
    
    if shuffle==1
        index_x = randsample(l,l);
        index_y = randsample(l,l);

        xvec_res = xvec(index_x);
        yvec_res = yvec(index_y);
    else
        xvec_res = max_x*rand(length(xvec),1);
        yvec_res = max_y*rand(length(yvec),1);
    end

% else if mask is used        
else
    
    ct = 1;
    while ct<=l
        
        if shuffle==1
            
            index_x = randsample(l,1);
            index_y = randsample(l,1);

            xvec_ct = xvec(index_x);
            yvec_ct = yvec(index_y);
            
        else
            xvec_ct = max_x*rand(1);
            yvec_ct = max_y*rand(1);
            
        end
            
        
        if maskImage(ceil(xvec_ct),ceil(yvec_ct))>0
            
            xvec_res(ct) = xvec_ct;
            yvec_res(ct) = yvec_ct;
            ct = ct+1;
        end
    end % of while
    
end % of if-else

xmat_res = 0*xmat;
ymat_res = 0*ymat;

xmat_res(upos) = xvec_res;
ymat_res(upos) = yvec_res;


mpmRand(:,1:2:sy) = xmat_res;
mpmRand(:,2:2:sy) = ymat_res;


end % of function
    
    