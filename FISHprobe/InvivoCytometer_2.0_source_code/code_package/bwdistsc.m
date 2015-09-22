function D=bwdistsc(bw,aspect)
% D=BWDISTSC(BW,ASPECT)
% BWDISTSC computes Euclidean distance transform of the binary 3D image 
% BW. For each pixel in BW, the distance transform assignes a number
% that is the distance between that pixel and the nearest nonzero pixel 
% of BW. BW may be a single 2D image, 3D array or a cell array of 
% 2D slices. ASPECT is 3-component vector defining aspect ratio in 
% the dataset BW. If ASPECT is not specified, isotropic aspect 
% ratio [1 1 1] is assumed.
%
% BWDISTSC uses fast optimized scanning algorithm and cell-arrays to 
% represent internal data, so that it is less demanding to physical 
% memory. In many cases BWDISTSC is actually faster than MATLAB's 
% optimized kd-tree algorithm used for Euclidean distance 
% transform in 3D. 
%
% BWDISTSC tries to use MATLAB bwdist for 2D scans if possible, which 
% is significantly faster. Otherwise BWDISTSC uses internal algorithm 
% to perform 2D scans.
%
%     Yuriy Mishchenko  JFRC HHMI Chklovskii Lab  JUL 2007

% This code is free for use or modifications, just please give credit 
% where appropriate. And if you modify code or fix bugs, please drop 
% me a message at gmyuriy@hotmail.com.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scan algorithms below use following Lema:                     %
% LEMA: let F(X,z) be lower envelope of a family of parabola:   %
% F(X,z)=min_{i} [G_i(X)+(z-k_i)^2];                            %
% and let H_k(X,z)=A(X)+(z-k)^2 be a parabola.                  %
% Then for H_k(X,z)==F(X,z) at each X there exist at most       %
% two solutions k1<k2 such that H_k12(X,z)=F(X,z), and          %
% H_k(X,z)<F(X,z) is restricted to at most k1<k2.               %
% Here X is any-dimensional coordinate.                         %
%                                                               %
% Thus, simply scan away from any z such that H_k(X,z)<F(X,z)   %
% in either direction as long as H_k(X,z)<F(X,z) and update     %
% F(X,z). Note that need to properly choose starting point;     %
% starting point is any z such that H_k(X,z)<F(X,z); z==k is    %
% usually, but not always the starting point!!!                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse inputs
if(nargin<2 || isempty(aspect)) aspect=[1 1 1]; end

% determine geometry of data
if(iscell(bw)) shape=[size(bw{1}),length(bw)]; else shape=size(bw); end

% fix to handle 1D & 2D images
if(length(shape)<3) shape(length(shape)+1:3)=1; end
if(length(aspect)<3) aspect(length(aspect)+1:3)=1; end

% allocate space
D=cell(1,shape(3)); for k=1:shape(3) D{k}=zeros(shape(1:2)); end

%%%%%%%%%%%%% scan along XY %%%%%%%%%%%%%%%%
for k=1:shape(3)    
    if(iscell(bw)) bwXY=bw{k}; else bwXY=bw(:,:,k); end
        
    % initialize arrays
    DXY=zeros(shape(1:2));
    D1=zeros(shape(1:2));
    DK=zeros(shape(1:2));

    % if can, use 2D bwdist from image processing toolbox    
    if(exist('bwdist') && aspect(1)==aspect(2))
        D1=aspect(1)^2*bwdist(bwXY).^2;
    else    % if not, use full XY-scan
        %%%%%%%%%%%%%%% X-SCAN %%%%%%%%%%%%%%%        
        % reference nearest bwXY "on"-pixel in x direction downward:
        
        %  scan bottow-up, copy x-reference from previous row unless 
        %  there is bwXY "on"-pixel in that point in current row
        xlower=repmat(Inf,shape(1:2)); 
        
        xlower(1,find(bwXY(1,:)))=1;    % fill in first row
        for i=2:shape(1)
            xlower(i,:)=xlower(i-1,:);  % copy previous row
            xlower(i,find(bwXY(i,:)))=i;% unless there is pixel
        end
        
        % reference nearest bwXY "on"-pixel in x direction upward:
        xupper=repmat(Inf,shape(1:2));
        
        xupper(end,find(bwXY(end,:)))=shape(1);
        for i=shape(1)-1:-1:1
            xupper(i,:)=xupper(i+1,:);
            xupper(i,find(bwXY(i,:)))=i;
        end
                
        % find points for which distance needs to be updated
        idx=find(~bwXY); [x,y]=ind2sub(shape(1:2),idx);
        
        % set distance as the shortest to upward or to downward
        DXY(idx)=aspect(1)^2*min((x-xlower(idx)).^2,(x-xupper(idx)).^2);
        
        %%%%%%%%%%%%%%% Y-SCAN %%%%%%%%%%%%%%%
        % this will be the envelop 
        % envelop is initialized at Inf to ensure single scan direction, 
        % otherwise may end up in infinite loop when trying to find 
        % starting point
        D1=repmat(Inf,shape(1:2));
        % these will be the references to parabolas defining the envelop
        DK=repmat(Inf,shape(1:2));
        % starting points
        i0=zeros(shape(1),1);
        % convenience x-coords array 
        x=(1:shape(1))'; 
        
        for i=1:shape(2)
            % need to select starting point for each X:
            % * starting point should be below current envelop
            % * i0==i is not necessarily a starting point
            % * there is at most one starting point
            % * there may be no starting point
            
            % i0 is the starting points for each X: i0(X) is the first 
            % y-index such that parabola from line i is below the envelop
            
            % first guess is the current y-line
            i0(:)=i;
            
            % some auxiliary datasets
            d0=DXY(:,i);

            % L0 indicates for which X starting point had been fixed
            L0=isinf(d0) | (d0==0);
            
            while(~isempty(find(~L0,1)))
                % reference starting points in DXY
                idx=sub2ind(shape(1:2),x(~L0),i0(~L0));
                
                % reduce out trivial points (DXY==0)
                L=(DXY(idx)==0);
                L0(~L0)=L;
                idx=idx(~L);
                
                if(isempty(idx)) continue; end
                
                % these are current best parabolas for starting points
                ik=DK(idx);
                
                % these are new values from parabola from line #i
                dtmp=d0(~L0)+aspect(2)^2*(i0(~L0)-i).^2;
                
                % these starting points are OK - below the envelop
                L=D1(idx)>dtmp; D1(idx(L))=dtmp(L);
                
                % points which are still above the envelop but ik==i0,
                % will not get any better, so fix them as well
                L=L | (ik==i0(~L0));
                                
                % all other points are not OK, need new starting point:
                % starting point should be at least below parabola 
                % beating us at current choice of i0
                
                % solve quadratic equation to find where this happens
                ik=(ik-i); 
                di=(D1(idx(~L))-dtmp(~L))./ik(~L)/2/aspect(2)^2;

                % should select next highest index to the equality
                di=fix(di)+sign(di);
                
                % the new starting points
                idx=find(~L0); 
                i0(idx(~L))=i0(idx(~L))+di;

                % update L0 to indicate which points we've fixed
                L0(~L0)=L; L0(idx(~L))=(di==0);
                
                % points that went out of boundaries can't get better; 
                % fix them as well
                idx=idx(~L);
                idx=idx((i0(idx)<1) | (i0(idx)>shape(2))); 
                i0(idx)=i;
                L0(idx)=1;                
            end

            % reduce out trivial points DXY(idx)<DXY(:,i)
            idx=sub2ind(shape(1:2),x,i0);
            L=(DXY(idx)>0) | (i0==i);
            idx=idx(L);

            % these will keep track along which X should 
            % keep updating distances            
            map_lower=L;
            map_upper=L;            
            idx_lower=idx;
            idx_upper=idx;
            
            % set trivial pixels D==0 in line #i:
            % this has to be done b/s we manually discarded them from L0
            D1(d0==0,i)=0;

            % scan from starting points for each X,i0 in increments of 1
            di=0;       % distance from current y-line
            eols=2;     % end-of-line-scan flag
            totlen=prod(shape(1:2));
            while(eols)
                eols=2;
                di=di+1;
                
                % select X which can be updated for di<0;
                % i.e. X which had been below envelop all way till now
                if(~isempty(idx_lower))
                    % shift y by -1
                    idx_lower=idx_lower-shape(1);
                    
                    % prevent index dropping below 1st
                    L=(idx_lower>=1);
                    map_lower(map_lower)=L;
                    idx_lower=idx_lower(L);
                    
                    if(~isempty(idx_lower))
                        dtmp=d0(map_lower)+...
                            aspect(2)^2*(i0(map_lower)-di-i).^2;
                    
                        % these pixels are to be updated with i0-di
                        L=D1(idx_lower)>dtmp & DXY(idx_lower)>0;
                        map_lower(map_lower)=L;
                        idx_lower=idx_lower(L);
                        D1(idx_lower)=dtmp(L);
                        DK(idx_lower)=i;
                    end
                else    % if this is empty, get ready to quit
                    eols=eols-1;
                end

                % select X which can be updated for di>0;
                % i.e. X which had been below envelop all way till now
                if(~isempty(idx_upper))
                    % shift y by +1
                    idx_upper=idx_upper+shape(1);
                    
                    % prevent index from going over array limits
                    L=(idx_upper<=totlen);
                    map_upper(map_upper)=L;
                    idx_upper=idx_upper(L);
                    
                    if(~isempty(idx_upper))                        
                        dtmp=d0(map_upper)+...
                            aspect(2)^2*(i0(map_upper)+di-i).^2;
                    
                        % check which pixels are to be updated with i0+di
                        L=D1(idx_upper)>dtmp & DXY(idx_upper)>0;
                        map_upper(map_upper)=L;
                        idx_upper=idx_upper(L);
                        D1(idx_upper)=dtmp(L);
                        DK(idx_upper)=i;
                    end
                else    % if this is empty, get ready to quit
                    eols=eols-1;
                end  
            end
        end
    end
    D{k}=D1; 
end

%%%%%%%%%%%%% scan along Z %%%%%%%%%%%%%%%%
% this will be the envelop: 
% envelop has to be initialized at Inf to ensure single direction of scan, 
% otherwise may end up in infinite loop when trying to find starting point
D1=cell(size(D));
for k=1:shape(3) D1{k}=repmat(Inf,shape(1:2)); end
% these will be the references to parabolas defining the envelop
DK=cell(size(D));
for k=1:shape(3) DK{k}=repmat(Inf,shape(1:2)); end

% start building the envelope 
for k=1:shape(3)
    % need to select starting point for each X:
    % * starting point should be below current envelop
    % * k0==k is not necessarily a starting point
    % * there may be no starting point
    
    % k0 is the starting points for each XY: k0(XY) is the first
    % z-index such that parabola from line k is below the envelop

    % initial starting point guess is current slice
    k0=repmat(k,shape(1:2));
    
    % L0 indicates which starting points had been fixed
    L0=isinf(D{k}) | (D{k}==0);
    idxtot=find(~L0);
    
    while(~isempty(idxtot))
        % because of using cells need to explicitly scan in Z
        % to avoid repetitious searches in k0, parse first
        ss=getregions(k0(idxtot));
        sslen=length(ss);
        
        for kk=1:sslen
            % these are starting points @kk which had not been set
            idx=idxtot(ss(kk).PixelIdxList);
            
            % reduce out trivial points (D==0)
            if(kk~=k)
                L=(D{kk}(idx)==0);
                L0(idx)=L;
                idx=idx(~L);
            end

            if(isempty(idx)) continue; end
            
            % these are currently best parabolas for slice kk
            ik=DK{kk}(idx);            
            
            % these are new values for slice kk from parabola from k
            dtmp=D{k}(idx)+aspect(3)^2*(kk-k)^2;
            
            % these points are OK - below current envelop
            L=D1{kk}(idx)>dtmp; D1{kk}(idx(L))=dtmp(L);            
            
            % these points are not OK, but since ik==k0
            % can't get any better
            L=L | (ik==kk);            
                
            % all other points are not OK, need new starting point:
            % starting point should be at least below parabola
            % beating us at current choice of k0
            
            % solve quadratic equation to find where this happens
            ik=(ik-k);
            dk=(D1{kk}(idx(~L))-dtmp(~L))./ik(~L)/2/aspect(3)^2;
            dk=fix(dk)+sign(dk);
            k0(idx(~L))=k0(idx(~L))+dk;
    
            % update starting points that had been set
            L0(idx)=L;
            L0(idx(~L))=(dk==0);
    
            % points that went out of boundaries can't get better
            idx=idx(~L);
            idx=idx((k0(idx)<1) | (k0(idx)>shape(3)));
            L0(idx)=1;
            k0(idx)=k;
        end

        idxtot=find(~L0);
    end
    
    % map_lower/map_upper keeps track of which pixels can be yet updated
    % with new distance, i.e. all such XY that had been below envelop for
    % all dk up to now for dk<0/dk>0 respectively
    map_lower=true(shape(1:2));
    map_upper=true(shape(1:2));

    % parse different values in k0 to avoid repetitious searching below
    ss=getregions(k0); 
    sslen=length(ss);

    % reduce out trivially faulty starting points
    for kk=1:sslen
        if(kk==k) continue; end
        
        idx=ss(kk).PixelIdxList;
        
        L=D{kk}(idx)>D{k}(idx);
        map_lower(idx)=L;
        map_upper(idx)=L;
    end
    
    % these are maintained to keep fast track of whether maps are empty
    idx_lower=find(map_lower);
    idx_upper=find(map_upper);
    
    % set trivial pixels D==0 in slice k:
    % this has to be done b/s we manually discarded them from L0
    D1{k}(D{k}==0)=0;        

    % scan away from starting points in increments of 1
    dk=0;       % distance from current xy-slice
    eols=2;     % end-of-scan flag
    while(eols)
        eols=2;
        dk=dk+1;

        if(~isempty(idx_lower))
            % prevent index from going over the boundaries
            L=(k0(map_lower)-dk>=1);
            map_lower(map_lower)=L;
            % need to explicitly scan in Z because of using cell-arrays
            for kk=1:sslen-dk
                % get all XY such that k0-dk==kk
                idx=ss(kk+dk).PixelIdxList; 
                L=map_lower(idx);
                idx=idx(L);                

                if(~isempty(idx))
                    dtmp=D{k}(idx)+aspect(3)^2*(kk-k)^2;
                    
                    % these pixels are to be updated with k0-dk
                    L=D1{kk}(idx)>dtmp & D{kk}(idx)>0;
                    map_lower(idx)=L;
                    D1{kk}(idx(L))=dtmp(L);
                    
                    % ridiculously, but this is faster than
                    % direct assignment
                    dtmp=idx(L);
                    dtmp(:)=k;
                    DK{kk}(idx(L))=k;
                end
            end
            idx_lower=idx_lower(map_lower(idx_lower));
        else
            eols=eols-1;
        end

        if(~isempty(idx_upper))
            % prevent index from going over the boundaries            
            L=(k0(map_upper)+dk<=shape(3));
            map_upper(map_upper)=L;
            % need to explicitly scan in Z because of using cell-arrays
            for kk=dk+1:min(shape(3),sslen+dk)
                % get all XY such that k0+dk==kk
                idx=ss(kk-dk).PixelIdxList;
                L=map_upper(idx);
                idx=idx(L);                

                if(~isempty(idx))                    
                    dtmp=D{k}(idx)+aspect(3)^2*(kk-k)^2;
    
                    % these pixels are to be updated with k0+dk
                    L=D1{kk}(idx)>dtmp & D{kk}(idx)>0;
                    map_upper(idx)=L;
                    D1{kk}(idx(L))=dtmp(L);
                    
                    dtmp=idx(L); 
                    dtmp(:)=k;
                    DK{kk}(idx(L))=dtmp;
                end
            end
            idx_upper=idx_upper(map_upper(idx_upper));
        else
            eols=eols-1;
        end
    end
end

% the answer
if(iscell(bw))
    D=cell(size(bw));
    for k=1:shape(3) D{k}=sqrt(D1{k}); end
else
    D=zeros(shape);
    for k=1:shape(3) D(:,:,k)=sqrt(D1{k}); end
end


function s=getregions(map)
% this function is replacer for regionprops(map,'PixelIdxList');
% it produces the list of different values along with list of 
% indexes of pixels in map with these values; 's' is struct-array 
% such that s(i).PixelIdxList contains list of pixels in map 
% with value i.

% enable using regionprops if available, faster on 7.3
fregionprops=1;

% version control for using regionprops
v=version; v=str2num(v(1:3)); 
fregionprops=fregionprops & v>=7.3;

% in later matlab regionprops is actually faster than this code
if(exist('regionprops') & fregionprops)
    s=regionprops(map,'PixelIdxList');
    return
end

idx=(1:prod(size(map)))';
dtmp=double(map(:));

[dtmp,ind]=sort(dtmp); 
idx=idx(ind);
ind=[0;find([diff(dtmp(:));1])];

s=[];
for i=2:length(ind)
    if(dtmp(ind(i)))==0 continue; end
    s(dtmp(ind(i))).PixelIdxList=idx(ind(i-1)+1:ind(i));
end