function D=bwdistX(bw,aspect)
% D=BWDISTX(BW,ASPECT)
% BWDISTX computes Euclidean distance transform of the binary 3D image 
% BW. For each pixel in BW, the distance transform assignes a number
% that is the distance between that pixel and the nearest nonzero pixel 
% of BW. BW may be a single 2D image, 3D array or a cell array of 
% 2D slices. ASPECT is 3-component vector defining voxel-aspect-ratio in 
% the dataset BW. If ASPECT is not specified, [1 1 1] isotropic aspect 
% ratio is assumed.
%
% BWDISTX uses fast optimized scan algorithm and cell-arrays to 
% represent internal data, and is less demanding to physical memory and 
% in many cases up to 10 times faster than MATLAB's native bwdist.
%
% This is an experimental version using an optimized forward-backward 
% scan version of the algorithm used in the original bwdistsc. It will 
% be developed into a full version and replace the older bwdistsc in the 
% package in the near future. However, currently there are following 
% limitations:
% - ASPECT should be of the form [a a b], that is only anisotropy in 
%   z-direction is currently supported;
% - image processing toolbox should be already present on your machine.
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
if(nargin<2) aspect=[1 1 1]; end

% geometry of data
shape=size(bw);

% allocate space
D=cell(1,shape(3)); 
for k=1:shape(3) 
  D{k}=zeros(shape(1:2)); 
end

%%%%%%%%%%%%% scan along XY %%%%%%%%%%%%%%%%
for k=1:shape(3)    
    D{k}=aspect(1)^2*bwdist(bw(:,:,k)).^2;
end

%%%%%%%%%%%%% scan along Z %%%%%%%%%%%%%%%%
D1=cell(size(D));
for k=1:shape(3) 
  D1{k}=repmat(Inf,shape(1:2)); 
end

% start building the envelope 
p=shape(3);
for k=1:shape(3)
    % if there are no objects in this slice, nothing to do
    if(isinf(D{k}(1,1)))
      continue;
    end
    
    % selecting starting point for (x,y):
    % * if parabolas are incremented in increasing order of k, then all 
    %   intersections are necessarily at the right end of the envelop, 
    %   and so the starting point can be always chosen as the right end
    %   of the axis
    
    % check which points are valid starting points, & update the envelop
    dtmp=D{k}+aspect(3)^2*(p-k)^2;
    L=D1{p}>dtmp; 
    D1{p}(L)=dtmp(L);    
    
    % map_lower keeps track of which pixels can be yet updated with the 
    % new distance, i.e. all such XY that had been under the envelop for
    % all Deltak up to now, for Deltak<0
    map_lower=L;
        
    % these are maintained to keep fast track of whether map is empty
    idx_lower=find(map_lower);
    
    % scan away from the starting points in increments of -1
    for kk=p-1:-1:1
        % new values for D
        dtmp=D{k}(idx_lower)+aspect(3)^2*(kk-k)^2;
                    
        % these pixels are to be updated
        L=D1{kk}(idx_lower)>dtmp;
        map_lower(idx_lower)=L;
        D1{kk}(idx_lower(L))=dtmp(L);
                    
        % other pixels are removed from scan
        idx_lower=idx_lower(L);
        
        if(isempty(idx_lower)) break; end
    end
end

% the answer
D=zeros(shape);
for k=1:shape(3) D(:,:,k)=sqrt(D1{k}); end