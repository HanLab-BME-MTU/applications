function [cord, mnp] = spotfind(fImg,dataProperties,verbose)
%SPOTFIND locates fluorescent tags in 3D data
%
% SYNOPSIS cord = spotfind(img)
%
% INPUT img   : stack time series
%       dataProperties: structure with movie properties
%       verbose : (optional) If 1 (default), waitbar is displayed
%
% OUTPUT cord : center coordinates

% c: 5/3/01	dT

% optional input arguments
if nargin < 3 || isempty(verbose)
    verbose = 1;
end

%CONST DEFINITIONS
%global PATCHSIZE;
%PATCHSIZE=dataProperties.PATCHSIZE;
FILTERSIZE = dataProperties.FILTERPRM(4:6);
PATCHSIZE = FILTERSIZE;
DEBUG = 0;

% init vars
d=floor(PATCHSIZE/2);
inTestD = floor(FILTERSIZE/2); %number of pixels a spot has to be away from the border to be accepted
movieSize = size(fImg);
tsteps=movieSize(5);

if verbose
    h= mywaitbar(0,[],tsteps,'Finding spots...');
end

%preassign mnp. provide for 100 points
mnpRows = 100;
mnpRowIncrement = 100;
mnp = zeros(mnpRows,tsteps);

if DEBUG
    figure
end

%loop through all time points
for t=1:tsteps

    %intialize counter
    ct=1;
    % current time point:
    mnplist=[];    %'spotiness'
    lst=[];            % list of local maxs
    k=[];             % curvature of local maxs

    pt=fImg(:,:,:,1,t);
    mn(t)=mean(pt(:));

    %norm to 0..1
    pt=100*pt/max(pt(:));

    %find all local max
    b=loc_max3Df(fImg(:,:,:,1,t),[3 3 3]);

    [FXX,FXY,FXZ,FYX,FYY,FYZ,FZX,FZY,FZZ]=hessian(pt); % hessian matrix of full intensity dist.

    %loop through all local maxs
    for i=1:size(b,1)
        %ignore pixels close to border
        if(all((b(i,:)-inTestD)>0) & all((b(i,:)+inTestD)...
                <=[movieSize(1:3)]))

            %cut pixels belonging to this local maximum
            patch=pt(b(i,1)-d(1):b(i,1)+d(1),b(i,2)-d(2):b(i,2)+d(2),b(i,3)-d(3):b(i,3)+d(3));

            %curvature filter
            %k(ct)=curvature3D(patch,[d d d]+1);
            k(ct)=det([FXX(b(i,1),b(i,2),b(i,3)) FXY(b(i,1),b(i,2),b(i,3)) FXZ(b(i,1),b(i,2),b(i,3));...
                FYX(b(i,1),b(i,2),b(i,3)) FYY(b(i,1),b(i,2),b(i,3)) FYZ(b(i,1),b(i,2),b(i,3));...
                FZX(b(i,1),b(i,2),b(i,3)) FZY(b(i,1),b(i,2),b(i,3)) FZZ(b(i,1),b(i,2),b(i,3))]);

            % only convex shapes allowed
            if k(ct) < 0
                if ct > mnpRows
                    % reassign mnp-matrix
                    mnpTmp = mnp;
                    newMnpRows = mnpRows + mnpRowIncrement;
                    mnp = zeros(newMnpRows,tsteps);
                    mnp(1:mnpRows,:) = mnpTmp;
                    mnpRows = newMnpRows;
                    
                    clear mnpTmp newMnpRows
                end
                
                mnp(ct,t)=-k(ct)*mean(patch(:));
                centp(ct,:)=centroid3D(patch);
                lm(ct,:)=b(i,:);
                ct=ct+1;
            end;
        end;
    end;

    % NEW: Use dataProperties.MAXSPOTS number of points. Because I don't
    % want to go into the editPropertiesGUI at the moment (and because
    % we're not working with mammalian stuff right now), I'll just use the
    % default 5. Later, I'll use something like ceil(maxSpots * 1.3).
    
    [mnpSorted,sortIdx] = sort(mnp(:,t),1,'descend');
    % cut at either MAXSPOTS or how many we have if it's less
    cps = sortIdx(1:min(dataProperties.MAXSPOTS,length(sortIdx)));
    
    
%     %cumulative histogram spot separation
%     %run only if more than 2 max found
%     if length(find(mnp(:,t))) > 2
%         cps=cutcumhist(mnp,t,dataProperties);
%     elseif isempty(mnp)
%         cps = 0;
%     else
%         cps = find(mnp(:,t)>100); %as in cutcumhist: spottiness has to be at least 100
%     end
    
   
    
    
    if cps~=0
        lst=[lm(cps,2) lm(cps,1) lm(cps,3)]-ones(length(cps),1)*(d+1)+centp(cps,:);
        mnplist=mnp(cps,t);
        %    lst=[lst(:,2) lst(:,1) lst(:,3)];
        for i=1:size(lst,1)
            %initalize and store info in struct
            spots(t).sp(i).type='';
            spots(t).sp(i).mult=0;
            spots(t).sp(i).cord=lst(i,:);
            spots(t).sp(i).mnint=mnplist(i);
            % include curvature
            %spots(t).sp(i).k=k(i);
        end;
    else
        spots(t).sp=[];
    end
    spots(t).mnint=mn(t);
    
    %== ADDED FOR NEW LINKER ==
    % Find the "center of gravity" of the image
    % use int^10 to get good results (maybe we need more for mammalian
    % cells with their bigger frames?)
    [x,y,z] = meshgrid(1:movieSize(1),1:movieSize(2),1:movieSize(3));
    pt10 = pt(:).^10;
    sumFrame = sum(pt10);
    centerX = x(:)'*pt10/sumFrame;
    centerY = y(:)'*pt10/sumFrame;
    centerZ = z(:)'*pt10/sumFrame;
    
    spots(t).COM = [centerX, centerY, centerZ];
    

    if verbose
        mywaitbar(t/tsteps,h,tsteps);
    end
end;
cord=spots;
if verbose
    close(h);
end
