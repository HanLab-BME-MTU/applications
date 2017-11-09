function spots = spotfind_new(fImg,dataProperties)
%SPOTFIND_NEW locates fluorescent tags in 3D data
%
% The creteria of spots selection is based on mnp threshold, instead of
% total number limitation (MAXSPOTS)
% SYNOPSIS cord = spotfind(img)
%
% INPUT img   : stack time series
%       dataProperties: structure with movie properties
%
% OUTPUT spots : nTimepoints-by-1 structure with fields
%                .sp  structure with fields
%                   .cord coordinates
%                   .mnint spottiness
%                .COM center of mass of image

% 04/2015 Ning

%CONST DEFINITIONS
%global PATCHSIZE;
%PATCHSIZE=dataProperties.PATCHSIZE;
FILTERSIZE = dataProperties.FILTERPRM(4:6);
PATCHSIZE = FILTERSIZE;
% DEBUG = 0;

% init vars
d=floor(PATCHSIZE/2);
% inTestD = floor(FILTERSIZE/2); %number of pixels a spot has to be away from the border to be accepted
inTestD = [3,3,1];
movieSize = size(fImg);
tsteps=size(fImg,5); % just in case there's only one frame loaded


%preassign mnp. provide for 100 points
mnpRows = 100;
mnpRowIncrement = 100;
mnp = zeros(mnpRows,tsteps);

% preassign spots
spots(1:tsteps,1) = struct('sp',[],'COM',[]);

%loop through all time points
for t=1:tsteps

    %intialize counter
    ct=1;
    % current time point:
    mnplist=[];    %'spottiness'
    lst=[];            % list of local maxs
    k=[];             % curvature of local maxs

    pt=fImg(:,:,:,1,t);

    %norm to 0..100
    pt=100*pt/max(pt(:));
    
    

    %find all local max
    b=loc_max3Df(fImg(:,:,:,1,t),[3 3 3]);

    [FXX,FXY,FXZ,FYX,FYY,FYZ,FZX,FZY,FZZ]=hessian(pt); % hessian matrix of full intensity dist.

    %loop through all local maxs
    for i=1:size(b,1)
        %ignore pixels close to border
        if(all((b(i,:)-inTestD)>0) && all((b(i,:)+inTestD)...
                <=[movieSize(1:3)]))

            %cut pixels belonging to this local maximum
            patch = stamp3d(pt,PATCHSIZE,b(i,:),0);
%             patch=pt(b(i,1)-d(1):b(i,1)+d(1),b(i,2)-d(2):b(i,2)+d(2),b(i,3)-d(3):b(i,3)+d(3));


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
                
                % mnp, the spottiness criterion, is the product of
                % curvature and mean intensity of the local mask.
                % The cutoff might be a bit nicer if we transformed
                % curvature and intensity to [0,1] first, but I leave it
                % for the moment.
                mnp(ct,t)=-k(ct)*mean(patch(:));
                centp(ct,:)=centroid3D(patch);
                lm(ct,:)=b(i,:);
                ct=ct+1;
                
                
            end;
        end;
    end;
        
    % Take MAXNUMSPOTS plus a few spots - we want to be sure that we don't
    % accidentially throw out a good spot, and we need a few bad apples to
    % make the amplitude cutoff work fine. We take between 2 and 10 more
    % spots, depending on MAXNUMSPOTS
    
%     additionalSpots = dataProperties.MAXSPOTS * 0.3;
%     additionalSpots = ceil(additionalSpots);
%     additionalSpots = max(additionalSpots,3);
%     additionalSpots = min(additionalSpots,10);
%     numberOfSpots = dataProperties.MAXSPOTS + additionalSpots;

    % Choose qualified spots number based on mnp value
    [mnpSorted,sortIdx] = sort(mnp(1:ct-1,t),1,'descend');
    
    % Plot cumulative histogram for mnp
    LM = zeros(size(mnpSorted,1),1);
    for i=1:size(LM)
        LM(i) = size(mnpSorted,1)-i+1;
    end
    plot(mnpSorted,LM,'r+');
    
    
%     mnpThreshold = input('Enter spottiness threshold > ');
% Carefully choose and modify mnpThreshold 01/2016
%     mnpThreshold = 0.001;
    mnpThreshold = 1e5;
    
    for qualifiedNum = 1:size(mnpSorted,1)
        if mnpSorted(qualifiedNum) < mnpThreshold
            qualifiedNum = qualifiedNum - 1;
            break
        else
            qualifiedNum = qualifiedNum + 1;
        end
    end
    numberOfSpots = qualifiedNum;
    % cut at either MAXSPOTS+1 or how many we have if it's less
    cps = sortIdx(1:min(numberOfSpots,length(sortIdx)));
    
    
    if cps~=0
        lst=[lm(cps,2) lm(cps,1) lm(cps,3)]-ones(length(cps),1)*(d+1)+centp(cps,:);
        mnplist=mnp(cps,t);
        j = 1;
        for i=1:size(lst,1)
            % Set boundary for spots center
            if lst(i,:) > PATCHSIZE/3
                if lst(i,:) < [size(fImg,2),size(fImg,1),size(fImg,3)]-PATCHSIZE/3
                    % store coordinates and spottiness
                    spots(t).sp(j).cord=lst(i,:);
                    spots(t).sp(j).mnint=mnplist(i);
                    j = j+1;
                end
            end
        end;
    end
%     spots(t).mnint=mn(t);
    
    %== ADDED FOR NEW LINKER ==
    % Find the "center of gravity" of the image
    % use int^10 to get good results (maybe we need more for mammalian
    % cells with their bigger frames?)
    
    % Seems useless
    spots(t).COM = 0;
    
%     spots(t).COM = centroid3D(pt,10);
        
    % clean memory
    clear FXX FXY FXZ FYX FYY FYZ FZX FZY FZZ
    
end;