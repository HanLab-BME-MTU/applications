function [cord, mnp] = spotfind(fImg,dataProperties)
%SPOTFIND locates fluorescent tags in 3D data
%
% SYNOPSIS cord = spotfind(img)
%
% INPUT img   : stack time series
%
% OUTPUT cord : center coordinates 

% c: 5/3/01	dT

%CONST DEFINITIONS
%global PATCHSIZE;
%PATCHSIZE=dataProperties.PATCHSIZE;
FILTERSIZE = dataProperties.FILTERPRM(4:6);
PATCHSIZE = FILTERSIZE;

% init vars
d=floor(PATCHSIZE/2); 
inTestD = floor(FILTERSIZE/2); %number of pixels a spot has to be away from the border to be accepted

%wb_hdl=waitbar(0,'Finding spots...');
tsteps=size(fImg,5);
h= mywaitbar(0,[],tsteps,'Finding spots...');

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
                <=[size(pt,1) size(pt,2) size(pt,3)]))
            
            %cut pixels belonging to this local maximum
            patch=pt(b(i,1)-d(1):b(i,1)+d(1),b(i,2)-d(2):b(i,2)+d(2),b(i,3)-d(3):b(i,3)+d(3));
            
            %curvature filter
            %k(ct)=curvature3D(patch,[d d d]+1);
            k(ct)=det([FXX(b(i,1),b(i,2),b(i,3)) FXY(b(i,1),b(i,2),b(i,3)) FXZ(b(i,1),b(i,2),b(i,3));...
                            FYX(b(i,1),b(i,2),b(i,3)) FYY(b(i,1),b(i,2),b(i,3)) FYZ(b(i,1),b(i,2),b(i,3));...
                            FZX(b(i,1),b(i,2),b(i,3)) FZY(b(i,1),b(i,2),b(i,3)) FZZ(b(i,1),b(i,2),b(i,3))]);

            % only convex shapes allowed
            if k(ct) < 0
                mnp(ct,t)=-k(ct)*mean(patch(:));
                centp(ct,:)=centroid3D(patch);
                lm(ct,:)=b(i,:);
                ct=ct+1;
            end;
        end;
    end;
    
    %cumulative histogram spot separation
    %run only if more than 2 max found
    if length(find(mnp(:,t))) > 2
        cps=cutcumhist(mnp,t,dataProperties);
    elseif isempty(mnp)
        cps = 0;
    else
        cps = find(mnp(:,t)>100); %as in cutcumhist: spottiness has to be at least 100
    end
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
    
    mywaitbar(t/tsteps,h,tsteps);
    %waitbar(t/tsteps,wb_hdl);
end;
cord=spots;
close(h);
