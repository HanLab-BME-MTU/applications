function [params,fluorData,not2keep]=spcklMovModel2(params,fluorData,not2keep)
% speckles converge with depoly such that avg I is constant, or with no
% depoly such that I gets brighter at pole

nmPerSecFlowSpeed=not2keep.nmPerSecFlowSpeed; % actin flow nm/s

if isempty(params.poleYX) % let default pole location be in center of image
    y0=not2keep.imgLnm/2;
    x0=not2keep.imgWnm/2;
else % use location given by user  
    y0pix=params.poleYX(1)+params.border.top;
    x0pix=params.poleYX(2)+params.border.left;
    y0=y0pix*params.pixNM;
    x0=x0pix*params.pixNM;
end
    
% initialize temporary variables for state and coodinates
state=fluorData.state(:,1); % nFluor x 2
py1=fluorData.py(:,1); % nFluor-vector
px1=fluorData.px(:,1); % nFluor-vector

for timePt=1:not2keep.nTmPtsInMovie
    % t is empty if current time point is not one of the ones to store;
    % otherwise it equals the current timePt.
    t=intersect(timePt,not2keep.pts2Store);  
    if ~isempty(t) % store the state/position info if t isn't empty
        col=find(not2keep.pts2Store==t); % column in which to store data
        fluorData.state(:,col)=state;
        fluorData.py(:,col)=py1;
        fluorData.px(:,col)=px1;
    end
    
    if params.depoly==1 % depoly should occur to keep avg I constant
        r=sqrt((px1-x0).^2+(py1-y0).^2); % distance of each pt from the pole (nm)
        P=zeros(size(r));
        P=(nmPerSecFlowSpeed./r)*params.dT;
        R=rand(size(P));
        PoverR=P./R; % ratio of prob. to random numbers
        state(PoverR>1)=0; % if prob > random number, let depoly happen
    elseif params.depoly~=0 % if depoly=0, don't do anything
        disp('depoly must be 1 or 0');
        return
    end
    
    % calculate velocity of each fluorophore based on position relative to
    % pole (all have the same speed, but the direction points inwards)
    vx=nmPerSecFlowSpeed*(-(px1-x0))./sqrt((px1-x0).^2+(py1-y0).^2);
    vy=nmPerSecFlowSpeed*(-(py1-y0))./sqrt((px1-x0).^2+(py1-y0).^2);
    vx(px1==x0 & py1==y0)=0; %vel=0 when fluorophore falls into pole
    vy(px1==x0 & py1==y0)=0;

    % new coordinates of fluorophores
    px1=px1+vx*params.dT; 
    py1=py1+vy*params.dT;
end

