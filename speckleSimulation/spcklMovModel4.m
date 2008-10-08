function [params,fluorData,not2keep]=spcklMovModel4(params,fluorData,not2keep)
% antiparallel microtubules in bundles
% initial positions and veloctiy assignments handled by spckleMovGenFluor

dT=params.dT;
vx=params.mod4vx;

% initialize temporary variables for state and coodinates
state=[fluorData.state(:,1) zeros(size(fluorData.state(:,1)))]; % nFluor x 2
py1=fluorData.py(:,1); % nFluor-vector
px1=fluorData.px(:,1); % nFluor-vector

% loop through every time point in the movie
for timePt=1:not2keep.nTmPtsInMovie
    % t is empty if current time point is not one of the ones to store;
    % otherwise it equals the current timePt.
    t=intersect(timePt,not2keep.pts2Store);  
    if ~isempty(t) % store the state/position info if t isn't empty
        col=find(not2keep.pts2Store==t); % column in which to store data
        fluorData.state(:,col)=state(:,1);
        fluorData.py(:,col)=py1;
        fluorData.px(:,col)=px1;
    end

   px1=vx*dT+px1;
   state(:)=1; 
end