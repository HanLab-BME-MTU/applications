function [transToEdgeDistNearest,transToEdgeDistMech]=fsmTransTransitionToEdgeDistances(spTrans_y,spTrans_x,spEdge_y,spEdge_x,numSplinePoints,MECH)
% fsmTransTransitionToEdgeDistances calculates the distances between the transition splines and the edge splines
%
% SYNOPSIS    [transToEdgeDistNearest,transToEdgeDistMech]=
%                  fsmTransTransitionToEdgeDistances(spTrans_y,spTrans_x,spEdge_y,spEdge_x,numSplinePoints,MECH)
%
% INPUT       spTrans_y       : spline structure for the y coordinates of the transition spline (to be used with fnval)
%             spTrans_x       : spline structure for the s coordinates of the transition spline (to be used with fnval)
%             spEdge_y        : spline structure for the y coordinates of the edge spline (to be used with fnval)
%             spEdge_x        : spline structure for the s coordinates of the edge spline (to be used with fnval)
%             numSplinePoints : number of knots along the splines
%             MECH            : [ 0 | 1 ] Turns off | on the use of a mechanical model to calculate distances between the splines 
%
% OUTPUT      transToEdgeDistNearest : structure (t = 1:number of frames).y  : y coordinate of the knot of the spline at time t
%                                                                 .x  : x coordinate of the knot of the spline at time t
%                                                                 .dy : y displacement at .y of the knot of the spline from time t to t+1
%                                                                 .dx : x displacement at .x of the knot of the spline from time t to t+1
%                                      The distances are calculated by simple nearest neighbor.
%             transToEdgeDistMech    : same as transToEdgeDistNearest, but calculated with a mechanical model.
%
% DEPENDENCES   fsmTransTransitionToEdgeDistances uses { }
%               fsmTransTransitionToEdgeDistances is used by { fsmTransition }
%
% Aaron Ponti, September 6th, 2004

% Constants
rmPoints=0.05;  % Fraction of points along the spline removed from both ends
points=1:numSplinePoints;

% Initialize waitbar
wHandle=waitbar(0,'Please wait...');

% Number of frames to be processed
numIterations=length(spEdge_y);

% Initialize output
transToEdgeDistNearest(1:numIterations)=struct('y',[],'x',[],'dy',[],'dx',[]);
transToEdgeDistMech=transToEdgeDistNearest;

% Go through frames
for i=1:numIterations
    
    % Calculate distances transition to leading edge - first 'nearest'
    offset=fix(rmPoints*length(points));      % Remove (rmPoints*100)% of the points at both boundaries
    rn=(offset+1):length(points)-offset;      % (to prevent boundary problems) 
    [baseN,dispN,r0]=prGetDispNearest(spTrans_y(i),spTrans_x(i),spEdge_y(i),spEdge_x(i),rn,'robust_min',1,'tol',0);

    % Store distances
    transToEdgeDistNearest(i).y=baseN(:,1);
    transToEdgeDistNearest(i).x=baseN(:,2);
    transToEdgeDistNearest(i).dy=dispN(:,1);
    transToEdgeDistNearest(i).dx=dispN(:,2);

    % Check r0 (in case there were still boundary problems)
    [r0,rn]=checkNearest(r0,rn,max(rn));
    
    if MECH==1
        
        % Then with a mechanical model (use output r0 of 'nearest' as staring point)
        try
            [base,disp]=prGetDispMech(spTrans_y(i),spTrans_x(i),spEdge_y(i),spEdge_x(i),rn,r0,0);
        catch
            
            % The mechanical model failed: use the result from 'nearest'
            fprintf(1,'Frame %3d: prGetDispMech failed - using the distance obtained from prGetDispNearestMulti.\n',i);
            base=baseN;
            disp=dispN;
            
        end
        
    else
        
        % Store the results of 'nearest'
        base=baseN;
        disp=dispN;
        
    end    
    
    % Store distances
    transToEdgeDistMech(i).y=base(:,1);
    transToEdgeDistMech(i).x=base(:,2);
    transToEdgeDistMech(i).dy=disp(:,1);
    transToEdgeDistMech(i).dx=disp(:,2);

    % Update waitbar
    waitbar(i/numIterations,wHandle)
    
end

close(wHandle);

% warning on MATLAB:divideByZero

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r0,rn]=checkNearest(r0,rn,maxPix)

lenRn=length(r0);
if r0(lenRn)>maxPix
    i=lenRn;
    while r0(i)>maxPix & length(r0)>0
        r0(i)=[];
        rn(i)=[];
        i=i-1;
    end                      
end  
if r0(1)<1
    i=1;
    while r0(i)<1 & length(r0)>0
        r0(i)=[];
        rn(i)=[];
    end                      
end  

