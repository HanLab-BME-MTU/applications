function [distancesNearest,distancesMech]=fsmTransProtrusions(spl_y,spl_x,points,MECH)
% fsmTransProtrusions calculates the distances (-> protrusion) between two splines 
%
% SYNOPSIS    [distancesNearest,distancesMech]=fsmTransProtrusions(spl_y,spl_x,points,MECH)
%
% INPUT       spl_y           : spline structure for the y coordinates of the curve described by the spline (to be used with fnval)
%             spl_x           : spline structure for the x coordinates of the curve described by the spline (to be used with fnval)
%             points          : indices of the knots along the spline
%             MECH            : [ 0 | 1 ] Turns off | on the use of a mechanical model to calculate distances between the splines 
%
% OUTPUT      distanceNearest : structure (t = 1:number of frames).y  : y coordinate of the knot of the spline at time t
%                                                                 .x  : x coordinate of the knot of the spline at time t
%                                                                 .dy : y displacement at .y of the knot of the spline from time t to t+1
%                                                                 .dx : x displacement at .x of the knot of the spline from time t to t+1
%                               The distances are calculated by simple nearest neighbor.
%             distanceMech    : same as distanceNearest, but calculated with a mechanical model.
%
% DEPENDENCES   fsmTransProtrusions uses { }
%               fsmTransProtrusions is used by { fsmTransition }
%
% Aaron Ponti, September 6th, 2004

% Check input
if nargin~=4
    error('Four input parameters expected.');
end
if length(spl_y)<2
    error('At least two frames are required to calculate a protrusion.');
end

% Constant
rmPoints=0.05;  % Fraction of points along the spline removed from both ends

% Initialize outputs
distancesNearest=struct('y',[],'x',[],'dy',[],'dx',[]);
distancesMech=distancesNearest;

% Set points along the spline to calculate distances
if nargin==2
    n=max(unique(spl_y(1).knots));
    points=1:n;
            
else
    if size(points)==[1 1]
        points=1:points;   
    end
    rmPoints=0;
       
end
offset=fix(rmPoints*length(points)); % Number of points to be removed from bpth ends

% Initialize waitbar 
hWait=waitbar(0,'Please wait...');

% Go through frames
numFrames=length(spl_x)-1;

for i=1:numFrames
    
    rn=points(offset+1:end-offset);  % Remove (rmPoints)% of the points at both boundaries
                                     % (to prevent boundary problems) 

    % Calculate distances transition to leading edge - first 'nearest'
    [baseN,dispN,r0]=prGetDispNearest(spl_y(i),spl_x(i),spl_y(i+1),spl_x(i+1),rn,'robust_min',1,'tol',0);

    % Store distances
    distancesNearest(i).y=baseN(:,1);
    distancesNearest(i).x=baseN(:,2);
    distancesNearest(i).dy=dispN(:,1);
    distancesNearest(i).dx=dispN(:,2);

    % Check r0 (in case there were still boundary problems)
    [r0,rn]=checkNearest(r0,rn,max(rn));
    
    if MECH==1
        
        % Then with a mechanical model (use output r0 of 'nearest' as starting point)
        try
            [base,disp]=prGetDispMech(spl_y(i),spl_x(i),spl_y(i+1),spl_x(i+1),rn,r0,0);
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
    distancesMech(i).y=base(:,1);
    distancesMech(i).x=base(:,2);
    distancesMech(i).dy=disp(:,1);
    distancesMech(i).dx=disp(:,2);
    
    % Display some info
    waitbar(i/numFrames,hWait)
end

% Close waitbar
close(hWait);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
