function [spTrans_y,spTrans_x,spEdge_y,spEdge_x,splineCoordsTrans,splineCoordsEdge,numSplinePoints]=fsmTransFitSplines(positionsTrans,edgePixels,splineTolEdge,splineTolTrans,medianFilterFrames,extractTransitionFrames)
% fsmTransFitSplines fits splines through edge and transition as extracted by prPanel and fsmTransExtractLpToLaTransition, respectively.
%
% SYNOPSIS      [spTrans_y,spTrans_x,spEdge_y,spEdge_x,splineCoordsTrans,splineCoordsEdge,numSplinePoints]= ...
%                    fsmTransFitSplines(positionsTrans,edgePixels,splineTolEdge,splineTolTrans,medianFilterFrames,extractTransitionFrames)
%
% INPUT         positionsTrans          : transition coordinates for all frames calculated by fsmTransExtractLpToLaTransition 
%                                         size = [ number of profiles x 2 x frames]
%               edgePixels              : structure(1:k).x : x coordinates of pixel edges 
%                                                       .y : y coordinates of pixel edges
%                                         this structure is returned by fsmReadPrPanelDatFiles (which imports prPanel data(
%               splineTolEdge           : tolerance for spline through edgePixels
%               splineTolTrans          : tolerance for spline through transition
%               medianFilterFrames      : order of the media filter (to smoothen the obtained transition line)
%               extractTransitionFrames : difference in frames due to the time-averaging of the transition; if frames 1 to 5 are used 
%                                         to average on transition, the resulting line should be then coupled to the edge line of frame 3.
%                                         In this specific case, if extractTransitionFrames is set to 5 (frames), the first spline for the 
%                                         transitions will correspond to the average of frames 1 - 5, the second to the average of frames 2 - 6,
%                                         and so on. The first spline for the edges will fit the edge of frame 3, the second will fit the edge 
%                                         of frame 4, and so on.
%                                         transitions
%                                         
% OUTPUT        spTrans_y              : spline structure for the y coordinates of the transitions (to be used with fnval)
%               spTrans_x              : spline structure for the x coordinates of the transitions (to be used with fnval)
%               spEdge_y               : spline structure for the y coordinates of the edges (to be used with fnval)
%               spEdge_x               : spline structure for the x coordinates of the edges (to be used with fnval)
%               splineCoordsTrans      : structure (1:number of frames).y = y coordinates of knots along the transition spline
%                                                                      .x = x coordinates of knots along the transition spline
%               splineCoordsEdge       : structure (1:number of frames).y = y coordinates of knots along the edge spline
%                                                                      .x = x coordinates of knots along the edge spline
%               numSplinePoints        : number of knots along the spline
%
% DEPENDENCES   fsmTransFitSplines uses { }
%               fsmTransFitSplines is used by { fsmTransition }
%
% Aaron Ponti, September 6th, 2004

% Turn on and off debugging mode
DEBUG=0;

if nargin~=6
    error('Six input parameters expected.');
end

% Initialize spline structures
spEdge_y=struct('form',[],'knots',[],'coefs',[],'number',[],'order',[],'dim',[]);
spEdge_x=spEdge_y;
spTrans_y=spEdge_y;
spTrans_x=spEdge_y;
splineCoordsTrans=struct('y',[],'x',[]);
splineCoordsEdge=splineCoordsTrans;

% Initialize waitbar
wHandle = waitbar(0,'Please wait...');

% Number of frames to be processed
numIterations=size(positionsTrans,3);

% Time offset between edge and transition (due to the averaging of the transition 
%    over 'extractTransitionFrames' frames)
dFrames=fix(extractTransitionFrames/2);

% Go through frames
for i=1:numIterations
    
    % Since the splines are not normalized, I need to subsample the spline with
    %   the largest number of points
    edgePoints=length(edgePixels(i+dFrames).y);
    transPoints=size(positionsTrans,1);
    
    if edgePoints>transPoints
        
        % Calculate steps for subsampling
        step=fix(edgePoints/transPoints);
        edgePoints=1:step:edgePoints;
        edgePoints=edgePoints(1:transPoints);
        
        % Subsample
        pixelEdgeY=edgePixels(i+dFrames).y(edgePoints);
        pixelEdgeX=edgePixels(i+dFrames).x(edgePoints);
        
        % Keep unchanged
        positionsTransY=positionsTrans(:,1,i);
        positionsTransX=positionsTrans(:,2,i);
        
    elseif transPoints>edgePoints
        
        % Calculate steps for subsampling
        step=fix(transPoints/edgePoints);
        transPoints=1:step:transPoints;
        transPoints=edgePoints(1:edgePoints);
        
        % Subsample
        positionsTransY=positionsTrans(:,1,i);
        positionsTransY=positionsTransY(transPoints);
        positionsTransX=positionsTrans(:,2,i);
        positionsTransX=positionsTransX(transPoints);
        
        % Keep unchanged
        pixelEdgeY=edgePixels(i+dFrames).y;
        pixelEdgeX=edgePixels(i+dFrames).x;
        
    else
        
        pixelEdgeY=edgePixels(i+dFrames).y;
        pixelEdgeX=edgePixels(i+dFrames).x;
        positionsTransY=positionsTrans(:,1,i);
        positionsTransX=positionsTrans(:,2,i);
        
    end        
    
    % Calculate spline through current edgePixels
    [spEdge_y(i),spEdge_x(i)]=imPixelChainSpline([pixelEdgeY pixelEdgeX],'tolerance',splineTolEdge);
    
    % Get spline coordinates through the edge
    points=1:length(pixelEdgeY);
    splineCoordsEdge(i).y=fnval(spEdge_y(i),points);
    splineCoordsEdge(i).x=fnval(spEdge_x(i),points);
    
    % Calculate spline through current positionTrans
    posY=medfilt1(positionsTransY,medianFilterFrames);
    posX=medfilt1(positionsTransX,medianFilterFrames);
    [spTrans_y(i),spTrans_x(i)]=imPixelChainSpline([posY posX],'tolerance',splineTolTrans);
    
    % Get spline coordinates through the transition
    points=1:size(positionsTrans,1);
    splineCoordsTrans(i).y=fnval(spTrans_y(i),points);
    splineCoordsTrans(i).x=fnval(spTrans_x(i),points);
    
    % Plot edge (spline) and transition: spline, raw, and median filtered
    if DEBUG==1
        h=figure;
        plot(splineCoordsEdge(i).x,splineCoordsEdge(i).y,'r-');
        hold on
        plot(splineCoordsTrans(i).x,splineCoordsTrans(i).y,'k-');
        plot(positionsTransX,positionsTransY,'g-');
        plot(posX,posY,'m-');
        axis equal
        legend('Spline edge','Spline trans','Raw trans','Raw trans, median filtered',0);
    end
    
    % Update waitbar
    waitbar(i/numIterations,wHandle)
    
end

close(wHandle);

% Return the number of points 
numSplinePoints=min(transPoints,length(edgePoints));


