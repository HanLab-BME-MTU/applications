function lpPolygons=fsmTransCreateLpPolygons(splineCoordsEdge,splineCoordsTrans,showPolygons)
% fsmTransCreateLpPolygons creates a series of lamellipodium polygon by connecting the edge and transition profiles
%
% SYNOPSIS   lpPolygons=fsmTransCreateLpPolygons(splineCoordsEdge,splineCoordsTrans,showPolygons)
%
% INPUT      splineCoordsEdge   : splineCoordsEdge as returned by fsmTransFitSplines
%            splineCoordsTrans  : splineCoordsTrans as returned by fsmTransFitSplines
%            showPolygons       : [ 0 | 1 ] toggled display of created polygons (default 0 : off)
%
% OUTPUT     lpPolygons          : lamellipodium polygons, with structure lpPolygons(1:number of frames).y : y coordinates
%                                                                                                      .x : x coordinates
%
% REMARK     THIS FUNCTION ASSUMES THAT THE DIFFERENCE IN FRAMES DUE TO THE TIME-AVERAGING OF THE TRANSITION HAS ALREADY BEEN TAKEN INTO
%            ACCOUNT WHEN THE STRUCTURE OF SPLINES HAS BEEN GENERATED (this is the case if the user interface fsmTransition has been used).
%
% DEPENDENCES   fsmTransCreateLpPolygons uses { }
%               fsmTransCreateLpPolygons is used by { fsmTransition }
%
% Aaron Ponti, September 6th, 2004

if nargin<2 | nargin>3
    error('Two or three input parameters expected');
end

if nargin==2
    showPolygons=0;
end

% Read number of frames
nm=length(splineCoordsTrans);

% Initialize lpPolygons structure
lpPolygons(1:nm)=struct('y',[],'x',[]);

% Initialize variables to store min and max y and x coordinates
MNY=[]; MXY=[]; MNX=[]; MXX=[];

% Create and store polygons
for i=1:nm
    
    y=[splineCoordsEdge(i).y,splineCoordsTrans(i).y(end:-1:1),splineCoordsEdge(i).y(1)]; % Close polygon
    x=[splineCoordsEdge(i).x,splineCoordsTrans(i).x(end:-1:1),splineCoordsEdge(i).x(1)]; % Close polygon
    
    % Update min and max coordinates for display
    if showPolygons==1
        mny=min(y);
        mxy=max(y);
        mnx=min(x);
        mxx=max(x);
        if isempty(MNY)
            MNY=mny; MXY=mxy; MNX=mnx; MXX=mxx;
        end
        if mny<MNY, MNY=mny; end
        if mxy>MXY, MXY=mxy; end
        if mnx<MNX, MNX=mnx; end
        if mxy>MXY, MXY=mxy; end
    end   
    
    % Store current polygon
    lpPolygons(i).y=y;
    lpPolygons(i).x=x;
    
end

% If requested by the user, plot polygons
if showPolygons==1
    fsmTransPlotLpPolygons(lpPolygons,[MXY MXX])
end
