function pSpPos=fsmTrackPropSpecklePos(spPos,Mi,mode,corLen)
% fsmTrackPropSpecklePos propagates speckle positions based on an interpolated vector field
%
% SYNOPSIS      pSpPos=fsmTrackPropSpecklePos(spPos,Mi,mode)
%
% INPUT         spPos      : n speckle positions saved in an (nx2) matrix [y x]n
%               Mi         : interpolated vector field saved in a (mx4) matrix
%                            [y0 x0 y x]m, where [y0 x0]m are the coordinates of 
%                            the interpolation points (see vectorFieldInterp)
%               mode       : (optional) either 'forward' or 'backward',
%                            defines whether the positions have to be forward 
%                            (next frame)- or backward (previous
%                            frame)-propagated.
%                            Default is 'forward'.
%
% OUTPUT        pSpPos     : propagated speckle positions [y x]n
%
% DEPENDENCES   fsmTrackPropSpecklePos uses { }
%               fsmTrackPropSpecklePos is used by { } 
%
% Aaron Ponti, January 3rd, 2003

% Check input parameters
if nargin<2
    error('At least two input parameters expected');
end
if size(spPos,2)~=2
    error('spPos must be an (nx2) matrix');
end
if size(Mi,2)~=4
    error('Mi must be an (mx4) matrix');
end
if nargin>=3
    mode=upper(mode);
    if ~strcmp(mode,'BACKWARD')
        mode='FORWARD';
    end
else
    mode='FORWARD';
end
if nargin < 4
    corLen = Inf;
end

% Calculate vectors
v=[Mi(:,3)-Mi(:,1) Mi(:,4)-Mi(:,2)];
if strcmp(mode,'BACKWARD')
    v=-v;
end

% Calculate all distances between speckles and interpolation points
D=createDistanceMatrix(spPos,Mi(:,1:2));

% Initialize output pSpPos
pSpPos=zeros(size(spPos));

% Calculate new speckle positions
for i=1:size(D,1)
    
    % For each speckle find the closest interpolation point
    minD = min(D(i,:));
    if minD <= corLen
        pos=find(D(i,:)==minD);
        pos=pos(1);

        % Update speckle position by adding the vector at the closest interpolation point
        pSpPos(i,1:2)=spPos(i,1:2)+v(pos,:);
    else %No propogation.
        pSpPos(i,1:2)= spPos(i,1:2);
    end

end

