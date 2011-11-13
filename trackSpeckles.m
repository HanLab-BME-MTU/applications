function [M vectors]=trackSpeckles(I,J,threshold, varargin)
% trackSpeckles uses the interpolated vector field to refine the tracking
%
% SYNOPSIS   M = trackSpeckles(I,J,threshold)
%            M = trackSpeckles(I,J,threshold,'initM',initM,'initCorLen',initCorLen)
%
% INPUT      I          : matrix [y x I]n of speckle coordinates and intensities for frame 1
%            J          : matrix [y x I]n of speckle coordinates and intensities for frame 2
%            threshold  : radius of the region searched by the tracker for matching speckles
%            gridSize   : (optional, default = 0) distance between two interpolation points on the grid [gy gx]
%                         if gridSize is set equal to zero, the field is interpolated onto
%                         the original vector positions
%
% OUTPUT     M          : matrix of matches [y1(1) x1(1) y1(2) x1(2)]n 
%
% References:
% A. Ponti et al., Biophysical J., 84 336-3352, 2003.
% A. Ponti et al., Biophysical J., 89 2459-3469, 2005.

% Aaron Ponti, September 8th, 2004
% Sebastien Besson, 5/2011
% Copied from fsmTrackTrackerMain

% Input Check
ip =inputParser;
ip.addRequired('I',@(x) isnumeric(x) && size(x,2)==3);
ip.addRequired('J',@(x) isnumeric(x) && size(x,2)==3);
ip.addRequired('threshold',@isscalar);
ip.addParamValue('initM',[],@isnumeric);
ip.addParamValue('initCorLen',Inf,@isscalar);
ip.addParamValue('enhanced',0,@isscalar);
ip.addParamValue('corrLength',0,@isscalar);
ip.parse(I,J,threshold, varargin{:})
initM = ip.Results.initM;
initCorLen = ip.Results.initCorLen;
enhanced = ip.Results.enhanced;
corrLength = ip.Results.corrLength;

% Track with initialization matrix initM as an initializer
M=fsmTrackTrackerIterative(initM,I,J,threshold,initCorLen);

% Returns if no hierarchical tracking
if ~enhanced, return; end

% Extract vector field from M (discard non-matched speckles)
raw=M(M(:,1)~=0 & M(:,3)~=0,:);

% In the very unlikely case that ALL particles are unmatched
if isempty(raw), vectors=zeros(1,4); return; end

% Extract vectors from M
raw=M(M(:,1)~=0 & M(:,3)~=0,:);
% Interpolate onto vector positions
grid=raw(:,1:2);

% Average returned M to be used to propagate I again
vectors=vectorFieldAdaptInterp(raw,grid,corrLength,[],'strain');

% Track with vectors as an initializer
M=fsmTrackTrackerIterative(vectors,I,J,threshold,initCorLen);


function M=fsmTrackTrackerIterative(initM,I,J,threshold,initCorLen)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If an initializer for the tracker exists, use it to propagate speckle positions at frame I
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(initM)
    
    % Sort I and make sure there are no empty lines
    I=sortrows(I(I(:,1)~=0,:));
    
    % Propagate speckle positions
    spPos=I(:,1:2); % Extract only positions
    pSpPos=propagateSpecklePositions(spPos,initM,'FORWARD',initCorLen);
    
    % Now I is the propagated version of the original I
    I=cat(2,pSpPos,I(:,3)); % Add original intensities to the propagated positions
    clear Itmp;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Track with linear assignment code
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NONLINK_MARKER = -1;    % value to indicate that two points cannot be linked.
extendedTesting = 0;
augmentCC = 1;          % automatically generate ghost marker to account for birth and deaths


% Extract speckle coordinates
posI(:,1)=I(:,1);
posI(:,2)=I(:,2);
posJ(:,1)=J(:,1);
posJ(:,2)=J(:,2);
% create cost matrix for the linear assignment code
% in the most general case that is equal to the distance matrix
cc = createSparseDistanceMatrix(posI, posJ, threshold);
%cc = createDistanceMatrix(posI, posJ);
[all_links_1, all_links_2] = lap(cc, NONLINK_MARKER, extendedTesting, augmentCC);

I_N = length(posI);
J_N = length(posJ);

% first get points that are in frame 1 and have a link in frame 2
all_links_1_2 = all_links_1(1:I_N);
frame_1_indices  = find(all_links_1_2 <= J_N);
frame_2_indices  = all_links_1(frame_1_indices);
Mone = [posI(frame_1_indices,1:2),  posJ(frame_2_indices,1:2)];

% Add non-paired speckles from frame 1
frame_1_indices_noLink  = find(all_links_1_2 > J_N);
lenNPI=size(frame_1_indices_noLink,1);
MNPI=zeros(lenNPI,4);
MNPI(1:lenNPI,1:2)=posI(frame_1_indices_noLink,1:2);
Mone=cat(1,Mone,MNPI);

% Add non-paired speckles from frame 2
all_links_2_1 = all_links_2(1:J_N);
frame_2_indices_noLink  = find(all_links_2_1 > I_N);
lenNPJ=size(frame_2_indices_noLink,1);
MNPJ=zeros(lenNPJ,4);
MNPJ(1:lenNPJ,3:4)=posJ(frame_2_indices_noLink,1:2);
Mone=cat(1,Mone,MNPJ);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If needed, correct back the propagated coordinates of the first frame
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(initM)

    % Create a copy of Mone
    copyMone=Mone;
    
    % Correct back the coordinates of the first frame (= currentM(:,1:2))
    for i=1:size(spPos,1)
        indx=find(Mone(:,1)==pSpPos(i,1) & Mone(:,2)==pSpPos(i,2));
        if length(indx)~=1
            error('Propagated position not univocally found in Mone');
        end
        copyMone(indx,1:2)=spPos(i,:);
    end
    M=copyMone;

else
    
    % Return the result of the first tracking
    M=Mone;
    
end