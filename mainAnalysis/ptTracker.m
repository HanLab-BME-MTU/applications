function M=ptTracker (I,J,radius,influence)
% ptTracker tracks particle postions resolving topological conflicts with neural networks
%
% SYNOPSIS   M = ptTracker (Lmax1,Lmax2,radius,influence)
%
% INPUT      I         :  matrix [y x I]n of speckle coordinates and intensities for frame 1
%            J         :  matrix [y x I]n of speckle coordinates and intensities for frame 2
%            radius    :  search radius
%            influence :  radius of influence for the neural network tracker. This is the initial search
%                         radius; the neural network can more effectively reconstruct flow if it can
%                         search over larger areas and therefore take more particles into account. The final
%                         matches will be constrained to have a maximum distance = 'radius', but initially 
%                         a larger search radius ('influence') will be used.
%
% OUTPUT     M         :  matrix of matches [y1(1) x1(1) y1(2) x1(2)]n 
%
% DEPENDENCES   ptTracker uses { createSparseDistanceMatrix ;
%                                           missingIndices ;
%                                           repeatingIndices ;
%                                           labonteNN }
%               ptTracker is used by { ptTrackCells }
%
% Aaron Ponti, March 14th, 2003
% Andre Kerstens, June 2004 - copied from FSM project

% Check input parameters
if nargin~=4
    error('Four input parameters expected');
end

% Initialize M
M=[];

% Options for LabonteNN (Neural Network)
options.r         = influence;
options.rf        = 0.005; % 1e-4; % Original value: 1e-5
options.alpha     = 5e-4; %5e-4; % Original value: 5e-4
options.beta      = 0.99; % Original value: 0.99
options.threshold = radius;
options.maxiter   = 1000;

% Extract speckle coordinates
posI=I(:,1:2);
posJ=J(:,1:2);

% Calculate distance map
D=createSparseDistanceMatrix(posI,posJ,radius);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extract non-paired speckles
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[y x]=find(D);
m1=missingIndices(y,size(D,1));
m2=missingIndices(x,size(D,2));
nonPairedI=posI(m1,:);      % List of coordinates of non-paired speckles from frame 1
nonPairedJ=posJ(m2,:);      % List of coordinates of non-paired speckles from frame 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extract speckles with topological conflicts
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[y x]=find(D);
matches=[y x];
[r1,n1]=repeatingIndices(y);
[r2,n2]=repeatingIndices(x);

% Store speckles with more than one partner
oneMatch=cat(2,r1,r2);
conflictMatches=matches(oneMatch,:);  % Speckles in conflict

% Store speckles univocally assigned
matches(oneMatch,:)=[];               % Speckles univocally assigned

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Add univocally assigned speckles to M
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty (matches)
   for c1=1:size(matches,1)
      M(c1,1:4)=[posI(matches(c1,1),:) posJ(matches(c1,2),:)];
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Solve topological conflicts with neural networks
%   and add the speckles to M
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(conflictMatches)
    nnPosI=posI(unique(conflictMatches(:,1)),:);
    nnPosJ=posJ(unique(conflictMatches(:,2)),:);
    conflict=1;
    while ~isempty(conflict)
        [Mnn,conflict,options]=labonteNN(nnPosI,nnPosJ,options);
        M=cat(1,M,Mnn);
        if ~isempty(conflict)
            disp('Some conflicts need to be resolved in a second run of the Neural Network tracker.'); 
            nnPosI=unique(conflict(:,1:2),'rows');
            nnPosJ=unique(conflict(:,3:4),'rows');
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Add non-paired speckles to M
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add non-paired speckles from frame 1
lenNPI=size(nonPairedI,1);
MNPI=zeros(lenNPI,4);
MNPI(1:lenNPI,1:2)=nonPairedI;
M=cat(1,M,MNPI);

% Add non-paired speckles from frame 2
lenNPJ=size(nonPairedJ,1);
MNPJ=zeros(lenNPJ,4);
MNPJ(1:lenNPJ,3:4)=nonPairedJ;
M=cat(1,M,MNPJ);
