function [velProfile,P,Mi,list]=fsmVelocityProfile(M,n,lineDescr,dist,d0,sampling)
% fsmVelocityProfile calculates a velocity profile along a line from M
%
% M is the stack of vector matches returned by fsmTrackMain
%
% [vel,P,Mi,list]=fsmVelocityProfile(M,n,lineDescr,dist,d0,sampling)
%
% INPUT
%
% M          : M stack as returned by fsmTrackMain
% n          : number of frames
%            : if n is set to 0, then all frames are considered
% lineDescr  : start and end points for the profile (line going from A to B)
%                          [ A(1) B(1)
%                            A(2) B(2) ]
% dist       : defines the distance from the profile where the contribution
%              of vectors is neglected
% d0         : (pixels) correlation length for interpolation
% sampling   : sampling rate, defined as: sampling*(1/vector length) 
%              defines the distance between two points calculated
%              along the profile defined by lineDescr
%              default: 1
%
% OUTPUT
% velProfile : velocity profile
% P          : coordinates of the points along the profile
% Mi         : interpolated vectors along the profile
% list       : vectors used to interpolate along the profile
%
% Aaron Ponti, 12/08/2003

if nargin==5
    sampling=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Collect all vectors into a long list
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=size(M,1);
if n==0
    d=size(M,3);
else
    d=n;
end

% Check d
sizeM=size(M,3);
if d>sizeM
    d=sizeM;
    fprintf(1,'Warning: there are not enough frames in M to process the requested %d frames. Using %d frames.\n',n,d);
end

% Initialize memory for the list
list=zeros(h*d,4);

% Fill list
for i=1:d
    list((i-1)*h+1:i*h,1:4)=M(:,:,i);
end

% Remove all zeros entries
list=list(find(list(:,1)~=0 & list(:,3)~=0),:);

% Extract pixel coordinates along line
P=lineCoords(lineDescr(:,1),lineDescr(:,2),sampling);

% Consider only vectors which are at most dist away from any of the profile points
D=createSparseDistanceMatrix(list(:,1:2),P,dist);

% Extract vectors from list
E=sum(D,2);
indx=find(E==0); % Those vectors which are not close to any
                 % of the profile points have "virtual" distance 0
                 % in the sparse distance matrix
list(indx,:)=[]; % Kill these entries

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Interpolate extracted vectors onto profile coordinates
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the divergence of the vector field along the profile line
% [div,d0]=vectorFieldDiv(list,P,d0,[]);

% Update d0 depending on divergence
% d0=updateD0FromDiv(div,d0,1,size(list,1),size(list,1));

% Interpolate along the profile line
Mi=vectorFieldInterp(list,P,d0,[]);

% Calculate velocity profile
vectors=[Mi(:,3)-Mi(:,1) Mi(:,4)-Mi(:,2)];
velProfile=sqrt(vectors(:,1).^2+vectors(:,2).^2);
