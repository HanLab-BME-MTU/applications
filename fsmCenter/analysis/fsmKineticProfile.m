function [scoreProfile,polyProfile,depolyProfile,P,posScores,img2C]=fsmKineticProfile(firstKinScore,imgSize,n,lineDescr,dist,d0,sampling)
% function [scoreProfile,P,posScores]=fsmKineticProfile(firstKinScore,imgSize,n,lineDescr,dist,d0,sampling)
% fsmKineticProfile calculates a kinetic profile along a line from polyMap and depolyMap
%
% polyMap and depolyMap are returned by fsmTransKineticMaps
%
% [scoreProfile,polyProfile,depolyProfile,P,posScores]=fsmKineticProfile(firstKinScore,M,n,lineDescr,dist,d0,sampling)
%
% INPUT
%
% firstKinScore : string containing name (with full path) of the first kinScore.mat structure
%                 or a net assembly map (img2C) as returned by fsmTransKineticMaps
% imgSize       : size of the analyzed image
% M             : M stack as returned by fsmTrackMain
% n             : number of frames for averaging
%               : if n is set to 0, then all frames are considered
% lineDescr     : start and end points for the profile (line going from A to B)
%                           [ A(1) B(1)
%                             A(2) B(2) ]
% dist          : defines the distance from the profile where the contribution
%                 of vectors is neglected
% d0            : (pixels) correlation length for interpolation
% sampling      : sampling rate, defined as: sampling*(1/vector length) 
%               defines the distance between two points calculated
%               along the profile defined by lineDescr
%               default: 1
%
% OUTPUT
% scoreProfile  : net kinetic profile
% polyProfile   : net poly profile
% depolyProfile : net depoly profile
% P             : coordinates of the points along the profile
% posScores     : coordinates of the selected scores
% img2C         : net assembly map (see fsmTransKineticMaps)
%
% Aaron Ponti, 12/08/2003

if nargin==6
    sampling=1;
end

% Extract pixel coordinates along line
P=lineCoords(lineDescr(:,1),lineDescr(:,2),sampling);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create poly and depoly maps
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ischar(firstKinScore) % String pointing to the fist kinScore###.mat
    [poly,depoly,img2C]=fsmTransKineticMaps(firstKinScore,imgSize,[n n],0); %5
else
    if size(firstKinScore,3)~=3
        error('Please make sure that the net assembly map has the correct dimensions - see fsmTransKineticMaps');
    else
        poly=firstKinScore(:,:,1);
        depoly=firstKinScore(:,:,2);
        img2C=firstKinScore;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extract only scores within a given distance from any of 
%   the points along the profile
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Poly scores
[y x pScores]=find(poly);
D=createSparseDistanceMatrix([y x],P,dist);

% Extract vectors from list
E=sum(D,2);
indx=find(E==0); % Those vectors which are not close to any
                 % of the profile points have "virtual" distance 0
                 % in the sparse distance matrix
y(indx,:)=[];    % Kill these entries
x(indx,:)=[];
pScores(indx,:)=[];

POLY=[y x pScores];
clear y x D E indx pScores;

% Depoly scores
[y x dScores]=find(depoly);
D=createSparseDistanceMatrix([y x],P,dist);

% Extract vectors from list
E=sum(D,2);
indx=find(E==0); % Those vectors which are not close to any
                 % of the profile points have "virtual" distance 0
                 % in the sparse distance matrix
y(indx,:)=[];    % Kill these entries
x(indx,:)=[];
dScores(indx,:)=[];

DEPOLY=[y x dScores];
clear y x D E indx dScores;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Interpolate extracted POLY scores onto profile coordinates
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate distances
D=createDistanceMatrix(P,POLY(:,1:2));

% Correlation matrix (d0 may be a scalar or a matrix)
G=exp(-D.^2./d0.^2); clear D;

% Interpolate
polyProfile=G*POLY(:,3);

% Filter
sG=sum(G,2);
sG(find(sG==0))=1; % Prevent division by zero
polyProfile=polyProfile./sG;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Interpolate extracted DEPOLY scores onto profile coordinates
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate distances
D=createDistanceMatrix(P,DEPOLY(:,1:2));

% Correlation matrix (d0 may be a scalar or a matrix)
G=exp(-D.^2./d0.^2); clear D;

% Interpolate
depolyProfile=G*DEPOLY(:,3);

% Filter
sG=sum(G,2);
sG(find(sG==0))=1; % Prevent division by zero
depolyProfile=depolyProfile./sG;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Interpolate extracted scores onto profile coordinates
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% All scores
posScores=cat(1,POLY(:,1:2),DEPOLY(:,1:2));
scores=cat(1,POLY(:,3),DEPOLY(:,3));

% Calculate distances
D=createDistanceMatrix(P,posScores);

% Correlation matrix (d0 may be a scalar or a matrix)
G=exp(-D.^2./d0.^2); clear D;

% Filter
scoreProfile=G*scores;

% Normalize
sG=sum(G,2);
sG(find(sG==0))=1; % Prevent division by zero
scoreProfile=scoreProfile./sG;

% figure
% plot(polyProfile,'r-');
% hold on
% plot(depolyProfile,'b-');


