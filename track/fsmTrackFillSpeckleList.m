function [specklePos,speckleMap]=fsmTrackFillSpeckleList(cands,imgSize)
% fsmTrackFillSpeckleList creates an (nx3) matrix [y x I]n of speckle coordinates and intensity and a 2D speckle map
%
% SYNOPSIS   [specklePos,speckleMap]=fsmTrackFillSpeckleListz(cands,imgSize)
%
% INPUT      cands     : structure with speckle information as saved in the PREPROCESSING MODULE
%            imgSize   : image size [y x] of the analyzed image
%
% OUTPUT     specklePos: matrix [y x I]n of speckle coordinates and intensities
%            speckleMap: 2D speckle map
%
% DEPENDENCES
%
% Aaron Ponti, October 8th, 2002

if nargin~=2
    error('Two parameters expected');
end

% Check cands fields
if (~isfield(cands,'status') | ...
        ~isfield(cands,'Lmax') | ...
        ~isfield(cands,'ILmax') | ...
        ~isfield(cands,'status'))
    error('The structure cands is not valid');
end

% Check imgSize
if size(imgSize)~=[1 2]
    error('imgSize must be a [1x2] row vector');
end

% Constant
nSp=length(cands);

% Initialize specklePos matrix with maximum possible dimensions
specklePos=zeros(nSp,3);

% Initialize empty (zero-filled) speckleMap
speckleMap=zeros(imgSize);

% Initialize speckle counter
numberOfSpeckles=0;

% Fill matrix
for counter1=1:nSp
    if cands(counter1).status==1
        numberOfSpeckles=numberOfSpeckles+1;
        specklePos(numberOfSpeckles,1:3)=[cands(counter1).Lmax, cands(counter1).ILmax];
        speckleMap(cands(counter1).Lmax(1),cands(counter1).Lmax(2))=cands(counter1).ILmax;
    end
end 

% Cut the not used rows of specklePos
if numberOfSpeckles~=nSp
    specklePos(numberOfSpeckles+1:nSp,:)=[];
end
