function C=candVectorsCirc(Lmax1,Lmax2,radius)
% candVectorsCirc fills the graph tracker input matrix with candidates from a circular region of given radius
%
% SYNOPSIS      C=candVectorsCirc(Lmax1,Lmax2,radius)
%
% INPUT         Lmax1  : local maximum map of frame 1
%               Lmax2  : local maximum map of frame 2
%                        (with original intensities at the local maximum position)
%               radius : radius of the circular region around the
%                        source speckle to look for candidate targets
% OUTPUT        C      : input matrix for the two graph-based trackers
%
% DEPENDENCES   candVectorsCirc uses { createDistanceMatrix , missingIndices }
%               candVectorsCirc is used by { } 
%
% Aaron Ponti, March 11th, 2003

% Find speckle positions
[y x]=find(Lmax1~=0); posI=[y x];
[y x]=find(Lmax2~=0); posJ=[y x];

% Calculate distance matrix -- to be replaced with createSparseDistanceMatrix!
D=createSparseDistanceMatrix(posI,posJ,radius);

% Fill candidate vectors matrix C
c0=1;
C = []; % inicialization
for c1=1:size(D,1)
    n=find(D(c1,:)>0);
    ln=length(n);
    count=0;
    for c2=c0:c0+ln-1
        count=count+1;
        C(c2,1:6)=[posI(c1,:),posJ(n(count),:),Lmax1(posI(c1,1),posI(c1,2)),Lmax2(posJ(n(count),1),posJ(n(count),2))];
    end
    c0=c0+ln;
end

% Add non-paired speckles
[y x]=find(D>0);
m=missingIndices(y,size(D,1));
if ~isempty(m)
    sz=size(C,1);
    for c1=1:length(m)
        C(sz+c1,1:2)=posI(m(c1),:);
        C(sz+c1,5)=Lmax1(posI(m(c1),1),posI(m(c1),2));
    end
end
n=missingIndices(x,size(D,2));
if ~isempty(n)
    sz=size(C,1);
    for c1=1:length(n)
        C(sz+c1,3:4)=posJ(n(c1),:);
        C(sz+c1,6)=Lmax2(posJ(n(c1),1),posJ(n(c1),2));
    end
end
