function M=fsmTrackTrackerP(I,J,threshold)
% fsmTrackTrackerP prepares the input and converts the output of track1063.dll
% into the Magic Position Matrix
%
% SYNOPSIS   M=fsmTrackTrackerP(I,J,threshold)
%
% INPUT      I          :   1st image loc-max Map
%            J          :   2nd image loc-max Map
%            threshold  :   radius of consideration for the tracker
%
% OUTPUT     M          :   the Magic Position Matrix
%
% DEPENDENCES   fsmTrackTrackerP uses { candvectorsCirc, track1063.dll, createDistanceMatrix}
%               fsmTrackTrackerP is used by { }

[sizey sizex]=size(I);
Ilin=find(I);
Jlin=find(J);

[y x]=find(I);posI=[y x];
[y x]=find(J);posJ=[y x];

% Candidate vectors
asifirst=candVectorsCirc(I,J,threshold);

% Call the magic function
clear track1063;
myasiout1=track1063(asifirst,asifirst);
% Reorganize the output of the magic function
l1=prod(size(myasiout1))-1;
s1=1;
s2=1;
a=0;
for (i1=0:l1)
    a(1+floor(i1/4),1+rem(i1,4))=myasiout1(i1+1);
end
myasiout1=a;

% Kill the entries which are not suitable
[y x]=find(myasiout1==0);
myasiout1(y,:)=[];
myasiout1=[myasiout1(:,2),myasiout1(:,1),myasiout1(:,4),myasiout1(:,3)];

% Convert linear index of nonzero positions in first image (I), resp. second image (J) to coordinate lists 
[y1l,x1l]=indexij(sizey,sizex,Ilin); 
[y2l,x2l]=indexij(sizey,sizex,Jlin); 
vect1=[y1l,x1l];
vect2=[y2l,x2l];

% Find non-matched entries in image 1
Imatched=myasiout1(:,1:2);
D1=createDistanceMatrix(Imatched,posI);
[y x]=find(D1==0);
% Remove matched entries from posI
posI(x,:)=[];

% Find non-matched entries in image 2
Jmatched=myasiout1(:,3:4);
D2=createDistanceMatrix(Jmatched,posJ);
[y x]=find(D2==0);
% Remove matched entries from posJ
posJ(x,:)=[];

% Fill M
lM=size(myasiout1,1);lI=size(posI,1);lJ=size(posJ,1);
totL=lM+lI+lJ;
M=zeros(totL,4);
M(1:lM,1:4)=myasiout1;
M(lM+1:lM+lI,1:2)=posI;
M(lM+lI+1:totL,3:4)=posJ;