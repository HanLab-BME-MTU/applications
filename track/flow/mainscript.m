%script demonstrating pascal's tracker
%define a maximum possible jump for speckles between two frames
global param
param.searchdist=4
%coonstruct candidate vectors from two images containing intensities of speckles at the positions of these speckles and zeroes otherwise.
asifirst=candvectorsIijwithint('C:\User\vallotton\matlab\data\fsm\retro1\ACTINMOVIEFEAT1.tif','C:\User\vallotton\matlab\data\fsm\retro1\ACTINMOVIEFEAT2.tif');
%save the notarget no source part
indexnotarget=find(asifirst(:,5)==0);
indexnosource=find(asifirst(:,3)==0);
%remove the first column which contains the distances between speckles
asifirst=asifirst(:,2:7); %chop off the distance column which is not used by the algo.
%call the magic function
%asifirst=asifirst(1:indexnotarget-1,:)
myasiout1=track1063(asifirst,asifirst);
%reorganize the output of the magic function
l1=prod(size(myasiout1))-1;
s1=1;
s2=1;
a=0;
for (i1=0:l1)
    a(1+floor(i1/4),1+rem(i1,4))=myasiout1(i1+1);
end
myasiout1=a;
%kill the entries you don't really like :)
myasiout1=(myasiout1(find (myasiout1(:,4)>0),:));
myasiout1=(myasiout1(find (myasiout1(:,3)>0),:));
myasiout1=(myasiout1(find (myasiout1(:,2)>0),:));
myasiout1=(myasiout1(find (myasiout1(:,1)>0),:));
myasiout1=[myasiout1(:,2),myasiout1(:,1),myasiout1(:,4),myasiout1(:,3)];


img1=imread('C:\User\vallotton\matlab\data\fsm\retro1\ACTINMOVIEFEAT1.tif');
img2=imread('C:\User\vallotton\matlab\data\fsm\retro1\ACTINMOVIEFEAT2.tif');
img1=double(img1);
img2=double(img2);
I=find(ne(img1,0));
J=find(ne(img2,0));
[sizey,sizex]=size(img1);

temp1=[]
temp2=[]
for i=transpose(I) %loop over all speckles in source image
    %[yJ,xJ]=indexij(sizey,sizex,I); %transfert out of the loop!
    [yi,xi]=indexij(sizey,sizex,i);
    if(size((find(myasiout1(:,1)==yi & myasiout1(:,2)==xi)),1)==0)
    temp1=[temp1;[0,yi,xi,0,0,0,0]];
    end
end
J=find(ne(img2,0));

for k=transpose(J) 
    %[yJ,xJ]=indexij(sizey,sizex,J); %transfert out of the loop!
    [yi,xi]=indexij(sizey,sizex,k);
    if(size((find(myasiout1(:,3)==yi & myasiout1(:,4)==xi)),1)==0)
    temp2=[temp2;[0,0,0,yi,xi,0,0]];
    end
end






%stick back
M=[myasiout1;asifirst(indexnotarget,1:4);temp1(:,2:5);asifirst(indexnosource,1:4);temp2(:,2:5)]