function candvectorsIijwithint=candvectorsIijwithint(img1,img2,ext)

%img1 source img containing points with intensities
%img2 target img
% global param;
% ext=param.searchdist %ext: half the lateral size of the patch containing candidate vectors
% 
% img1=imread(img1);
% img2=imread(img2);
% 
% img1=double(img1);
% img2=double(img2);
I=find(ne(img1,0));
J=find(ne(img2,0));

[sizey,sizex]=size(img1);

targindex=1; % index describing current position in the array notarg (for those speckles which have no target)
targsource=1; %index describing current position in the array nosource (for those speckles which have no source)
ct=1; %index for the current position in the assigned speckles array (those speckles which have a potential match)

%initialisation of an array to contain the result
candvectorsIijwithi(1,:)=[0.1,0.3,0.5,0.6,0.8,0.2,0.8];

for i=transpose(I) %loop over all speckles in source image
    
    [yJ,xJ]=indexij(sizey,sizex,J); %transfert out of the loop!
    [yi,xi]=indexij(sizey,sizex,i);
    
    
    
    filter=find((indexij1(sizey,sizex,J,1)>(yi-ext))&... % return those speckles in the target which are in a patch of lateral dim 2 ext centered on the current pixel i
        ( indexij1(sizey,sizex,J,1) <(yi+ext))&...
        (indexij1(sizey,sizex,J,2)>(xi-ext))&...
        (indexij1(sizey,sizex,J,2)<(xi+ext)));
        
    
    for j=transpose(filter) 
        [sourcey,sourcex]=indexij(sizey,sizex,i);
        intsourc=img1(sourcey,sourcex);
        [targety,targetx]=indexij(sizey,sizex,J(j));
        inttarg=img2(targety,targetx);
        candvectorsIijwithi(ct,:)=[((targetx-sourcex)^2+(targety-sourcey)^2)^0.5,sourcey,sourcex,targety,targetx,double(intsourc),double(inttarg)]; 
        ct=ct+1;
    end
end
op=1;
%Mplotarr(candvectorsIijwithi(:,2:5));

candvectorsIijwithint=candvectorsIijwithi;
N=candvectorsIijwithint;

for i=transpose(I) %loop over all speckles in source image
    %[yJ,xJ]=indexij(sizey,sizex,I); %transfert out of the loop!
    [yi,xi]=indexij(sizey,sizex,i);
    if(size((find(N(:,2)==yi & N(:,3)==xi)),1)==0)
    candvectorsIijwithint=[candvectorsIijwithint;[0,yi,xi,0,0,0,0]];
    end
end

for i=transpose(J) 
    %[yJ,xJ]=indexij(sizey,sizex,J); %transfert out of the loop!
    [yi,xi]=indexij(sizey,sizex,i);
    if(size((find(N(:,4)==yi & N(:,5)==xi)),1)==0)
    candvectorsIijwithint=[candvectorsIijwithint;[0,0,0,yi,xi,0,0]];
    end
end

crap=0;