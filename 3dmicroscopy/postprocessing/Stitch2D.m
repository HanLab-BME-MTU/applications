%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       small program to stitch axial MIPs 
%	manual selection of shift vector
%       primitive version with visual stitching for light-sheet imaging.
%       by Reto Fiolka  12-04-2015
%                    
%                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all
dataP1='/project/cellbiology/gdanuser/december/shared/EB3mNG/151203/Cell4'

savePath='/home2/rfiolka/Pictures/test';

zdim=67;ydim=128;xdim=512;
V1=zeros(zdim,ydim,xdim);V2=zeros(zdim,ydim,xdim);V2=zeros(zdim,ydim,xdim);
for k=10:200
for i=1:zdim; %%
    if k<100
   V1(i,:,:)=double(imread([dataP1 '/ch1/1_CAM01_0000' num2str(k) '.tif'],i));
    V2(i,:,:)=double(imread([dataP1 '/ch2/1_CAM02_0000' num2str(k) '.tif'],i));
     V3(i,:,:)=double(imread([dataP1 '/ch3/1_CAM03_0000' num2str(k) '.tif'],i));
    else 
          V1(i,:,:)=double(imread([dataP1 '/ch1/1_CAM01_000' num2str(k) '.tif'],i));
    V2(i,:,:)=double(imread([dataP1 '/ch2/1_CAM02_000' num2str(k) '.tif'],i));
     V3(i,:,:)=double(imread([dataP1 '/ch3/1_CAM03_000' num2str(k) '.tif'],i));
    end
        
end

%%
Mip1=squeeze(sum(V1,2));Mip2=squeeze(sum(V2,2));Mip3=squeeze(sum(V3,2));
Mip1=[zeros(2*zdim,xdim);Mip1];Mip2=[zeros(zdim,xdim);Mip2;zeros(zdim,xdim)];Mip3=[Mip3;zeros(2*zdim,xdim)];


Mip1=Mip1-min(Mip1(:));Mip2=Mip2-min(Mip2(:));Mip3=Mip3-min(Mip3(:));

% manually set shift vectors
zshift1=-7; xshift1=31;
zshift3=11; xshift3=-7;


Mip1s=circshift(Mip1,[zshift1,xshift1]);
scale1=mean2(Mip2(2*zdim+1+zshift1:2*zdim,:))/mean2(Mip1s(2*zdim+1+zshift1:2*zdim,:));
Mip1=Mip1*scale1;

Mip3s=circshift(Mip3,[zshift3,xshift3]);
scale3=mean2(Mip2(zdim+1:zdim+zshift3,:))/mean2(Mip3s(zdim+1:zdim+zshift3,:));
Mip3=Mip3*scale3;

MIP=circshift(Mip1,[zshift1,xshift1])+Mip2+circshift(Mip3,[zshift3,xshift3]);

MIP(zdim+1:zdim+zshift3,:)=MIP(zdim+1:zdim+zshift3,:)/2;
MIP(2*zdim+1+zshift1:2*zdim,:)=MIP(2*zdim+1+zshift1:2*zdim,:)/2;

%MIP(zdim+1,:)=Mip2(zdim+1,:);
MIP(2*zdim+1+zshift1,:)=Mip2(2*zdim+1+zshift1,:);

%figure;imshow(MIP(zshift3+1:end+zshift1,:),[])

imwrite(uint16(MIP(zshift3+1:end+zshift1,:)),[savePath '/test' num2str(k) '.tif'],'Compression','none');
end




