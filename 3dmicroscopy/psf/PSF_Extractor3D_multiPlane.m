%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       PSF_Extractor
%
%       This program extracts the 2DPSF-Profiles 
%       user selection. Maximum Position is recalibrated and for multiple 
%       Profiles they are averaged. It averages as many PSF as you select
%
%       
%
%       by Reto Fiolka 08 07 2006
%                      07.09.2007 added automatic FWHM measurement
%                      (via interpoaltion) in pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function [AvPSF]=PSF_Extractor3D(V)
%%

clear all 
clc
close all
pixX=160;
pixZ=160;

ydim=512; xdim=512;zdim=21;

%ydim=xdim;

%dataP='/project/cellbiology/gdanuser/shared/LouisXIVdata';
dataP='/project/cellbiology/gdanuser/shared/MultiPlane/TestDefocus/Beads10BmicronNormalLowNA/150717/Cell4'
V=zeros(zdim,ydim,xdim);
for i=1:zdim; 
   V(i,:,:)=double(imread([dataP '/1_CAM01_000000.tif'],i));
    %V(i,:,:)=double(imread([dataP '/deconEDFpsf160b4.tif'],i));
end
%%
%V=V(150:end-150,200:end-200,:);

%%
clear FWHMx FWHMy FWHMz X2 X
[zdim,ydim2,x2]=size(V);
A=squeeze(max(V));
layer=0;
figure('Name', ['Stack Layer' num2str(layer) ],'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);

imshow(A,[],'Border','tight');iptsetpref('ImshowAxesVisible','off'); 

[x,y]=getpts;

x=round(x);
y=round(y);

%Yav=(y(1)+y(2))/2;

counter=(length(x));

temp=zeros(32,32);
AvPSF=zeros(32,32);

count=1;
for i=1:counter
   

SubWin=A(y(i)-3:y(i)+3,x(i)-3:x(i)+3);
[maxi Indi]=max(SubWin(:));
[ypeak,xpeak]= ind2sub (size(SubWin),Indi);

xpeak=round(x(i)-4+xpeak); 
ypeak=round(y(i)-4+ypeak);

[maxi2,zpeak]=max(V(:,ypeak,xpeak));


lineX=squeeze(V(zpeak,ypeak,xpeak-10:xpeak+10));
lineY=squeeze(V(zpeak,ypeak-10:ypeak+10,xpeak));


lineX=lineX-min(lineX); lineX=lineX/max(lineX);
lineY=lineY-min(lineY); lineY=lineY/max(lineY);

param0=[1,2,10,0];


%figure;plot(gaussianFunc(param,[1:21]));hold on;plot(lineZ,'r')
[param,resnorm,residual] = lsqcurvefit(@gaussianFunc,param0,[1:21],lineX');
FWHMx(count)=param(2)*2.3548*pixX;

[param,resnorm,residual] = lsqcurvefit(@gaussianFunc,param0,[1:21],lineY);
FWHMy(count)=param(2)*2.3548*pixX;
PeakI(count)=maxi;

X(count)=xpeak;
Y(count)=ypeak;
%figure;plot(gaussianFunc(param,[1:21]));hold on;plot(lineZ,'r')
count=count+1;
end



FWHMx
FWHMy
maxi



%%
%[X2,sortind]=sort(X,'ascend');

%figure;plot(X2,FWHMz(sortind)*pixZ,'*');hold on;plot(X2,FWHMy(sortind)*pixX,'r*') ;%hold on;plot(X2,FWHMx(sortind)*160,'k*')
%%
% figure;surf(AvPSF);
% %
% figure;plot(AvPSF(:,15));
% figure;plot(AvPSF(15,:));
% 
% psf=AvPSF(:,15);
% 
% X=[1:length(psf)];
% XI=[1:0.1:length(psf)];
% psf=interp1(X,psf,XI);
% psf2=abs((psf)/max(psf)-0.5);
% psf2=1./psf2;
% 
% indmax= find(imregionalmax(psf2) == 1);        %indexes of maxima
% [peaks,sortind]= sort(psf2(indmax),2,'descend');
% FWHM=abs(indmax(sortind(1))-indmax(sortind(2)))/10
%  
% psf=AvPSF(15,:);
% 
% X=[1:length(psf)];
% XI=[1:0.1:length(psf)];
% psf=interp1(X,psf,XI);
% psf2=abs((psf)/max(psf)-0.5);
% psf2=1./psf2;
% 
% indmax= find(imregionalmax(psf2) == 1);        %indexes of maxima
% [peaks,sortind]= sort(psf2(indmax),2,'descend');
% FWHM2=abs(indmax(sortind(1))-indmax(sortind(2)))/10





