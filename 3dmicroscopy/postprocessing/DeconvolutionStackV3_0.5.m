%%
clear all 
clc
close all
render=0;
saving=1;
agg=-1;      % aggressivness: -1: gentle, 0: smooth, 1: sharp, 2: sharp but noisy
crop=0;

xdim=512;zdim=151;

dataP='/project/cellbiology/gdanuser/melanomaModel/PrimaryMelanoma/M405/TwoPhotonLightSheet/GFP/140219/Cell5/';
dataS='/project/cellbiology/gdanuser/shared/Deconvolution2';
load('/project/cellbiology/gdanuser/shared/AveragedPSF/PSFAVE501slicesZ1024compress.mat');

mkdir(dataS)

%%
PSF=PSFave(round((501-zdim)/2):round((501-zdim)/2)+zdim-1,round((1024-xdim)/2):round((1024-xdim)/2)-1+xdim,round((1024-xdim)/2):round((1024-xdim)/2)-1+xdim);

%%
%clear PSFave

OTFrot=fftshift(fftn(PSF));

OTFrot=abs(OTFrot);
figure;imagesc(squeeze(OTFrot(:,255,:)))

%%

for n=1:81



cell=zeros(zdim,xdim,xdim);
for i=1:zdim; 
    temp(i,:,:)=im2double(imread([dataP num2str(n) '.tif'],i));
end

if n==1
    if crop
        [z,x,y]=size(temp);
        for k=1:x
            for l=1:y
                mip(k,l)=sum(temp(:,k,l));
            end
        end
        
        
        figure;imshow(mip,[])
        disp('please click on center of cell')
        [x1,y1]=ginput(1);
        
        if x1+round(xdim/2)>x
            cx=x-(round(xdim/2)+1);
        elseif x1-round(xdim/2)<0
            cx=round(xdim/2)+1;
        else
            cx=x1;
        end
        x1+round(xdim/2)>x
        if y1+round(xdim/2)>y
            cy=y-(round(xdim/2)+1);
        elseif y1-round(xdim/2)<0
            cy=round(xdim/2)+1;
        else
            cy=y1;
        end
        cx=round(cx);cy=round(cy);
        
        
        cell=temp(:,cy-round(xdim/2):cy+round(xdim/2)-1,cx-round(xdim/2):cx+round(xdim/2)-1);
        
    else
        
        cell=temp;
    end
    
else
    if crop
        cell=temp(:,cy-round(xdim/2):cy+round(xdim/2)-1,cx-round(xdim/2):cx+round(xdim/2)-1);
    else
        cell=temp;
    end
end

disp('done loading data')
  %%  Apodizationan mask

  if n==1
      
      elli=zeros(zdim,xdim,xdim);
      
      if agg==-1
          a=round(xdim/6);b=round(xdim/6);c=round(zdim/6.4);
      elseif agg==0
          a=round(xdim/5.1);b=round(xdim/5.1);c=round(zdim/5.3);
      elseif agg==1
          a=round(xdim/4.65);b=round(xdim/4.65);c=round(zdim/4.7);
      elseif agg==2
          a=round(xdim/3.5);b=round(xdim/3.5);c=round(zdim/3.8);
          disp('aggressive')
      end
      disp('start creating apodization mask')
      
      for i=1:zdim;for k=1:xdim;for l=1:xdim;
                  z=i-floor(zdim/2)+1;
                  x=k-floor(xdim/2)+1;
                  y=l-floor(xdim/2)+1;
                  
                  r=x^2/(a*a)+y^2/(b*b)+z^2/(c*c);
                  
                  if r<1
                      elli(i,k,l)=1;
                  end
                  
                  %a=round(xdim/5.1);b=round(xdim/5.1);c=round(zdim/5.3);
                  
              end
          end
      end
      
      [u,v,w] = meshgrid(-round(xdim/2):round(xdim/2)-1,-floor(zdim/2):round(zdim/2)-1,-round(xdim/2):round(xdim/2)-1);
      
      Gauss =exp(-((1/3)^2*v.^2 + (1/3)^2*u.^2+ (1/3)^2*w.^2 )); Gauss=Gauss/max(Gauss(:));
      
      elli2=fftshift(fftn(elli));
      mask=Gauss.*elli2;
      mask=ifftn(ifftshift(mask));
      disp('done with apodization mask')
  end

%% Wiener Deconvolution
% disp('0')
cellMax = max(abs(cell(:)));
w=0.00001;
cell2=fftshift(fftn(cell));

% disp('1')
% max(abs(cell2(:)))
cell2=cell2.*OTFrot./((OTFrot.*OTFrot)+w);

% disp('2')
% max(abs(cell2(:)))
cell2=cell2.*mask;

% disp('3')
% max(abs(cell2(:)))
cell2=ifftn(ifftshift(cell2));

% disp('4')
% max(abs(cell2(:)))
%cell2 = cellMax*cell2/max(cell2(:));
cell2 = cell2*(10^7);
max(abs(cell2(:)))

% disp('5')
% max(abs(cell2(:)))

%% Maximum Intensity projection
disp('rendering data')

if render
    
    clear mip
    for i=1:zdim;for k=1:xdim;
            
            mip(i,k)=max(cell(i,:,k));
        end
        
    end
    
    clear mip4
    for i=1:xdim;for k=1:xdim;
            
            mip4(i,k)=max(abs(cell(:,i,k)));
        end
    end
    
    figure;imshow(imresize(mip,2),[]);%colormap jet
    figure;imshow(imresize(mip4,2),[]);%colormap jet
    clear mip2
    
    for i=1:zdim;for k=1:xdim;
            
            mip2(i,k)=max(abs(cell2(i,:,k)));
        end
    end
    cell2=abs(cell2);
    %cell2=cell2/max(cell2(:))*2^16;
    cell2 = cell2*(2^16-1);
    clear mip3
    for i=1:xdim;for k=1:xdim;
            
            mip3(i,k)=max(abs(cell2(:,i,k)));
        end
    end
    
    figure;imshow(imresize(mip2,2),[]);%colormap jet
    figure;imshow(imresize(mip3,2),[]);%colormap jet
    

    
end

%% saving data

cell2=abs(cell2);
%cell2=cell2-min(cell2(:));
%cell2=cell2/max(cell2(:))*20000;
cell2 = (2^16-1)*cell2;
%max(cell2(:))
if render==0
        clear mip3
    for i=1:xdim;
        for k=1:xdim;
            mip3(i,k)=max(abs(cell2(:,i,k)));
        end
    end
end

%mip3=abs(mip3); 
%mip3=mip3-min(mip3(:));
%mip3=mip3/max(mip3(:))*20000;
%mip3 = (2^16-1)*mip3;
%max(mip3(:))
if saving
    
for i=1:zdim
    if i == 1
        % the 5 is a hack
        imwrite(squeeze(uint16(5*cell2(i,:,:))),[dataS '/Deconvolved_' num2str(n-1) '.tif'],'Compression','none');
    else
        imwrite(squeeze(uint16(5*cell2(i,:,:))),[dataS '/Deconvolved_' num2str(n-1) '.tif'],'Compression','none','WriteMode','append');
    end
end

imwrite(squeeze(uint16(mip3)),[dataS '/DeconvolvedMIP_' num2str(n-1) '.tif'],'Compression','none');
end


end


