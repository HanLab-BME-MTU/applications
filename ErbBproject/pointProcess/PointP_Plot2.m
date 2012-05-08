% Plots Point Process
%
%Jeffrey L. Werbin
%Harvard Medical School
%
%Last Update: 9/6/2011
%
%(pos,ImgSize,sf,A,FN,subA)
%

function [img]=PointP_Plot2(pos,ImgSize,sf,A,FN,pixNum,subA)
%PointP_Plot takes a list of corrdinates (x,y) with their corresponding
%uncertainies (dx,dy) and displays them in an image. Each point is represented as
%an ellipse centered at (x,y) with major and minor axes of (dx,dy).
%ImgSize is the size of the original image. sf is a multipicitive scaling
%factor; 10x is a good starting point. A is the gaussian amp, FN is the
%filename to be saved. sub is an optional parameter [xMin,xMax,yMin,yMax]
%that uses a sub section of the points

%Adjust pos values if a sub
if nargin > 6
    pos = pos( pos(:,1)<=subA(2) & pos(:,1)>=subA(1) & pos(:,2)<=subA(4) & pos(:,2)>=subA(3),:);
    pos(:,1)=pos(:,1)-subA(1);
    pos(:,2)=pos(:,2)-subA(3);
    ImgSize=[subA(2)-subA(1),subA(4)-subA(3)];
end

%rescale positions
s =size(pos);
pos = (pos*sf);

%checks to make sure there are no zeros in the sigmas
zlist = find( pos(:,3:4)== 0.0);
m = numel(zlist);
if m > 0
    pos(zlist)= 30/130; %sets equal to 30 nm
end


%pads the positions for easy display
pos(:,1:2)=pos(:,1:2)+pixNum;

%Half the pixNum for correct off setting
pN = floor(pixNum/2);

%Gaussian Amplitude
%A = 10;

%set up the final image assumes that the 2*pixnum padding is already
%account
img = zeros([ImgSize,3]);

for i=1:s(1)
 x = pos(i,1);
 y = pos(i,2);
 xd = mod(x,1);
 yd = mod(y,1);
 x = fix(x);
 y = fix(y);

 %img(x-pN:x+pN,y-pN:y+pN) = img(x-pN:x+pN,y-pN:y+pN)+GaussSpotGen(xd,yd,pos(i,3),pos(i,4),A,pixNum);
 %Places guassians in the green channel
 temp=GaussSpotGen(xd,yd,pos(i,3),pos(i,4),A,pixNum);
 if sum(sum(isnan(temp)))
     temp(:,:)=0.0;
 end
 img(x-pN:x+pN,y-pN:y+pN,2) = img(x-pN:x+pN,y-pN:y+pN,2)+temp;
 
 %Places center indicators in the blue channel
 img(x,y,3)=255;
 
end
 %img(:,:,2) = img(:,:,2)*((65535)/(max(max(img(:,:,2))))); %rescales to 16 bit
 %img(:,:,2) = img(:,:,2)*((255)/(max(max(img(:,:,2))))); %rescales to 8 bit
 %img = img(pixNum:ImgSize(1)+pixNum,pixNum:ImgSize(2)+pixNum,:);
 %imwrite(uint8(img),FN,'TIFF');
 hImage = imshow(img);
    
end

