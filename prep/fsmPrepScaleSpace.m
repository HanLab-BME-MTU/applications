function [frGR,bgGR,yN,xN]=fsmPrepScaleSpace(I,strg,counter,sigmaOne,sigmaTwo,Q) 

% scaleSpace: segments an image to speckled and none-speckled part
% finds all the speckles in the image after applying statistical test
%
% SYNOPSIS   [frGR,bgGR,yN,xN]=scaleSpace(sigmaOne,sigmaTwo)
%
% INPUT      I        : raw data image 
%            strg     : format string for the correct file numbering
%            counter  : image number          
%            sigmaOne : sigma of initial filtering (ex. 1)
%            sigmaTwo : sigma of subsequatial filtering (ex. 1.06)
%            Q        : quantile (ex. 1.96)
%            
% OUTPUT     frGR     : speckled part of the cell   
%            bgGR     : none-speckled part of the cell (background)
%            yN       : Y speckle coordinates 
%            xN       : X speckle coordinates
%
% DEPENDENCES   scaleSpace uses { gauss2d, locmax2d, edge, imfill }
%               scaleSpace is used by { }
%
% REMARKS       11 debug figures
%
% Alexandre Matov, March 26th, 2003

DEBUG=0;
if nargin==0
    DEBUG=1;
    sigmaOne=1; 
    sigmaTwo=1.06;
    Q=1.96;   
    [fileName,dirName] = uigetfile('*.tif','Choose an image');
    I=imread([dirName,filesep,fileName]);
    I=double(I);   
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% image segmentation

If=Gauss2D(I,sigmaOne); % filter the image with sigmaOne
MN=mean(If(:)); % find the Mean Intensity over the whole filtered image
If1=If<=MN; % thresholded image - under the mean

BWs1 = double(edge(If1)); % edge detection based on the mean as a threshold
BWdfill = double(imfill(BWs1,'holes')); % fill in the holes found after the edge detection
BWdfill=~BWdfill; % the reverse of the filled image

% thresholding
MNN=mean(BWdfill(:)); % find the Mean Intensity over image with filled holes
BWdfill=BWdfill>MNN; % thresholded image - above the mean

BG=(If1)&(BWdfill); % logical And between the two thresholded images above 
BG=~BG; % reverse
bg1=double(imfill(double(BG),'holes')); % fill in holes
bg2=~bg1; % reverse
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scale space speckle detection

% filter sigma 2
I2=Gauss2D(I,sigmaTwo);
% substract
Isub=If-I2;
% clipping, for there are some negative values after substraction of the (more) filtered image
Ifsss=minZero(Isub); % local function
% set border (5pixels) to 0
Ifsss(1:6,:)=0;
Ifsss(end-5:end,:)=0;
Ifsss(:,1:6)=0;
Ifsss(:,end-5:end)=0;
% foreground
frGR=Ifsss.*bg1; 
% background
bgGR=Ifsss.*bg2; 
% detect the local maxima of the BG/noise speckles and find their coordinates
Imaxsss=locmax2d(bgGR,[5 5]);
u=find(ne(Imaxsss,0));
% noise speckles (vector of the local maxima intensities)
v=Ifsss(u); 
% detect the local maxima of the FG/real speckles and find their coordinates
ImaxFsss=locmax2d(frGR,[5 5]);
uF=find(ne(ImaxFsss,0));
% real speckles (vector of the local maxima intensities)
vF=Ifsss(uF);
% calculation of delta I critical
meanSig=mean(v); % E(x) of Noise Speckles
stdSig=std(v); % STD of Noise Speckles
sumMandSs=meanSig + Q*stdSig; % Significance Test
% reject the insignificant speckles
Mask=ImaxFsss>sumMandSs;
locMax=Mask.*ImaxFsss;
[yN xN]=find(ne(locMax,0)); % final result: confirmed speckles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize cands structure
cands=struct(...
    'Lmax',[0 0],...                 % Local maximum position - [y x]
    'Bkg1',[0 0],...                 % First local minimum position - [y x]
    'Bkg2',[0 0],...                 % Second local minimum position - [y x]
    'Bkg3',[0 0],...                 % Third local minimum position - [y x]
    'ILmax',0,...                    % Local maximum intensity
    'IBkg',0,...                     % Mean background intensity
    'deltaI',0,...                   % Intensity difference: ILmax-IBkg
    'deltaICrit',0,...               % Critical intensity difference as calculated with the noise model
    'sigmaLmax',0,...                % Error on local maximum intensity
    'sigmaBkg',0,...                 % Error on background intensity 
    'status',0,...                   % Significance of the local maximum: 1, speckle; 0, weak local maximum
    'speckleType',0);                % Describes the level of the speckle hierarchical structure
    
% fill in the cands
for i=1:length(yN)
    cands(i).Lmax=[yN(i) xN(i)];
    cands(i).ILmax=locMax(yN(i),xN(i));
    cands(i).deltaI=locMax(yN(i),xN(i));
    cands(i).status=1;
    cands(i).speckleType=1;
end
 
% Save speckle information (cands and locMax) to disk for future use
if DEBUG==0    
    indxStr=sprintf(strg,counter);
    eval(strcat('save cands',filesep,'cands',indxStr,'.mat cands;')); % Save speckle info
    eval(strcat('save locMax',filesep,'locMax',indxStr,'.mat locMax;')); % Save loc max positions
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% debug figures
if DEBUG==1
    
    close all
    
    figure,imshow(If1)
    title('Thresholded image (threshold - the Mean) after filtering gauss2d(I,1)');
    
    figure,hist(If(:),[min(If(:)):1:max(If(:))]);% yes
    title('Intensity Histogram of the Image')
    
    figure, imshow(BWdfill); 
    title('~binary image with filled holes');
    
    figure,imshow(BG)
    title('background');% yes
    
    figure,imshow(bg1)
    title('foreground is 1');% yes
    
    minI=min(Ifsss(:));
    maxI=max(Ifsss(:));
    
    figure,imshow(frGR,[minI,maxI])% yes
    title('foreground Only (after segmentation)')
     
    [y x]=find(ne(Imaxsss,0));
    
    figure,imshow(bgGR,[minI maxI])% yes
    hold on
    plot(x,y,'r.')
    hold off
    title('noise speckles in the BG')
    
    figure,hist(v);% yes
    title('histogram of Noise (BG) Speckles')
 
    figure,hist(vF); % yes
    title('histogram of ForeGround Speckles')
   
    [yF xF]=find(ne(ImaxFsss,0));
    [yN xN]=find(ne(ImaxN,0));
    
    figure,imshow(frGR,[minI maxI])% yes
    hold on
    plot(xF,yF,'r.')
    plot(xN,yN,'g.')
    hold off
    title('significant speckles GREEN, rejected RED')
   
    figure,imshow(If(5:end-4,5:end-4),[])% yes
    hold on
    plot(xN-4,yN-4,'g.')
    hold off
    title('Significant Speckles (overlaid on the original image)')
end

% local function 
function M=minZero(M)
m=size(M,1);
n=size(M,2);
for i=1:m
    for j=1:n
        if M(i,j)<0
            M(i,j)=0;
        end
    end
end