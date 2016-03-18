function bandsDetect(image)
%bandsDetect uses matched filter followed by watershed segmentation
%to identify lanes and bands on Universal STELA gel with mannual
%adjustment. It then calculates some statistics and give the output with
%each band's position and size.
%
%   [ratio,short20Size,avgBandSize,bandStat]=bandsDetect(image)
%
%   Input: 
%       image: 2D matrix of imput image (RGB can be converted to gray image)
% 
%   Output: 
%       ratio: the ratio of the short telomeres to the rest of the telomeres in a sample
%       short20Size: size of the shortest 20% telomeres in a sample (Kb)
%       avgBandSize: average telomere length in Kb
%       bandStat: a structure with each band's position and size (Kb)
%
%   For runing the program, image matrix is taken from:
%   im=imread('Path');
%
% Ning 10/2015

close all
% Convert input to gray scale image and apply local adaptive thresholding
if size(image,3)==4
    image(:,:,4)=[];
end
if size(image,3)==3

    image = rgb2gray(image);
end
close

imshow(image);
im=imcrop;

% Imcomplement the image if its background is white
% IM always has black background
IM = mat2gray(im);
if mean(mean(IM)) > 0.6
    IM=imcomplement(mat2gray(IM));
end

adjIM=imadjust(IM,[min(min(IM)),max(max(IM))],[0,1]);
[X,Y] = meshgrid(1:size(adjIM,2),1:size(adjIM,1));
[fineX,fineY] = meshgrid(1:.2:size(adjIM,2),1:.2:size(adjIM,1));
fineIM = interp2(X,Y,adjIM,fineX,fineY);

figure,imshow(imcomplement(fineIM),[])

% Apply matched filter f, gaussian kernel defined as sigmax=13, sigmay=6
f=anisoGaussian2Dkernel(0,0,1,13,6,0,-39:39,-18:18);
fIM = imfilter(fineIM,f);
bandMap=zeros(size(fIM));

% Project all pixels intensity to x axis
intensityProfile=zeros(1,size(fIM,2));
for j=1:size(fIM,2)
    for i=1:size(fIM,1)
        intensityProfile(j) = intensityProfile(j) + fIM(i,j);
    end
end
% figure, plot(intensityProfile),title('Intensity Profile of filtered image')

thres1=0.05*(max(intensityProfile)-min(intensityProfile));
laneCenterLoc=find(~watershed(imhmin(intensityProfile,thres1)));
disp(strcat(num2str(length(laneCenterLoc)), ' lane(s) detected'))

for k = 1:length(laneCenterLoc)
    % Take average intersity value of laneCenter +/- 5
    laneCenter = zeros(size(fIM,1),1);
    leftBoundary = max(0,laneCenterLoc(k)-10);
    rightBoundary = min(size(fIM,2),laneCenterLoc(k)+10);
    for rangeLoc = leftBoundary:1:rightBoundary
        laneCenter = laneCenter + fIM(:,rangeLoc);
    end
    %figure,plot(laneCenter)
    %findpeaks(laneCenter)

    % Needs smarter threshold to eliminate the noise without hurting peaks
    thres2=0.05*(max(laneCenter)-min(laneCenter));
    bandLoc=find(~watershed(imhmin(laneCenter,thres2)));
    bandsIntersity = laneCenter(bandLoc);
    avgIntensity = mean(laneCenter(bandLoc));
    for m = 1:length(bandLoc)
        % Mark double bands if the band intersity is higher than 1.5
        % average
        if bandsIntersity(m) > avgIntensity*1.5
            bandMap(bandLoc(m),laneCenterLoc(k)) = 2;
        else
            bandMap(bandLoc(m),laneCenterLoc(k)) = 1;
        end
    end
end

% Show preliminary detected result for manual modification
figure,imshow(imcomplement(fineIM),[]),hold on
% p goes row by row and q goes column by column
for p=1:size(fIM,1)
    for q=1:size(fIM,2)
        if bandMap(p,q)==1
            plot(q,p,'r.'),hold on
        else if bandMap(p,q)==2
                plot(q,p,'g.'),hold on
            end
        end
    end
end

% Manual adjustment of detected bands
bandsAdd = input('Enter the band(s) number you want to add > ');
for v=1:bandsAdd
    [x,y] = ginput(1);
    x = round(x,0);
    y = round(y,0);
    bandMap(y,x)=1;
end

bandsDel = input('Enter the band(s) number you want to delete > ');
for w=1:bandsDel
    [x,y] = ginput(1);
    x = round(x,0);
    y = round(y,0);
    if bandMap(y,x)==1
        bandMap(y,x)=0;
    else
        dist = 100;
        rowUp = max(1,y-dist/2);
        rowDown = min(size(fIM,1),y+dist/2);
        colLeft = max(1,x-dist/2);
        colRight = min(size(fIM,2),x+dist/2);
        adjMat = bandMap(rowUp:rowDown,colLeft:colRight);
        for r = 1:size(adjMat,1)
            for s = 1:size(adjMat,2)
                if adjMat(r,s)==1;
                    % Find the closest spot to the clicking point
                    newDist = (((rowDown-rowUp)/2+1-r)^2+((colRight-colLeft)/2+1-s)^2)^0.5;
                    if newDist < dist
                        dist = newDist;
                        % Need smarter index finding code!
                        bandX = r-((rowDown-rowUp)/2+1);
                        bandY = s-((colRight-colLeft)/2+1);
                    end
                end
            end
        end
        bandMap(y+bandX,x+bandY) = 0;
    end
end



% Plot identified band centers and count the ratio
impixelinfo;
% side = input('Is the marker lane on the left or right? > ','s');
limitX = input('Enter the left/right boundary to eliminate marker lane > ');
limitY = input('Enter the threshold marker size > ');

countShort = 0;
countTotal = 0;
bandStat = struct('bandPos',[],'bandSize',[],'inputIM',im);

markerNum = input('Please enter the number of marker bands (7 or 8) > ');
disp('Please click on position of each marker band')
[markerSize, markerPos] = ginput(markerNum);
switch markerNum
    case 7
        markerSize = [9.4, 6.1, 5.4, 4.4, 3.3, 1.6, 0.8]';
    case 8
        markerSize = [18.8, 9.4, 6.1, 5.4, 4.4, 3.3, 1.6, 0.8]';
end

close all

limitY = markerPos(markerSize==limitY);
% Linear regression y=ax+b for any two adjacent marker pair (x:markerPos,
% y:markerSize)
G(:,1) = markerPos;
G(:,2) = ones(markerNum,1);
para = struct('slope',[],'intercept',[]);
for num = 1:markerNum-1
    regPara = G(num:num+1,:)\markerSize(num:num+1);
    para(num).slope = regPara(1);
    para(num).intercept = regPara(2);
end
% para1 = G(1:2,:)\markerSize(1:2);
% para2 = G(2:3,:)\markerSize(2:3);
% para3 = G(3:4,:)\markerSize(3:4);
% para4 = G(4:5,:)\markerSize(4:5);
% para5 = G(5:6,:)\markerSize(5:6);
% para6 = G(6:7,:)\markerSize(6:7);


% Plot identification results
figure;
imshow(imcomplement(fineIM)),hold on

% Get bands statistics
% p goes row by row and q goes column by column
for p=1:size(fIM,1)
    for q=1:size(fIM,2)
        if bandMap(p,q)==1 && q > limitX
            plot(q,p,'r.'),hold on            
%             switch side
%                 case 'left'
                    countTotal=countTotal+1;
                    bandStat(countTotal).bandPos = p;

                    if p <= markerPos(1)
                        bandStat(countTotal).bandSize = para(1).slope*p+para(1).intercept;
                    else if p > markerPos(markerNum)
                            bandStat(countTotal).bandSize = para(markerNum-1).slope*p+para(markerNum-1).intercept;
                        else
                            for num = 1:markerNum-1
                                if p > markerPos(num) && p <= markerPos(num+1)
                                    bandStat(countTotal).bandSize = para(num).slope*p+para(num).intercept;
                                end
                            end
                        end
                    end
                    
                    if p > limitY
                        countShort=countShort+1;
                    end
                    
%                 case 'right'
%                     if q < limitX
%                         countTotal=countTotal+1;
%                         if p > limitY
%                             countShort=countShort+1;
%                         end
%                     end
%                 otherwise
%                     warning('There should be a marker lane on either left or right side.')
%             end
        
        else if bandMap(p,q)==2 && q > limitX
                plot(q,p,'g.'),hold on 
                for repeat = 1:2
                    countTotal=countTotal+1;
                    bandStat(countTotal).bandPos = p;

                    if p <= markerPos(1)
                        bandStat(countTotal).bandSize = para(1).slope*p+para(1).intercept;
                    else if p > markerPos(markerNum)
                            bandStat(countTotal).bandSize = para(markerNum-1).slope*p+para(markerNum-1).intercept;
                        else
                            for num = 1:markerNum-1
                                if p > markerPos(num) && p <= markerPos(num+1)
                                    bandStat(countTotal).bandSize = para(num).slope*p+para(num).intercept;
                                end
                            end
                        end
                    end

                    if p > limitY
                        countShort=countShort+2;
                    end
                end
                
            end
            
        end
    end
end

for ii = 1:size(bandStat,2)
    if bandStat(ii).bandSize<0
        bandStat(ii).bandSize=0;
    end
end

% short20Pos = bandStat((round((0.8*size(bandStat,2)),0))).bandPos;
short20Size = round(bandStat((round((0.8*size(bandStat,2)),0))).bandSize,2);
avgBandSize = round(mean([bandStat(:).bandSize]),2);

% Plot a threshold line
plot(1:q,limitY*ones(1,q),'b')
hold on
% plot(1:q,short20Pos*ones(1,q),'m')
% hold on
% impixelinfo;
fileName = input('Enter the file name to save statistics > ','s');
title(fileName)
g(1) = gcf;

% Calculate the ratio based on the input threshold
ratio = round(countShort/countTotal*100,2);

% Size distribution histogram
figure, 
h = histogram([bandStat(:).bandSize],'BinWidth',1);
% Normalization
bar(h.BinEdges(2:size(h.BinEdges,2))-0.5,h.Values./sum(h.Values)*100,1)
xlabel('Telomere size (Kb)')
ylabel('Percentage of detected bands (%)')
title(fileName)
g(2) = gcf;

save(fileName,'ratio','avgBandSize','bandStat','im','short20Size')
savefig(g,fileName)
