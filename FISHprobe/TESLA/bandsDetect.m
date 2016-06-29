function bandsDetect(varargin)
%BANDSDETECT uses watershed segmentation to identify lanes and bands on 
% Universal STELA gel with mannual adjustment. It can calculate and store 
% statistics of the bands in an efficient way.
%
%   bandsDetect(varargin)
%
%   Input(optional): pathName
% 
%   Output: 
%       ratio: the ratio of the short telomeres to the rest of the telomeres in a sample
%       short20Size: size of the shortest 20% telomeres in a sample (Kb)
%       avgBandSize: average telomere length in Kb
%       bandStat: a structure with each band's position and size (Kb)
%       imInput: input image (cropped) used for analysis
%       Figure: Bands detection results and telomere size distribution
%
% Ning 04/2016

%% Image reading and preprocessing
close all

if nargin == 0
    [fileName, pathName] = uigetfile('*.tif', 'Select the STELA image');
else
    pathName = varargin{1};
    [fileName, pathName] = uigetfile('*.tif', 'Select the STELA image', pathName);
end

% Try to store and retrieve path history in a better way
% save pathName pathName

dataFilePath = fullfile(pathName, fileName);
image = imread(dataFilePath);

% Convert input to single frame 2D gray scale image
if size(image,3) == 4
    image(:,:,4) = [];
else if size(image,3) == 3
        image = rgb2gray(image);
    end
end

figure, imshow(image,[]);
imInput=imcrop;

% Imcomplement the image if its background is white
% IM, fineIM, fIM always has black background
IM = mat2gray(imInput);
if mean(IM(:)) > 0.5
    IM=imcomplement(mat2gray(IM));
end

adjIM = imadjust(IM,[min(IM(:)),max(IM(:))], [0,1]);
[X,Y] = meshgrid(1:size(adjIM,2), 1:size(adjIM,1));
[fineX,fineY] = meshgrid(1:.2:size(adjIM,2), 1:.2:size(adjIM,1));
fineIM = interp2(X,Y,adjIM,fineX,fineY);

figure, imshow(imcomplement(fineIM),[]);
title('Imput Image');

% Anisogaussian filter is NOT applied!! 
% Apply matched filter f, gaussian kernel defined as sigmax=13, sigmay=6
% f=anisoGaussian2Dkernel(0,0,1,13,6,0,-39:39,-18:18);
% fIM = imfilter(fineIM,f);
% figure, imshow(imcomplement(fIM),[]);
% title('Filtered Image');

fIM = fineIM;

%% Bands detection

bandMap = zeros(size(fIM));
% Project all pixels intensity to x axis
intensityProfile = zeros(1,size(fIM,2));
for j = 1:size(fIM,2)
    for i = 1:size(fIM,1)
        intensityProfile(j) = intensityProfile(j) + fIM(i,j);
    end
end
% figure, plot(intensityProfile),title('Intensity Profile of filtered image')

% Should all range relavant calculations be normalized???
% Instead of pointing out the center, can you watershed out a range that
% indicates bands???

% Lanes and bands detection.
% Carefully adjust thresh1 and thresh2
thres1 = 0.05*(max(intensityProfile)-min(intensityProfile));
laneCenterLoc = find(~watershed(imhmin(intensityProfile,thres1)));
disp(strcat(num2str(length(laneCenterLoc)), ' lane(s) detected'))

halfBandRange = max(1, round(0.3*(median(diff(laneCenterLoc))/2)));
for k = 1:length(laneCenterLoc)
    % Take average intersity value of laneCenter +/- halfBandRange
    bandProfile = zeros(size(fIM,1),1);
    laneLeftBoundary = max(0,laneCenterLoc(k)-halfBandRange);
    laneRightBoundary = min(size(fIM,2),laneCenterLoc(k)+halfBandRange);
    count = 0;
    for bandRange = laneLeftBoundary:laneRightBoundary
        bandProfile = bandProfile + fIM(:,bandRange);
        count = count + 1;
    end
    bandProfile = bandProfile/count;
    %figure,plot(laneCenter)
    %findpeaks(laneCenter)

    % Needs smarter threshold to eliminate the noise without hurting peaks
    thres2 = 0.05*(max(bandProfile)-min(bandProfile));
    bandLoc = find(~watershed(imhmin(bandProfile,thres2)));
    
    for m = 1:length(bandLoc)
        bandMap(bandLoc(m),laneCenterLoc(k)) = 1;
    end

% Count overlapping bands vertically
%     bandsIntersity = laneCenter(bandLoc);
%     avgIntensity = mean(laneCenter(bandLoc));
%     for m = 1:length(bandLoc)
%         % Mark double bands if the band intersity is higher than 1.5
%         % average
%         if bandsIntersity(m) > avgIntensity*1.5
%             bandMap(bandLoc(m),laneCenterLoc(k)) = 2;
%         else
%             bandMap(bandLoc(m),laneCenterLoc(k)) = 1;
%         end
%     end
    
end


% Show preliminary detected result for manual modification
figure,imshow(imcomplement(fineIM),[]),hold on
% p goes row by row and q goes column by column
for p = 1:size(fIM,1)
    for q = 1:size(fIM,2)
        if bandMap(p,q) == 1
            plot(q,p,'r.'),hold on
%         else if bandMap(p,q)==2
%                 plot(q,p,'g.'),hold on
%             end
        end
    end
end

%% Manual adjustment of detected bands
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
    if bandMap(y,x) == 1
        bandMap(y,x) = 0;
    else
        dist = 100;
        rowUp = max(1,y-dist/2);
        rowDown = min(size(fIM,1),y+dist/2);
        colLeft = max(1,x-dist/2);
        colRight = min(size(fIM,2),x+dist/2);
        adjMat = bandMap(rowUp:rowDown,colLeft:colRight);
        for r = 1:size(adjMat,1)
            for s = 1:size(adjMat,2)
                if adjMat(r,s) == 1;
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

%% Band size annotation

impixelinfo;
% Marker lane is required to be on the left
limitX = input('Enter the left boundary to eliminate marker lane > ');

close all
% Plot identification results
figure;
imshow(imcomplement(fineIM)),hold on
for p=1:size(fIM,1)
    for q=1:size(fIM,2)
        if bandMap(p,q)==1 && q > limitX
            plot(q,p,'r.'),hold on
        end
    end
end

% Enter threshold to calculate shortest telomere ratio
% shortThresh = input('Enter the threshold marker size > ');
shortThresh = 1.6;
markerNum = input('Please enter the number of marker bands (7 or 8) > ');
disp('Please click on position of each marker band')
[markerSize, markerPos] = ginput(markerNum);
switch markerNum
    case 7
        markerSize = [9.4, 6.1, 5.4, 4.4, 3.3, 1.6, 0.8]';
    case 8
        markerSize = [18.8, 9.4, 6.1, 5.4, 4.4, 3.3, 1.6, 0.8]';
end

limitY = markerPos(markerSize==shortThresh);
% Linear regression y=ax+b for any two adjacent marker pair (x:markerPos, y:markerSize)
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


% Band size annotation
bandStat = struct('index', [], 'bandPos', [], 'bandSize', [], 'bandIntensity', [], 'count', []);
bandIndex = 0;
for q=1:size(fIM,2)
    for p=1:size(fIM,1)
        if bandMap(p,q)==1 && q > limitX
            bandIndex = bandIndex + 1;
            bandStat(bandIndex).index = bandIndex;
            % Bands positions are recorded in cartesian coordinate
            bandStat(bandIndex).bandPos = [q, p]; 
            bandStat(bandIndex).count = 1;
            
            if p <= markerPos(1)
                bandStat(bandIndex).bandSize = para(1).slope*p+para(1).intercept;
            else if p > markerPos(markerNum)
                    bandStat(bandIndex).bandSize = para(markerNum-1).slope*p+para(markerNum-1).intercept;
                else
                    for num = 1:markerNum-1
                        if p > markerPos(num) && p <= markerPos(num+1)
                            bandStat(bandIndex).bandSize = para(num).slope*p+para(num).intercept;
                        end
                    end
                end
            end
            
            % Take average intensity around a band seed
            halfBandHeight = round(0.2*halfBandRange);
            bandIntensity = 0;
            bandLeftBoundary = max(0, q - halfBandRange);
            bandRightBoundary = min(size(fIM,2), q + halfBandRange);
            bandTopBoundary = max(0, p - halfBandHeight);
            bandBotBoundary = min (size(fIM,1), p + halfBandHeight);
            count = 0;
            for horiRange = bandLeftBoundary:bandRightBoundary
                for vertiRange = bandTopBoundary:bandBotBoundary
                    bandIntensity = bandIntensity + fIM(vertiRange, horiRange);
                    count = count + 1;
                end
            end
            
            bandStat(bandIndex).bandIntensity = bandIntensity/count;
        end
    end
end

for trivialBand = 1:size(bandStat,2)
    if bandStat(trivialBand).bandSize<0
        bandStat(trivialBand).bandSize=0;
    end
end


%% Overlapping bands detection and counting
% Interval setup by percentile of bandsize
minBandSize = min([bandStat(:).bandSize]);
maxBandSize = max([bandStat(:).bandSize]) + 0.01;
step = 10;
interDistance = (maxBandSize - minBandSize) / step; 
bottom = minBandSize;
for i = 1:step
    top = bottom + interDistance;
    subBandStat = bandStat([bandStat(:).bandSize] >= bottom & [bandStat(:).bandSize] < top);
    % Choice between median intensity and average intensity
    medIntensity = median([subBandStat(:).bandIntensity]);
    
    % Set threshold to select highly intensed bands within the top/bottom interval
    intensityThreshold = medIntensity*1.5;
    overlappingBand = subBandStat([subBandStat(:).bandIntensity] > intensityThreshold);    
    for addCount = 1:numel(overlappingBand)
        overlappingBand(addCount).count = 2;
        bandStat([bandStat(:).index] == overlappingBand(addCount).index).count = 2;
        bandIndex = bandIndex + 1;
        bandStat(bandIndex) = bandStat([bandStat(:).index] == overlappingBand(addCount).index);
        bandStat(bandIndex).count = 'Overlapping Band';
    end
    
    bottom = top;
end

%% Calculate bands statistics

% Plot identification results
close all
figure;
imshow(imcomplement(fineIM)), hold on

countShort = 0;
countTotal = 0;

for i = 1:numel(bandStat)
    if bandStat(i).count == 1
        plot(bandStat(i).bandPos(1), bandStat(i).bandPos(2), 'r.')
        countTotal = countTotal+1;        
        if bandStat(i).bandPos(2) > limitY
            countShort = countShort+1;
        end        
    else if bandStat(i).count == 2
            plot(bandStat(i).bandPos(1), bandStat(i).bandPos(2), 'g.')
            countTotal = countTotal+2;
            if bandStat(i).bandPos(2) > limitY
                countShort = countShort+2;
            end 
        end
    end
end

avgBandSize = round(mean([bandStat(:).bandSize]),2);
fprintf('The average telomere length is %.2f.\n', avgBandSize)

% Calculate the ratio based on the input threshold
ratio = round(countShort/countTotal*100,2);
fprintf('The ratio of shortest telomere below %f is %.2f.\n',shortThresh, ratio)

bandSorted = sort([bandStat(:).bandSize]);
short20Size = round(bandSorted((round((0.2*size(bandSorted,2)),0))),2);
fprintf('The shortest 20%% telomere threshold is %.2f.\n', short20Size)


% Plot a threshold line
plot(1:q,limitY*ones(1,q),'b')
hold on
% plot(1:q,short20Pos*ones(1,q),'m')
% hold on
% impixelinfo;
fileName = input('Enter the file name to save statistics > ','s');
title(fileName)
g(1) = gcf;

% Size distribution histogram
figure, 
h = histogram([bandStat(:).bandSize],'BinWidth',1);
% Normalization
bar(h.BinEdges(2:size(h.BinEdges,2))-0.5,h.Values./sum(h.Values)*100,1)
xlabel('Telomere size (Kb)')
ylabel('Percentage of detected bands (%)')
title(fileName)
g(2) = gcf;

save(fileName,'ratio','avgBandSize','bandStat','imInput','short20Size')
savefig(g,fileName)
