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
% Ning 07/2016

%% Image loading and preprocessing
close all

% load history
if ~isdeployed
    historyFile = fullfile(fileparts(mfilename('fullpath')),'imagePathHistory.mat');
    if exist(historyFile, 'file')
        history = load(historyFile);
        pathName = history.pathName;
        [fullFileName, pathName] = uigetfile({'*.tif'; '*.jpg'}, 'Select the image for analysis', pathName);
    else
        [fullFileName, pathName] = uigetfile({'*.tif'; '*.jpg'}, 'Select the image for analysis');
    end
    
    dataFilePath = fullfile(pathName, fullFileName);
    [pathName, fileName, ext] = fileparts(dataFilePath);
    
    % save history
    if ~isempty(pathName) && all(pathName ~= 0)
        historyFile = fullfile(fileparts(mfilename('fullpath')), 'imagePathHistory.mat');
        save(historyFile, 'pathName');
    end
end

% Try to store and retrieve path history in a better way
% save pathName pathName

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
title('Input Image');

% Anisogaussian filter is NOT applied!!
% Apply matched filter f, gaussian kernel defined as sigmax=13, sigmay=6
% f=anisoGaussian2Dkernel(0,0,1,13,6,0,-39:39,-18:18);
% fIM = imfilter(fineIM,f);
% figure, imshow(imcomplement(fIM),[]);
% title('Filtered Image');

fIM = fineIM;

%% Lanes and bands detection

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
% Carefully adjust laneThresh and bandThresh
laneThresh = 0.01*(max(intensityProfile)-min(intensityProfile));
laneCenterLoc = find(~watershed(imhmin(intensityProfile,laneThresh)));
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
    % Determine the sensitivity of bands detection
    bandThresh = 0.03*(max(bandProfile)-min(bandProfile));
    bandLoc = find(~watershed(imhmin(bandProfile,bandThresh)));
    
    for m = 1:length(bandLoc)
        bandMap(bandLoc(m),laneCenterLoc(k)) = 1;
    end
    
    clear bandLoc

end

% Show preliminary detected result for manual modification
bandMapPlot(fIM, bandMap);

impixelinfo;
% Marker lane is required to be on the left
disp('Click to define the boundary of marker lane.');
[markerBoundary, limitXX] = ginput(1);
hold on
yAxis = ylim;
plot(markerBoundary * ones(1,2), [yAxis(1), yAxis(2)], 'b-')

% Reset marker lane in case markers are not detected??


%% Manual adjustment of detected bands
for bandsAddCound = 1:100
    bandsAdd = input('Do you want to manually add a band (y/n) > ', 's');
    if lower(bandsAdd) == 'y'
        disp('Please click on the image to add a single band.');
        [x,y] = ginput(1);
        x = round(x,0);
        y = round(y,0);
        bandMap(y,x)=1;
        bandMapPlot(fIM, bandMap);
        hold on
        plot(markerBoundary * ones(1,2), [yAxis(1), yAxis(2)], 'b-')
    else if lower(bandsAdd) == 'n'
            break
        else
            continue
        end
    end
end


for bandsDelCound = 1:100
    bandsDel = input('Do you want to manually delete bands (y/n) > ', 's');
    if lower(bandsDel) == 'y'
        disp('Please choose the region with bands you want to delete and then double click.');
        [delRegion, boundaryInfo] = imcrop;
        for xCord = max(1, floor(boundaryInfo(1))) : min(size(fIM, 2), ceil(boundaryInfo(1)+boundaryInfo(3)))
            for yCord = max(1, floor(boundaryInfo(2))) : min(size(fIM, 1), ceil(boundaryInfo(2)+boundaryInfo(4)))
                if bandMap(yCord, xCord) == 1
                    bandMap(yCord, xCord) = 0;
                end
            end
        end
        bandMapPlot(fIM, bandMap);
        hold on
        plot(markerBoundary * ones(1,2), [yAxis(1), yAxis(2)], 'b-')
    else if lower(bandsAdd) == 'n'
            break
        else
            continue
        end
    end
end

% Simple repeat to make correction
for bandsAddCound = 1:100
    bandsAdd = input('Do you want to manually add a band (y/n) > ', 's');
    if lower(bandsAdd) == 'y'
        disp('Please click on the image to add a single band.');
        [x,y] = ginput(1);
        x = round(x,0);
        y = round(y,0);
        bandMap(y,x)=1;
        bandMapPlot(fIM, bandMap);
        hold on
        plot(markerBoundary * ones(1,2), [yAxis(1), yAxis(2)], 'b-')
    else if lower(bandsAdd) == 'n'
            break
        else
            continue
        end
    end
end


for bandsDelCound = 1:100
    bandsDel = input('Do you want to manually delete bands (y/n) > ', 's');
    if lower(bandsDel) == 'y'
        disp('Please choose the region with bands you want to delete and then double click.');
        [delRegion, boundaryInfo] = imcrop;
        for xCord = max(1, floor(boundaryInfo(1))) : min(size(fIM, 2), ceil(boundaryInfo(1)+boundaryInfo(3)))
            for yCord = max(1, floor(boundaryInfo(2))) : min(size(fIM, 1), ceil(boundaryInfo(2)+boundaryInfo(4)))
                if bandMap(yCord, xCord) == 1
                    bandMap(yCord, xCord) = 0;
                end
            end
        end
        bandMapPlot(fIM, bandMap);
        hold on
        plot(markerBoundary * ones(1,2), [yAxis(1), yAxis(2)], 'b-')
    else if lower(bandsAdd) == 'n'
            break
        else
            continue
        end
    end
end

%% Band size annotation

% Enter threshold to calculate shortest telomere ratio
% shortThresh = input('Enter the threshold marker size > ');
shortThresh = 1.6;
markerNum = 0;
for markerY = 1:size(fIM,1)
    for markerX = 1:round(markerBoundary)
        if bandMap(markerY, markerX) == 1
            markerNum = markerNum + 1;
            markerPos(markerNum) = markerY;
        end
    end
end
markerPos = markerPos';

switch markerNum
    case 7
        markerSize = [9.4, 6.1, 5.4, 4.4, 3.3, 1.6, 0.8]';
    case 8
        markerSize = [18.8, 9.4, 6.1, 5.4, 4.4, 3.3, 1.6, 0.8]';
    otherwise
        error('Unusal number of marker bands detected.')
end

limitY = markerPos(markerSize == shortThresh);
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
        if bandMap(p,q)==1 && q > markerBoundary
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

for negativeBand = 1:size(bandStat,2)
    if bandStat(negativeBand).bandSize<0
        bandStat(negativeBand).bandSize=0;
    end
end

clear bandMap

%% Overlapping bands detection and counting
% Interval setup by percentile of bandsize
minBandSize = min([bandStat(:).bandSize]);
maxBandSize = max([bandStat(:).bandSize]) + 0.0001;
step = 5;
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
        overlappingBand(addCount).count = round(overlappingBand(addCount).bandIntensity/medIntensity);
        bandStat([bandStat(:).index] == overlappingBand(addCount).index).count = overlappingBand(addCount).count;
        %         bandIndexOld = bandIndex;
        %         bandIndex = bandIndex + overlappingBand(addCount).count - 1;
        %         for indexAdd = 1:(bandIndex - bandIndexOld)
        %             bandStat(indexAdd) = bandStat([bandStat(:).index] == overlappingBand(addCount).index);
        %             bandStat(indexAdd).count = 'Overlapping Band';
        %         end
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
    else if bandStat(i).count == 2
            plot(bandStat(i).bandPos(1), bandStat(i).bandPos(2), 'g.')
        else if bandStat(i).count > 2
                plot(bandStat(i).bandPos(1), bandStat(i).bandPos(2), 'm.')
            end
        end
    end
    
    countTotal = countTotal + bandStat(i).count;
    if bandStat(i).bandPos(2) > limitY
        countShort = countShort + bandStat(i).count;
    end
end
% Plot a threshold line
plot(1:q,limitY*ones(1,q),'b')
fprintf('Input file name: %s%s\n', fileName, ext)
avgBandSize = [bandStat(:).bandSize]*[bandStat(:).count]'/sum([bandStat(:).count]);
fprintf('The average telomere length is %.2f kb.\n', avgBandSize)

% Calculate shortest telomeres (below 1.6kb) ratio
ratio = countShort/countTotal*100;
fprintf('The ratio of shortest telomere below %.1fkb is %.2f%%.\n',shortThresh, ratio)

% Calculate the shortest 20% band size threshold
[sortBandSize, order] = sort([bandStat(:).bandSize]);
sortCount = [bandStat(order).count];
short20Point = sum(sortCount) * 0.2;
cumulativeSumCount = cumsum(sortCount);
boundaryPoint = find(cumulativeSumCount <= short20Point, 1, 'last');
if cumulativeSumCount(boundaryPoint) == short20Point
    short20Size = mean(sortBandSize([boundaryPoint, boundaryPoint+1]));
else
    short20Size = sortBandSize(boundaryPoint+1);
end

fprintf('The shortest 20%% telomere threshold is %.2f kb.\n', short20Size)

save([fullfile(pathName, fileName), '.mat'],'ratio', 'avgBandSize', 'bandStat', 'imInput', 'short20Size', 'dataFilePath')

% plot(1:q,short20Pos*ones(1,q),'m')
% hold on
% impixelinfo;
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

savefig(g,fullfile(pathName, fileName))


function bandMapPlot(image, bandMap)
close all
figure,
set(gcf,'position',get(0,'screensize'))
subplot(1,2,1)
imshow(imcomplement(image),[])
title('Input image')
subplot(1,2,2)
imshow(imcomplement(image),[]),
title('Detected bands'), hold on
% p goes row by row and q goes column by column
for p = 1:size(image,1)
    for q = 1:size(image,2)
        if bandMap(p,q) == 1
            plot(q, p, 'r.'),hold on
        end
    end
end
hold off