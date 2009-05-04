function [movieInfo]=eb1SpotDetector(runInfo,timeRange,bitDepth,savePlots)
% EB1SPOTDETECTOR locates EB1/EB3 comets in a movie stack
%
%SYNOPSIS [movieInfo]=eb1SpotDetector(runInfo,timeRange,bitDepth,savePlots)
%
%INPUT  runInfo           : structure containing fields .anDir, which gives
%                           the full path to the roi_x directory
%                           and .imDir, which gives the full path to the
%                           folder containing the images for overlay.
%                           if given as [], program will query user for
%                           roi_x directory.
%       timeRange         : row vector of the form [startFrame endFrame]
%                           indicating time range to plot. if not given or 
%                           given as [], tracks from the whole movie will
%                           be displayed. 
%       bitDepth          : bit depth of the images - should be 12, 14, or 16
%       savePlots         : 1 to save overlay plots of detection results, 
%                           0 if not
%
%OUTPUT movieInfo         : nFrames-structure containing x/y coordinates
%       stdList           : nFrames-vector containing the standard
%                           deviation of the difference of Gauss-filtered
%                           images corresponding to each frame. this is
%                           based on either the user-selected ROI (if
%                           provided) or a region estimated to be within the
%                           cell from the background point (if ROI wasn't
%                           provided). both the ROI and bg point are saved
%                           during setupRoiDirectories.m

warningState = warning;
warning('off','MATLAB:divideByZero')

% CHECK INPUT AND SET UP DIRECTORIES

% get runInfo in correct format
if nargin<1 || isempty(runInfo)
    % if not given as input, ask user for ROI directory
    % assume images directory is at same level
    runInfo.anDir=uigetdir(pwd,'Please select ROI directory');
    homeDir=pwd;
    cd(runInfo.anDir);
    cd('..');
    runInfo.imDir=[pwd filesep 'images'];
    cd(homeDir)
else
    % adjust for OS
    if ~isfield(runInfo,'imDir') || ~isfield(runInfo,'anDir')
        error('--eb1SpotDetector: first argument should be a structure with fields imDir and anDir');
    else
        [runInfo.anDir] = formatPath(runInfo.anDir);
        [runInfo.imDir] = formatPath(runInfo.imDir);
    end
end

% count number of images in image directory
[listOfImages] = searchFiles('.tif',[],runInfo.imDir,0);
nImTot = size(listOfImages,1);

% check timeRange input, assign start and end frame
if nargin<2 || isempty(timeRange)
    startFrame = 1;
    endFrame = nImTot;
elseif isequal(unique(size(timeRange)),[1 2])
    if timeRange(1)<=timeRange(2) && timeRange(1)>0 && timeRange(2)<=nImTot
        startFrame = timeRange(1);
        endFrame = timeRange(2);
    else
        startFrame = 1;
        endFrame = nImTot;
    end

else
    error('--eb1SpotDetector: timeRange should be [startFrame endFrame] or [] for all frames')
end
nFrames = endFrame-startFrame+1;

% get image dimensions, max intensity from first image
fileNameIm = [char(listOfImages(1,2)) filesep char(listOfImages(1,1))];
img = double(imread(fileNameIm));
[imL,imW] = size(img);
maxIntensity = max(img(:));

% get bit depth if not given
if nargin < 3 || isempty(bitDepth)
    imgData = imfinfo(fileNameIm);
    bitDepth = imgData.BitDepth;
    disp(['bitDepth estimated to be' bitDepth])
end

% check bit depth to make sure it is 12, 14, or 16 and that its dynamic
% range is not greater than the provided bitDepth
if sum(bitDepth==[12 14 16])~=1 || maxIntensity > 2^bitDepth-1
    error('--eb1SpotDetector: bit depth should be 12, 14, or 16');
end

% check input for savePlots
if nargin<4 || isempty(savePlots)
    savePlots = 1;
end

% make feat directory if it doesn't exist from batch
featDir = [runInfo.anDir filesep 'feat'];
if isdir(featDir)
    rmdir(featDir,'s')
end
mkdir(featDir)
mkdir([featDir filesep 'filterDiff']);    
if savePlots==1
    mkdir([featDir filesep 'overlayImages']);
end


% look for region of interest info from project setup step
if ~exist([runInfo.anDir filesep 'roiMask.tif'],'file')
    % not roi selected; use the whole image
    roiMask = ones(imL,imW);
    roiYX=[1 1; imL 1; imL imW; 1 imW; 1 1];
else
    % get roi edge pixels and make region outside mask NaN
    roiMask = double(imread([runInfo.anDir filesep 'roiMask.tif']));
    roiYX=load([runInfo.anDir filesep 'roiYX']);
    roiYX=roiYX.roiYX;
end


% string for number of files
s1 = length(num2str(length(startFrame:endFrame)));
strg1 = sprintf('%%.%dd',s1);


% START DETECTION

% initialize structure to store info for tracking
movieInfo(nFrames,1).xCoord = 0;
movieInfo(nFrames,1).yCoord = 0;
movieInfo(nFrames,1).amp = 0;
movieInfo(nFrames,1).int = 0;

% get difference of Gaussians image for each frame and standard deviation
% of the cell background, stored in stdList
stdList=zeros(endFrame-startFrame+1,1);
for iFrame = startFrame:endFrame

    % load image and normalize to 0-1
    fileNameIm = [char(listOfImages(iFrame,2)) filesep char(listOfImages(iFrame,1))];
    img = double(imread(fileNameIm))./((2^bitDepth)-1);

    % create kernels for gauss filtering
    blurKernelLow  = fspecial('gaussian', 21, 1);
    blurKernelHigh = fspecial('gaussian', 21, 4);

    % use subfunction that calls imfilter to take care of edge effects
    lowPass = filterRegion(img,roiMask,blurKernelLow);
    highPass = filterRegion(img,roiMask,blurKernelHigh);

    % get difference of gaussians image
    filterDiff = lowPass-highPass;

    % if bg point was chosen and saved, get bgMask from first frame
    if iFrame==startFrame && exist([runInfo.anDir filesep 'bgPtYX.mat'])~=0
        bgPtYX=load([runInfo.anDir filesep 'bgPtYX.mat']);
        bgPtYX=bgPtYX.bgPtYX;
        [bgMask]=eb3BgMask(filterDiff,bgPtYX);
        saveas(gcf,[featDir filesep 'filterDiff' filesep 'bgMask.tif']);
    end
    % if bg point wasn't chosen, use ROI
    if iFrame==startFrame && exist([runInfo.anDir filesep 'bgPtYX.mat'])==0
        bgMask=logical(roiMask);
    end

    stdList(iFrame)=std(filterDiff(bgMask));
    
    indxStr1 = sprintf(strg1,iFrame);
    save([featDir filesep 'filterDiff' filesep 'filterDiff' indxStr1],'filterDiff')
    save([featDir filesep 'stdList'],'stdList')
end


sF=max([startFrame:endFrame]-1,ones(1,nFrames));
eF=min([startFrame:endFrame]+1,nFrames*ones(1,nFrames));


% loop thru frames and detect
for iFrame = startFrame:endFrame

    if iFrame==startFrame
        tic
    end

    indxStr1 = sprintf(strg1,iFrame);
    filterDiff=load([featDir filesep 'filterDiff' filesep 'filterDiff' indxStr1]);
    filterDiff=filterDiff.filterDiff;

    % cut postive component of histogram; make below that nan
    %filterDiff(filterDiff<0) = nan;
    %[cutoffInd, cutOffValueInitInt] = cutFirstHistMode(filterDiff(:),0);
    %noiseFloor = 1.0 * cutOffValueInitInt;
    %filterDiff(filterDiff<noiseFloor) = nan;

    % thickness of intensity slices is average std from filterDiffs over
    % from one frame before to one frame after
    stepSize=mean(stdList(sF(iFrame):eF(iFrame)));
    thresh=3*stepSize;
    
    % we assume each step size down the intensity profile should be on
    % the order of the size of the background std; here we find how many
    % steps we need and what their spacing should be. we also assume peaks
    % should be taller than 3*std
    nSteps = round((nanmax(filterDiff(:))-thresh)/(stepSize));
    threshList = linspace(nanmax(filterDiff(:)),thresh,nSteps);

    % compare features in z-slices startest from the highest one
    for p = 1:length(threshList)-1

        % slice1 is top slice; slice2 is next slice down
        % here we generate BW masks of slices
        if p==1
            slice1 = filterDiff>threshList(p);
        else
            slice1 = slice2;
        end
        slice2 = filterDiff>threshList(p+1);

        % now we label them
        featMap1 = bwlabel(slice1);
        featMap2 = bwlabel(slice2);
        featProp2 = regionprops(featMap2,'PixelIdxList');

        % loop thru slice2 features and replace them if there are 2 or
        % more features from slice1 that contribute
        for iFeat = 1:max(featMap2(:))
            pixIdx = featProp2(iFeat,1).PixelIdxList; % pixel indices from slice2
            featIdx = unique(featMap1(pixIdx)); % feature indices from slice1 using same pixels
            featIdx(featIdx==0) = []; % 0's shouldn't count since not feature
            if length(featIdx)>1 % if two or more features contribute...
                slice2(pixIdx) = slice1(pixIdx); % replace slice2 pixels with slice1 values
            end
        end

    end

    % label slice2 again and get region properties
    featMap2 = bwlabel(slice2);
    featProp2 = regionprops(featMap2,filterDiff,'PixelIdxList','Area');
    %[cutoffInd, cutoffValueArea] = cutFirstHistMode(vertcat(featProp2(:,1).Area),1);

    %imshow(featMap2)

    %figure(2); hist(vertcat(featProp2(:,1).Area),50);

    %         % fill label matrix with area value
    %         areaLabel = zeros(imL,imW);
    %         for iFeat = 1:max(featMap2(:))
    %             pixIdx = featProp2(iFeat,1).PixelIdxList;
    %             areaLabel(pixIdx) = featProp2(iFeat,1).Area;
    %         end
    %         figure(1); hist(vertcat(featProp2(:,1).Area),100)
    %         figure(2); hist(vertcat(featProp2(:,1).MaxIntensity),100)
    %         figure(3); imshow(featMap2>0,[])
    %         figure(4); imshow(areaLabel>0,[])
    %         figure(5); imshow(areaLabel>20 & areaLabel<80,[])


    % here we sort through features and retain only the "good" ones
    % we assume the good features have area > 2 pixels
    goodFeatIdx = find(vertcat(featProp2(:,1).Area)>2);
%    goodFeatIdxI = find(vertcat(featProp2(:,1).MaxIntensity)>2*cutOffValueInitInt);
%    goodFeatIdx = intersect(goodFeatIdxA,goodFeatIdxI);

    % make new label matrix and get props
    featureMap = zeros(imL,imW);
    featureMap(vertcat(featProp2(goodFeatIdx,1).PixelIdxList)) = 1;
    [featMapFinal,nFeats] = bwlabel(featureMap);
    featPropFinal = regionprops(featMapFinal,filterDiff,'PixelIdxList','Area','WeightedCentroid','MaxIntensity'); %'Extrema'

    %             % loop thru features and calc total intensity and peak pixel
    %             for iFeat = 1:nFeats
    %                 % intensity values for pixels in feature iFeat
    %                 featIntens = filterDiff(featPropFinal(iFeat,1).PixelIdxList);
    %
    %                 % integrate intensity over feature
    %                 featPropFinal(iFeat).IntSum = sum(featIntens);
    %
    %                 % make centroid to be intensity maximum of feature
    %                 maxIntPixIdx = featPropFinal(iFeat,1).PixelIdxList(featIntens==max(featIntens));
    %                 [r c] = ind2sub([imL,imW],maxIntPixIdx);
    %                 featPropFinal(iFeat).Centroid = [r c];
    %
    %             end

    % centroid coordinates with 0.5 uncertainties for Khuloud's tracker
    yCoord = 0.5*ones(nFeats,2); xCoord = 0.5*ones(nFeats,2);
    temp = vertcat(featPropFinal.WeightedCentroid);
    yCoord(:,1) = temp(:,2);
    xCoord(:,1) = temp(:,1);

    % area
    featArea = vertcat(featPropFinal(:,1).Area);
    amp = zeros(nFeats,2);
    amp(:,1) = featArea;

    % intensity
    featInt = vertcat(featPropFinal(:,1).MaxIntensity);
    featI = zeros(nFeats,2);
    featI(:,1) = featInt;

    % make structure compatible with Khuloud's tracker
    movieInfo(iFrame,1).xCoord = xCoord;
    movieInfo(iFrame,1).yCoord = yCoord;
    movieInfo(iFrame,1).amp = amp;
    movieInfo(iFrame,1).int = featI;


    indxStr1 = sprintf(strg1,iFrame); % frame

    %save([featDir filesep 'featureInfo' filesep 'featPropFinal' indxStr1],'featPropFinal');

    if savePlots==1
        % use extrema to draw polygon around feats - here we get
        % coordinates for polygon
        %             outline = [featPropFinal.Extrema]; outline = [outline; outline(1,:)];
        %             outlineX = outline(:,1:2:size(outline,2));
        %             outlineY = outline(:,2:2:size(outline,2));



        %plot feat outlines and centroid on image
        fileNameIm = [char(listOfImages(iFrame,2)) filesep char(listOfImages(iFrame,1))];
        img = double(imread(fileNameIm))./((2^bitDepth)-1);

        figure(1);
        imagesc(img);
        hold on
        scatter(xCoord(:,1),yCoord(:,1),'c.'); % plot centroid in cyan
        colormap gray
        plot(roiYX(2),roiYX(1),'w')
        %plot(outlineX,outlineY,'r'); % plot feat outlines in red

        saveas(gcf,[featDir filesep 'overlayImages' filesep 'overlay' indxStr1 '.tif']);
        %saveas(gcf,[featDir filesep 'overlay' filesep 'overlay_g' indxStr1 '_f' indxStr3])

        %             % save image of colored feats
        %             RGB = label2rgb(featMapFinal, 'jet', 'k','shuffle');
        %             figure(2);
        %             imshow(RGB,'Border','tight');
        %
        %             save([featDir filesep 'featMaps' filesep 'feats_g' indxStr1 '_f' indxStr3], 'featMapFinal');
        %             saveas(gcf,[featDir filesep 'featMaps' filesep 'feats_g' indxStr1 '_f' indxStr3 '.tif']);

    end

end
save([featDir filesep 'movieInfo'],'movieInfo');


close(gcf)
warning(warningState);



function filteredIm = filterRegion(im, mask, kernel)

im(mask~=1) = 0;
filteredIm = imfilter(im, kernel);
W = imfilter(double(mask), kernel);
filteredIm = filteredIm ./ W;
filteredIm(~mask) = nan;


function [bgMask]=eb3BgMask(filterDiff,bgPtYX)

% local max detection
fImg=locmax2d(filterDiff,[20 20],1);

% get indices of local maxima
idx=find(fImg);
[r c]=find(fImg);

% calculate percentiles of max intensities to use for rough idea of cell
% region
p1=prctile(fImg(idx),80);
p2=prctile(fImg(idx),90);

% get indices of those maxima within the percentile range ("good" features)
goodIdx=find(fImg(idx)>p1 & fImg(idx)<p2);

% get indices for nearest fifty points to user-selected point
D=createDistanceMatrix([bgPtYX(1) bgPtYX(2)],[r(goodIdx) c(goodIdx)]);
[sD,closeIdx]=sort(D);
closeIdx=closeIdx(1:min(50,length(closeIdx)));

% get convex hull and create ROI from that
K = convhull(c(goodIdx(closeIdx)),r(goodIdx(closeIdx)));
[bgMask,xi,yi]=roipoly(fImg,c(goodIdx(closeIdx(K))),r(goodIdx(closeIdx(K))));


figure(1); imagesc(filterDiff); colormap gray;
hold on
scatter(bgPtYX(2),bgPtYX(1),'*y') % user-selected point
scatter(c(goodIdx),r(goodIdx),'.g') % all "good" features in green
scatter(c(goodIdx(closeIdx)),r(goodIdx(closeIdx)),'r') % nearest fifty to point in red
plot(xi,yi) % plot mask outline in blue


