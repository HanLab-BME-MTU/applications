function [movieInfo]=eb1SpotDetector(runInfo,frameRange,bitDepth,savePlots,overwriteData)

% runInfo: structure containing image and analysis directories
% frameRange: [startFrame endFrame] or [] for whole movie

% sigma for lowpass gauss filtering (if want to try multiple sigmas, use a
% vector such as [4 8 10], for example.  4 seems to work best in most
% cases.
gaussSigma = [4];

warningState = warning;
warning('off','MATLAB:divideByZero')

% CHECK INPUT AND SET UP DIRECTORIES

% check runInfo format and directory names between file systems
if ~isfield(runInfo,'imDir') || ~isfield(runInfo,'anDir')
    error('eb1SpotDetector: runInfo should contain fields imDir and anDir');
else
    [runInfo.anDir] = formatPath(runInfo.anDir);
    [runInfo.imDir] = formatPath(runInfo.imDir);
end

% count number of images in image directory
[listOfImages] = searchFiles('.tif',[],runInfo.imDir,0);
nImTot = size(listOfImages,1);

% check frameRange input, assign start and end frame
if nargin<2 || isempty(frameRange)
    startFrame = 1;
    endFrame = nImTot;
elseif isequal(unique(size(frameRange)),[1 2])
    if frameRange(1)<=frameRange(2) && frameRange(1)>0 && frameRange(2)<=nImTot
        startFrame = frameRange(1);
        endFrame = frameRange(2);
    else
       startFrame = 1;
    endFrame = nImTot;
    end

else
    error('eb1SpotDetector: frameRange should be [startFrame endFrame] or [] for all frames')
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
end

% check bit depth to make sure it is 12, 14, or 16 and that its dynamic
% range is not greater than the provided bitDepth
if sum(bitDepth==[12 14 16])~=1 || maxIntensity > 2^bitDepth-1
    error('eb1SpotDetector: bit depth should be 12, 14, or 16');
end

% check input for savePlots
if nargin<4 || isempty(savePlots)
    savePlots = 1;
end

% check input for overwriteData
if nargin<5 || isempty(overwriteData)
    overwriteData = 1;
end

% make feat directory if it doesn't exist from batch
featDir = [runInfo.anDir filesep 'feat'];
if overwriteData == 1
    if isdir(featDir)
        rmdir(featDir,'s')
    end
end
if ~isdir(featDir)
    mkdir([featDir filesep 'featureInfo']);
    mkdir([featDir filesep 'overlayImages']);
end

% look for region of interest info from project setup step
if ~exist([runInfo.anDir filesep 'roiMask_cropped.tif'])
    % not roi selected; use the whole image
    roiMask = ones(imL,imW);
    edgePix = [];
else
    % get roi edge pixels and make region outside mask NaN
    roiMask = double(imread([runInfo.anDir filesep 'roiMask_cropped.tif']));
    polyEdge = bwmorph(roiMask,'remove');
    edgePix = find(polyEdge);
end

% string for sigma in gauss filtered
s1 = length(num2str(max(gaussSigma)));
strg1 = sprintf('%%.%dd',s1);

% string for percent max cutoff
s2 = length(num2str(2));
strg2 = sprintf('%%.%dd',s2);

% string for number of files
s3 = length(num2str(length(startFrame:endFrame)));
strg3 = sprintf('%%.%dd',s3);


% START DETECTION

% initialize structure to store info for tracking
movieInfo(nFrames,1).xCoord = 0;
movieInfo(nFrames,1).yCoord = 0;
movieInfo(nFrames,1).amp = 0;
movieInfo(nFrames,1).int = 0;



% loop thru frames and detect
for i = startFrame:endFrame
    
    if i==startFrame
        tic
    end

    % load image and normalize to 0-1
    fileNameIm = [char(listOfImages(i,2)) filesep char(listOfImages(i,1))];
    img = double(imread(fileNameIm))./((2^bitDepth)-1);

    for sig = 1:length(gaussSigma) % loop thru sigma sizes

        % create kernels for gauss filtering
        blurKernelLow = fspecial('gaussian', 21, 1);
        blurKernelHigh = fspecial('gaussian', 21, gaussSigma(sig));

        % use subfunction that calls imfilter to take care of edge effects
        lowPass = filterRegion(img,roiMask,blurKernelLow);
        highPass = filterRegion(img,roiMask,blurKernelHigh);

        % get difference of gaussians image
        filterDiff = lowPass-highPass;

        % cut postive component of histogram; make below that nan
        filterDiff(filterDiff<0) = nan;
        [cutoffInd, cutOffValueInitInt] = cutFirstHistMode(filterDiff(:),1);
        noiseFloor = 1.0 * cutOffValueInitInt;
        filterDiff(filterDiff<noiseFloor) = nan;

        % get min/max of filterDiff
        minFiltDiff = nanmin(filterDiff(:));
        maxFiltDiff = nanmax(filterDiff(:));

        % we assume each step size down the intensity profile should be on
        % the order of the size of the noise level; here we find how many
        % steps we need and what their spacing should be 
        nSteps = round((maxFiltDiff-minFiltDiff)/(noiseFloor));
        threshList = linspace(maxFiltDiff,minFiltDiff,nSteps);

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
        featProp2 = regionprops(featMap2,filterDiff,'PixelIdxList','Area','MaxIntensity');
        [cutoffInd, cutoffValueArea] = cutFirstHistMode(vertcat(featProp2(:,1).Area),0);
        
        
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
        % we assume the good features have area > cutoffValueArea pixels
        % also good features have a max intensity > 2*cutoff
        goodFeatIdxA = find(vertcat(featProp2(:,1).Area)>=cutoffValueArea);
        goodFeatIdxI = find(vertcat(featProp2(:,1).MaxIntensity)>2*cutOffValueInitInt);
        goodFeatIdx = intersect(goodFeatIdxA,goodFeatIdxI);

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

        % centroid coordinates with zero uncertainties for Khuloud's tracker
        yCoord = zeros(nFeats,2); xCoord = zeros(nFeats,2);
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
        movieInfo(i,1).xCoord = xCoord;
        movieInfo(i,1).yCoord = yCoord;
        movieInfo(i,1).amp = amp;
        movieInfo(i,1).int = featI;

        if i==1
            tm = toc;
            expectedTime = (tm*nFrames)/60; % in minutes
            disp(['first frame took ' num2str(tm) ' seconds; ' num2str(nFrames) ' frames will take ' num2str(expectedTime) ' minutes.'])
        end
            indxStr1 = sprintf(strg1,gaussSigma(sig)); % gaussSigma
            %indxStr2 = sprintf(strg2,lowCutPercent);   % percentMaximum
            indxStr3 = sprintf(strg3,i);               % frame
            
            save([featDir filesep 'featureInfo' filesep 'featPropFinal' indxStr3],'featPropFinal');
            
        if savePlots==1
            % use extrema to draw polygon around feats - here we get
            % coordinates for polygon
%             outline = [featPropFinal.Extrema]; outline = [outline; outline(1,:)];
%             outlineX = outline(:,1:2:size(outline,2));
%             outlineY = outline(:,2:2:size(outline,2));

            

            %plot feat outlines and centroid on image
            figure(1);
            im2show = (lowPass-nanmin(lowPass(:)))./(nanmax(lowPass(:))-nanmin(lowPass(:))); % image for overlay
            im2show(edgePix) = 1; % ROI in white
            imshow(im2show,[],'Border','tight');
            hold on
            scatter(xCoord(:,1),yCoord(:,1),'c.'); % plot centroid in cyan
            %plot(outlineX,outlineY,'r'); % plot feat outlines in red
           
            saveas(gcf,[featDir filesep 'overlayImages' filesep 'overlay_g' indxStr1 '_f' indxStr3 '.tif']);
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

end

save([featDir filesep 'movieInfo'],'movieInfo');


close all
warning(warningState);


function filteredIm = filterRegion(im, mask, kernel)

im(mask~=1) = 0;
filteredIm = imfilter(im, kernel);
W = imfilter(double(mask), kernel);
filteredIm = filteredIm ./ W;
filteredIm(~mask) = nan;




