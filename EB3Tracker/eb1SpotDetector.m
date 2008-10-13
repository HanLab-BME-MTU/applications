function [movieInfo]=eb1SpotDetector(runInfo,nFrames,savePlots,debug,overwriteData,bitDepth)

% runInfo: structure containing image and analysis directories
% nFrames: number of frames on which to detect; []=all (default)
% debug: 0 to detect feats on all images; 1 to test one image

warningState=warning;
warning('off','MATLAB:divideByZero')

if nargin<5
    error('eb1SpotDetector: need 4 input arguments')
end

% check runInfo format and directory names between file systems
if ~isfield(runInfo,'imDir') || ~isfield(runInfo,'anDir')
    error('eb1SpotDetector: runInfo should contain fields imDir and anDir');
else
    [runInfo.anDir]=formatPath(runInfo.anDir);
    [runInfo.imDir]=formatPath(runInfo.imDir);
end

if bitDepth~=12 & bitDepth~=14 & bitDepth~=16
    warning('Choose 12, 14, or 16 for camera bit depth')
end


% make feat directory if it doesn't exist from batch
featDir=[runInfo.anDir filesep 'feat'];
if overwriteData==1
    if isdir(featDir)
        rmdir(featDir,'s')
    end
end
if ~isdir(featDir)
    mkdir([featDir filesep 'featMaps']);
    mkdir([featDir filesep 'overlay']);
end

% count number of images in image directory
[listOfImages] = searchFiles('.tif',[],runInfo.imDir,0);
nImTot=size(listOfImages,1);

if debug==1
    % let user select an image to test
    currentDir=pwd;
    cd(runInfo.imDir)
    [fileName,dirName] = uigetfile('*.tif','Choose a .tif file');
    cd(currentDir)

    for i=1:nImTot
        if strcmp(listOfImages{i,1},fileName)
            break
        end
    end
    m=i; n=i;

    gaussSigma=[2 5 10]; % sigma for gauss filtering

elseif debug==0
    % do detection on nFrames
    m=1;
    if ~isempty(nFrames) && nFrames>0 && nFrames<nImTot
        n=nFrames;
    else
        n=nImTot;
    end

    gaussSigma=4; % sigma for gauss filtering

end

% get image dimensions
[imL,imW]=size(double(imread([char(listOfImages(1,2)) filesep char(listOfImages(1,1))])));

if ~exist([runInfo.anDir filesep 'roiMask_cropped.tif'])
    % not roi selected; use the whole image
    roiMask=ones(imL,imW);
    edgePix=[];
else
    % get roi edge pixels and make region outside mask NaN
    roiMask=double(imread([runInfo.anDir filesep 'roiMask_cropped.tif']));
    polyEdge=bwmorph(roiMask,'remove');
    edgePix=find(polyEdge);
    %roiMask=swapMaskValues(roiMask,0,NaN);
end


% string for sigma in gauss filtered
s1=length(num2str(max(gaussSigma)));
strg1=sprintf('%%.%dd',s1);

% string for percent max cutoff
s2=length(num2str(2));
strg2=sprintf('%%.%dd',s2);

% string for number of files
s3=length(num2str(length(m:n)));
strg3=sprintf('%%.%dd',s3);


% initialize structure to store info for tracking
movieInfo(n-m+1,1).xCoord=0;
movieInfo(n-m+1,1).yCoord=0;
movieInfo(n-m+1,1).amp=0;
movieInfo(n-m+1,1).int=0;

% create kernels for gauss filtering
blurKernelLow=fspecial('gaussian', 21, 1);
blurKernelHigh=fspecial('gaussian', 21, 4);

if debug==1
    % load image and normalize to 0-1
    fileNameIm=[char(listOfImages(m,2)) filesep char(listOfImages(m,1))];
    img=double(imread(fileNameIm))./((2.^bitDepth)-1);

    minMax=zeros(length(gaussSigma),2);
    for sig=1:length(gaussSigma) % loop thru sigma sizes

        % use subfunction that calls imfilter to take care of edge effects
        lowPass=filterRegion(img,roiMask,blurKernelLow);
        highPass=filterRegion(img,roiMask,blurKernelHigh);

        % get difference of gaussians image
        filterDiff=roiMask.*(lowPass-highPass);

        % record min/max values for this sigma size
        minFiltDiff=nanmin(filterDiff(:));
        maxFiltDiff=nanmax(filterDiff(:));
        minMax(sig,:)=[minFiltDiff maxFiltDiff];

        % cut off negative values
        filterDiff(filterDiff<0)=0;

        % find gradient and magnitude
        [x y]=meshgrid(1:imW,1:imL);
        [px,py] = gradient(filterDiff,1,1);
        magF=sqrt(px.^2+py.^2);

        figure(sig); imshow(filterDiff,[]);
        hold on; quiver(x,y,px,py,'r')
        figure(2*sig); imshow(magF,[]);
    end
    figure(sig+1); imshow(lowPass,[])

else
    % not debug mode...loop thru frames and detect
    for i=m:n
if i==1
    tic
end
        % load image and normalize to 0-1
        fileNameIm=[char(listOfImages(i,2)) filesep char(listOfImages(i,1))];
        img=double(imread(fileNameIm))./((2.^bitDepth)-1);

        for sig=1:length(gaussSigma) % loop thru sigma sizes

            % use subfunction that calls imfilter to take care of edge effects
            lowPass=filterRegion(img,roiMask,blurKernelLow);
            highPass=filterRegion(img,roiMask,blurKernelHigh);

            % get difference of gaussians image
            filterDiff=lowPass-highPass;

            % cut postive component of histogram; make below that nan
            filterDiff(filterDiff<0)=nan;
            [cutoffInd, cutoffValue] = cutFirstHistMode(filterDiff,0);
            filterDiff(filterDiff<1.5*cutoffValue)=nan;                     % hardwired


            % find spots that bigger than 6 pixels in area
            featMap=bwlabel(filterDiff>0);
            featProp=regionprops(featMap,'PixelIdxList','Area');
            goodFeatIdx=find(vertcat(featProp(:,1).Area)>=8);               % hardwired

            % concat list of all pixels from good features
            % make new mask
            bwMask=nan*zeros(imL,imW);
            bwMask(vertcat(featProp(goodFeatIdx,1).PixelIdxList))=1;

            filterDiff=filterDiff.*bwMask;
            %figure(3); imshow(filterDiff,[]);                               % figure

            % get min/max of filterDiff
            minFiltDiff=nanmin(filterDiff(:));
            maxFiltDiff=nanmax(filterDiff(:));

            % we assume each step should be a near the inital
            % cutoffValue from the histogram
            nSteps=round((maxFiltDiff-minFiltDiff)/(1.0*cutoffValue));      % hardwired
            threshList=linspace(maxFiltDiff,minFiltDiff,nSteps);

            % compare features in z-slices startest from the highest one
            for p=1:length(threshList)-1

                % slice1 is top slice; slice2 is next slice down
                % here we generate BW masks of slices
                if p==1
                    slice1=filterDiff>threshList(p);
                else
                    slice1=slice2;
                end
                slice2=filterDiff>threshList(p+1);

                % now we label them
                featMap1=bwlabel(slice1);
                featMap2=bwlabel(slice2);
                featProp2=regionprops(featMap2,'PixelIdxList');

                % loop thru slice2 features and replace them if there are 2 or
                % more features from slice1 that contribute
                for iFeat=1:max(featMap2(:))
                    pixIdx=featProp2(iFeat,1).PixelIdxList; % pixel indices from slice2
                    featIdx=unique(featMap1(pixIdx)); % feature indices from slice1 using same pixels
                    featIdx(featIdx==0)=[]; % 0's shouldn't count since not feature
                    if length(featIdx)>1 % if two or more features contribute...
                        slice2(pixIdx)=slice1(pixIdx); % replace slice2 pixels with slice1 values
                    end
                end

            end

            % label slice2 again and get region properties
            featMap2=bwlabel(slice2);
            featProp2=regionprops(featMap2,filterDiff,'PixelIdxList','Area','MaxIntensity');

            % here we sort through features and retain only the "good" ones
            % we assume the good features have area > n pixels
            % also good features have a max intensity > 2*cutoff
            goodFeatIdxA=find(vertcat(featProp2(:,1).Area)>=8);             % hardwired
            goodFeatIdxI=find(vertcat(featProp2(:,1).MaxIntensity)>2*cutoffValue);
            goodFeatIdx=intersect(goodFeatIdxA,goodFeatIdxI);
            
            % make new label matrix and get props
            featureMap=zeros(imL,imW);
            featureMap(vertcat(featProp2(goodFeatIdx,1).PixelIdxList))=1;
            [featMapFinal,nFeats]=bwlabel(featureMap);
            featPropFinal=regionprops(featMapFinal,filterDiff,'PixelIdxList','Area','WeightedCentroid','Extrema','MaxIntensity');

%             % loop thru features and calc total intensity and peak pixel
%             for iFeat=1:nFeats
%                 % intensity values for pixels in feature iFeat
%                 featIntens=filterDiff(featPropFinal(iFeat,1).PixelIdxList);
% 
%                 % integrate intensity over feature
%                 featPropFinal(iFeat).IntSum=sum(featIntens);
% 
%                 % make centroid to be intensity maximum of feature
%                 maxIntPixIdx=featPropFinal(iFeat,1).PixelIdxList(featIntens==max(featIntens));
%                 [r c]= ind2sub([imL,imW],maxIntPixIdx);
%                 featPropFinal(iFeat).Centroid=[r c];
%                 
%             end
          
            % centroid coordinates with zero uncertainties for Khuloud's tracker
            yCoord=zeros(nFeats,2); xCoord=zeros(nFeats,2);
            temp=vertcat(featPropFinal.WeightedCentroid);
            yCoord(:,1)=temp(:,2);
            xCoord(:,1)=temp(:,1);

            % area
            featArea=vertcat(featPropFinal(:,1).Area);
            amp=zeros(nFeats,2);
            amp(:,1)=featArea;
            
            % intensity
            featInt=vertcat(featPropFinal(:,1).MaxIntensity);
            featI=zeros(nFeats,2);
            featI(:,1)=featInt;

            % make structure compatible with Khuloud's tracker
            movieInfo(i,1).xCoord=xCoord;
            movieInfo(i,1).yCoord=yCoord;
            movieInfo(i,1).amp=amp;
            movieInfo(i,1).int=featI;

if i==1
    toc
end

            if savePlots==1
                % use extrema to draw polygon around feats - here we get
                % coordinates for polygon
                outline=[featPropFinal.Extrema]; outline=[outline; outline(1,:)];
                outlineX=outline(:,1:2:size(outline,2));
                outlineY=outline(:,2:2:size(outline,2));

                indxStr1=sprintf(strg1,gaussSigma(sig)); % gaussSigma
                %indxStr2=sprintf(strg2,lowCutPercent);   % percentMaximum
                indxStr3=sprintf(strg3,i);               % frame

                %plot feat outlines and centroid on image
                figure(1);
                im2show=(lowPass-nanmin(lowPass(:)))./(nanmax(lowPass(:))-nanmin(lowPass(:))); % image for overlay
                im2show(edgePix)=1; % ROI in white
                imshow(im2show,[],'Border','tight');
                hold on
                scatter(xCoord(:,1),yCoord(:,1),'c.'); % plot centroid in cyan
                %plot(outlineX,outlineY,'r'); % plot feat outlines in red

                %saveas(gcf,[featDir filesep 'overlay' filesep 'overlay_g' indxStr1 '_f' indxStr3])
                saveas(gcf,[featDir filesep 'overlay' filesep 'overlay_g' indxStr1 '_f' indxStr3 '.tif']);

                % save image of colored feats
                RGB=label2rgb(featMapFinal, 'jet', 'k','shuffle');
                figure(2);
                imshow(RGB,'Border','tight');

                save([featDir filesep 'featMaps' filesep 'feats_g' indxStr1 '_f' indxStr3], 'featMapFinal');
                saveas(gcf,[featDir filesep 'featMaps' filesep 'feats_g' indxStr1 '_f' indxStr3 '.tif']);

            end

        end

    end
end

save([featDir filesep 'movieInfo'],'movieInfo');
close all
warning(warningState);


function filteredIm = filterRegion(im, mask, kernel)

im(mask~=1)=0;
filteredIm = imfilter(im, kernel);
W = imfilter(double(mask), kernel);
filteredIm = filteredIm ./ W;
filteredIm(~mask) = nan;




