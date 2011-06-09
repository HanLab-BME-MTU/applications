function [filterDiff,img,detectCoord,thresh] = plusTipTroubleShootDetect( projData,scales,bitDepth,useLastImg,multFactor4Thresh,doDetect,doJustFilterDiff, hfig, haxes)

% plusTipTroubleShootDetect: If doDetect = 1 and doJustFilterDiff ~= 1
%                            performs a detection on either the 
%                            first or last image using the set parameters
%                            1).Creates a grayscale figure of the original 
%                            unfiltered image with detected particles 
%                            overlayed in cyan. 
%                            2).Creates a 2-D surface contour figure of the difference of
%                            the diff gaussian filtered image applying the threshhold calculated from the 
%                            standard deviation of the image intensity.
%
%                            If doDetect ~= 1 and doJustFilterDiff ~= 1 
%                            will just create the 2-D surface contour figure
%                            of the diff gaussian filtered image (will not 
%                            perform detection). 
%
%                            If doDetect ~= 1 and doJustFilterDiff == 1,
%                            Just filter the image using the
%                            difference of the gaussian filter. Will 
%                            produce a imagesc fig of the filtered image 
%                            (this option is to be used when calling the 
%                             sigma parameter slider in the plusTipGetTracks GUI). 
%          
%                           
%                            
%                            
% %
% INPUT  projData        : structure containing fields .projData.anDir, which gives the full path 
%                          the full path to the roi_x directory
%                          and .imDir, which gives the full path to the
%                          folder containing the images for overlay.
%
%                          if given as [], program will query user for roi_x directory.
%        scales           : [low high] std for diff of Gaussians (input from slider)
%
%        bitDepth         : bit depth of the images - should be 12, 14, or 16
%
%        useLastImg       : if this value equals 1 use last image in movie else 
%
%        multFactor4Thresh: multiplication factor. Ultimate Thresh for detection 
%                           will be this factor*std of the image. Default = 3. 
%                           parameter to be manipulated by slider    
%
%        doDetect         : if doDetect == 1, run though the detection using 
%                           the set parameters and output a grayscale
%                           overlay figure
%
%        doJustFilterDiff : if doJustFilterDiff == 1, calculate just the 
%                           difference filter of either the first or the 
%                           last image, do not do detection or calculate
%                           the thresh
%   
%
% OUTPUT: 
% Calculated for all:
%       filterDiff : the difference of the gaussians filtered image (for either
%                    first or last of movie depending on useLastImg)
%
% (Only calculated if doJustFilterDiff ~= 1)
%       thresh     : the threshold for detection (calulated from the mean STD
%                    of the two images (for example if first image used
%                    calculate 
%
% (Only calculated if doJustFilterDiff ~= 1 AND doDetect == 1)
%       img        : the original image to be plotted for overlay
%       
%       detectCoord: coordinates for detected particles detectCoord(:,1) = x,
%                    detectCoord(:,2)= y
%
%
%                    
%       
%    
%

%% ---- Check Input and Load All Relevant Files ----
detectCoord = []; %%%
thresh = []; %%%

if nargin<8
   hfig = []; 
   haxes = [];
end

%GET IMAGE FILENAMES 
[listOfImages] = searchFiles('.tif',[],projData.imDir,0);
nImTot = size(listOfImages,1);
%Check if image file names are padded by zeros, If they are not 
% padded sort list to avoid error
% (Note there may be a better way to do this 
% quick fix by MB 07/2010)
imageName1 = [char(listOfImages(1,2)) filesep char(listOfImages(1,1))];
[path body no ext ] = getFilenameBody(imageName1);

if length(no)>1
    padded = 1;
else 
    padded = 0;
end

if padded == 0
    %Initialize Cells for Sorting
    %path = path of file, body = body of filename, ex/home/mb228/orchestra/groups/lccb-comet/Pellman/Mijung/7-28-10_EB_for_Maria/EB1_4_Cropped/roi_2t = extension of filename 
    %(tif etc) (all of these require a cell because they are strings)
    % num = number of filename (do not want in cell so can sort)
    pathCell = cell(nImTot,1);
    bodyCell = cell(nImTot,1);
    extCell = cell(nImTot,1);
    num = zeros(nImTot,1);

        %Sort List
        % For each frame get the image name from listOfImages
        for iFrame =  1:nImTot;
            imageName = [char(listOfImages(iFrame,2)) filesep char(listOfImages(iFrame,1))];

            %Call "getFilenameBody" from common dir to split filename into path,
            %body, no, and ext. Put path, body and ext into appropriate cell/number vector for that
            %frame
            [path body no ext ] = getFilenameBody(imageName);


            pathCell(iFrame) = cellstr(path);
            bodyCell(iFrame) = cellstr(body);
            extCell(iFrame) = cellstr(ext);

            % output "no" is a string so convert to number to sort
            num(iFrame)  = str2double(no);
 
        end
    %Sort number vector numerically
    sortednum = sort(num);

    %Convert number vector to cell
    sortednum_cell = num2cell(sortednum);

    %Create Sorted Image List
    sortedImages = [pathCell, bodyCell, sortednum_cell, extCell];
else
end % end sort image list 
 
%LOAD FIRST IMAGE FILE

fileNameIm = [char(listOfImages(1,2)) filesep char(listOfImages(1,1))];
img = double(imread(fileNameIm))./((2^bitDepth)-1);
[imL,imW] = size(img);


%[path body no ext ] = getFilenameBody(projData.imDir);
%projData.anDir = [path filesep 'roi_2'];

%LOAD MASK 
% look for region of interest info from project setup step
if ~exist([projData.anDir filesep 'roiMask.tif'],'file')
    % not roi selected; use the whole image
    roiMask = ones(imL,imW);
    roiYX=[1 1; imL 1; imL imW; 1 imW; 1 1];
else
    % get roi edge pixels and make region outside mask NaN
    roiMask = double(imread([projData.anDir filesep 'roiMask.tif']));
    
    % for some reason the donut mask will open not as 1 and 0s but as 
    % 255 (an 8 bit image) and 0s, Check for this and correct so the donut mask will 
    % run through detection (added MB 09/2010)
    if isempty(find(roiMask==1,1))
        roiMask(roiMask==255) = 1;
    else
    end 
    roiYX=load([projData.anDir filesep 'roiYX']);
    roiYX=roiYX.roiYX;
end


% CHECK ALL OTHER INPUT 

% check input for sigma values of gaussians 
if nargin<2 || isempty(scales)
    scales=[1 4]; 
end   
    
% get bit depth if not given
if nargin < 3 || isempty(bitDepth)
    imgData = imfinfo(fileNameIm);
    bitDepth = imgData.BitDepth;
    disp(['bitDepth estimated to be' bitDepth])
end

% check bit depth to make sure it is 12, 14, or 16 and that its dynamic
% range is not greater than the provided bitDepth
img = double(imread(fileNameIm));
[imL,imW] = size(img);
maxIntensity = max(img(:));

if sum(bitDepth==[12 14 16])~=1 || maxIntensity > 2^bitDepth-1
    error('--plusTipCometDetector: bit depth should be 12, 14, or 16');
end


% check input for useLastImg
if nargin<4 || isempty(useLastImg)
   useLastImg = 0; % default is to use the first image
end

% check input for multFactor4Thresh
if nargin<5 || isempty(multFactor4Thresh)
    multFactor4Thresh = 3; 
end 
    
% check input for doDetect
if nargin<6 || isempty(doDetect)
    doDetect = 1; % default is to do the detection
end 


% check input for doJustFilterDiff
if nargin<7 || isempty(doJustFilterDiff)
    doJustFilterDiff = 0; % default is to run through the detection
end 
 
    
%% ----Initialize ----

% Initiate frame list, if only do filterDiff will only need to filter one
% image, if proceed with thresh calc/detection need two images to get STD 
% for step size and thresh 

if doJustFilterDiff == 1 ;
    frameList = zeros(1,1); % need to filter one image
    stdList = zeros(1,1) ;
else  
    
    frameList = zeros(2,1); % need to filter two images
    stdList = zeros(2,1); % will need to calc STD for both images
end

%if useFirstImg option selected use the first image of list, else use
% the last image of list 

if useLastImg == 1 % use 2nd to last and last image of movie
    
    if doJustFilterDiff == 1
            frameList(1,1) = size(listOfImages,1); % filter just the last image
        else % set up for thresh calc (need two filter diff images)
            frameList(2,1) =   size(listOfImages,1); % final image filtered
            % will be last image (so can go on plotting without saving)
            frameList(1,1)  =  frameList(2,1) - 1 ; % next to last image in movie
    end 

        
else % use 1st and 2nd image of movie
    if doJustFilterDiff  == 1
       frameList(1,1) = 1; % filter just the first image
        else % set up for thresh calc (need two filter diff images)
        frameList(1,1) = 2;
        frameList(2,1) = 1; %last image filtered will be 1st image (so can go
        % on to plot without saving)
    end % end if firstImg
end % if use last image

%Initiate iteration counter   
    iter = 1; 

%Name for figures 
if useLastImg == 1
    name = 'Last Image of Movie';
else 
    name = 'First Image of Movie';
end 

%% ---- MAIN BODY ---- %%
 
for i = 1:length(frameList)
    
    iFrame = frameList(i,1); % get frame of interest
    
    if padded == 1;% If image filename padded use listOfImages 
            fileNameIm = [char(listOfImages(iFrame,2)) filesep char(listOfImages((iFrame),1))];
        else % use sortedImages 
            fileNameIm = [char(sortedImages(iFrame,1)) filesep char(sortedImages(iFrame,2)),...;
            num2str(sortednum(iFrame)) char(sortedImages(iFrame,4))];
    end
    
    % Load and normalize image based on bitdepth
    img = double(imread(fileNameIm))./((2^bitDepth)-1); 
    
    
    %FILTER IMAGE

    % create kernels for gauss filtering
    blurKernelLow  = fspecial('gaussian', 21, scales(1));
    blurKernelHigh = fspecial('gaussian', 21, scales(2));

    % use internal subfunction that calls imfilter to take care of edge effects
    
    lowPass = filterRegion(img,roiMask,blurKernelLow); 
    highPass = filterRegion(img,roiMask,blurKernelHigh);

    % get difference of gaussians image
    filterDiff = lowPass-highPass;
    
    % if bg point was chosen and saved, get bgMask from first frame
    if i==1 && exist([projData.anDir filesep 'bgPtYX.mat'])~=0
        bgPtYX=load([projData.anDir filesep 'bgPtYX.mat']);
        bgPtYX=bgPtYX.bgPtYX;
        [bgMask]=eb3BgMask(filterDiff,bgPtYX);
        %saveas(gcf,[featDir filesep 'filterDiff' filesep 'bgMask.tif']);
        close(gcf)
    end
    
    % if bg point wasn't chosen, use ROI
    if i==1 && exist([projData.anDir filesep 'bgPtYX.mat'])==0
        bgMask=logical(roiMask); %Note: not sure why she has logical here (it doesn't change anything as far as I can tell)
    end
    
    %CALC: std of image intensity
    stdList(iter) = std(filterDiff(bgMask));
   
    iter = iter + 1;
end  % end for i
    
if doJustFilterDiff == 1 % if just want to do filter diff stop here and plot
    
    %CREATE IMAGESC FIGURE OF FILTERED IMAGE    
    himage = [];
    if isempty(haxes)
        figure
    else
        set(hfig, 'CurrentAxes', haxes)
        limits = get(haxes,{'XLim','YLim'});
        
        if strcmp(get(get(haxes,'Children'),'Type'),'image');
            himage=get(haxes,'Children');
        end
    end 
    
   
    if ~isempty(himage) && ishandle(himage)
        set(himage,'CData',filterDiff)
        set(haxes,{'XLim','YLim'},limits);
    else
        imagesc(filterDiff);
        axis([1,imW,1,imL])
    end

    colorbar; 
    title(name);
    
else % Proceed with Calc
        % Calculate step size/thresh based on standard deviation of image
        % intensity
        stepSize = mean(stdList(1:2));       
        thresh = multFactor4Thresh*stepSize;

if doDetect == 1 
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

        % now we label them using the "bwlabel" function from matlab which 
        % labels connected components in a 2-D binary image
        featMap1 = bwlabel(slice1);
        featMap2 = bwlabel(slice2);
        
        % get the regionproperty 'PixelIdxList' using "regionprops" function in matlab 
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
    featProp2 = regionprops(featMap2,'PixelIdxList','Area');

    % here we sort through features and retain only the "good" ones
    % we assume the good features have area > 2 pixels
    goodFeatIdx = find(vertcat(featProp2(:,1).Area)>2);
%    goodFeatIdxI = find(vertcat(featProp2(:,1).MaxIntensity)>2*cutOffValueInitInt);
%    goodFeatIdx = intersect(goodFeatIdxA,goodFeatIdxI);

    % make new label matrix and get props
    featureMap = zeros(imL,imW);
    featureMap(vertcat(featProp2(goodFeatIdx,1).PixelIdxList)) = 1;
    [featMapFinal,nFeats] = bwlabel(featureMap);
    
    verDate=version('-date');
    if str2double(verDate(end-3:end))>=2008
        featPropFinal = regionprops(featMapFinal,filterDiff,'PixelIdxList','Area','WeightedCentroid','MaxIntensity'); %'Extrema'
    else
        featPropFinal = regionprops(featMapFinal,'PixelIdxList','Area','Centroid');
        for iFeat=1:length(featPropFinal)
            featPropFinal(iFeat,1).WeightedCentroid=featPropFinal(iFeat,1).Centroid; % centroid's close enough...
            featPropFinal(iFeat,1).MaxIntensity=max(filterDiff(featPropFinal(iFeat,1).PixelIdxList)); % find maximum intensity
        end
    end

    if nFeats==0
        yCoord = [];
        xCoord = [];
        amp = [];
        featI = [];
        
    else
        % centroid coordinates with 0.5 uncertainties for Khuloud's tracker
        yCoord = 0.5*ones(nFeats,2);
        xCoord = 0.5*ones(nFeats,2);
        temp = vertcat(featPropFinal.WeightedCentroid);
        yCoord(:,1) = temp(:,2);
        xCoord(:,1) = temp(:,1);
        detectCoord = temp;
        
        % area
        featArea = vertcat(featPropFinal(:,1).Area);
        amp = zeros(nFeats,2);
        amp(:,1) = featArea;

        % intensity
        featInt = vertcat(featPropFinal(:,1).MaxIntensity);
        featI = zeros(nFeats,2);
        featI(:,1) = featInt;
    end

    



    %plot feat outlines and centroid on image
    if useLastImg == 1 
        plotFrame = size(listOfImages,1); 
    else 
        plotFrame = 1;
    end
        
        if padded == 1;% If image filename padded use listOfImages 
            fileNameIm = [char(listOfImages(plotFrame,2)) filesep char(listOfImages(plotFrame,1))];
        else % use sortedImages 
            fileNameIm = [char(sortedImages(plotFrame,1)) filesep char(sortedImages(plotFrame,2)),...;
            num2str(sortednum(plotFrame)), char(sortedImages(plotFrame,4))];
        end
        
        %MAKE OVERLAY PLOT OF ORIGNAL IMAGE
        img = double(imread(fileNameIm))./((2^bitDepth)-1);
        figure
        imagesc(img)
        hold on
        scatter(xCoord(:,1),yCoord(:,1),'c.'); % plot centroid in cyan
        colormap gray
        plot(roiYX(2),roiYX(1),'w')
        axis equal   
        title(name);

        %MAKE SURF CONTOUR FIG OF DIFF OF GAUSSIAN 
        hold off
        forTitle = ['Thresh = ', num2str(thresh)];
        figure;
        surf(filterDiff)
        colormap;
        axis([1,imL,1,imW,thresh,max(filterDiff(:))])
        view(2)
        set(gca,'YDir','reverse');
        title({name ; forTitle});
        colorbar;
   
    
    else  % skip detection and just make surf contour figure
        forTitle = ['Thresh = ', num2str(thresh)];
        if isempty(haxes)
            figure('Toolbar', 'figure');

        else
            set(hfig, 'CurrentAxes', haxes)
            limits = get(haxes,{'XLim','YLim'});
            % For uniform color bar in preview GUI
            minvalue = min(filterDiff(:)); 
            filterDiff(filterDiff(:)<thresh) = NaN;
            filterDiff(1,1) = minvalue;
        end
        
        surf(filterDiff)
        axis([1,imW,1,imL,thresh,max(filterDiff(:))])
        if isempty(haxes)
            %shading flat            
        else
            %shading(haxes, 'flat')
            set(haxes, 'YAxisLocation', 'right')
            set(haxes,{'XLim','YLim'},limits);
        end
        
        
        view(2)
        set(gca,'YDir','reverse');
        title(forTitle);
        colorbar;
        
end % end if doDetect == 1  
end % end if doJustFilterDiff == 1 
    
end

% internal function called above
function filteredIm = filterRegion(im, mask, kernel)

im(mask~=1) = 0;
filteredIm = imfilter(im, kernel);
W = imfilter(double(mask), kernel);
filteredIm = filteredIm ./ W;
filteredIm(~mask) = nan;

end 

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


figure 
imagesc(filterDiff); 
colormap gray;
axis equal
hold on
scatter(bgPtYX(2),bgPtYX(1),'*y') % user-selected point
scatter(c(goodIdx),r(goodIdx),'.g') % all "good" features in green
scatter(c(goodIdx(closeIdx)),r(goodIdx(closeIdx)),'r') % nearest fifty to point in red
plot(xi,yi) % plot mask outline in blue
end
