function [frameInfo imgDenoised] = detectSpotsWT(img, S, dthreshold, postProcLevel)

if nargin<2
    S = 4;
end
if nargin<3
    dthreshold = 5;
end
if nargin<4
    postProcLevel = 1;
end


maxI = max(img(:));
minI = min(img(:));
%img = img*450/max(img(:));
[ny nx] = size(img);

%===================================================
% Iterative filtering from significant coefficients
%===================================================

%imgDenoised = significantCoefficientDenoising(img, S);
%imgDenoised = awtDenoising(img,S,0);
imgDenoised = wstarck2(img, S, 0);




res = img - imgDenoised; % residuals
sigma_res0 = std(res(:));

delta = 1;
while delta > 0.002
    %resDenoised = significantCoefficientDenoising(res, S);
    %resDenoised = awtDenoising(res,S,0);
    resDenoised = wstarck2(res, S, 0);
    imgDenoised = imgDenoised + resDenoised; % add significant residuals
    res = img - imgDenoised;
    sigma_res1 = std(res(:));
    delta = abs(sigma_res0/sigma_res1 - 1);
    sigma_res0 = sigma_res1;
end

%===================================================
% Multiscale product of wavelet coefficients
%===================================================
% The support of the objects is given by the multiscale product in the wavelet domain.
% Since only positive objects are of interest, negative coefficients are set to zero.
W = awt(imgDenoised, S);
%imgDenoised = imgDenoised - min(imgDenoised(:)); %%%%% This generates false detections, will be removed in next release.

% W(W<0) = 0;
% imgMSP = abs(prod(W(:,:,1:S),3));
imgMSP = wstarck222(imgDenoised, S, 1);

%===================================================
% Binary mask
%===================================================
% Establish thresholds
%[imAvg imStd] = localAvgStd2D(imgDenoised, 9);
[imAvg imStd] = wlocav(imgDenoised, 4); %%%%%%%%%%%%%%%%%%%%%%%%%

mask = zeros(ny,nx);
mask((imgDenoised >= imAvg+0.5*imStd) & (imgDenoised.*imgMSP >= mean(imgDenoised(:)))) = 1;



% Morphological postprocessing
mask = bwmorph(mask, 'clean'); % remove isolated pixels
mask = bwmorph(mask, 'fill'); % fill isolated holes
mask = bwmorph(mask, 'thicken');
mask = bwmorph(mask, 'spur'); % remove single pixels 8-attached to clusters
mask = bwmorph(mask, 'spur');
mask = bwmorph(mask, 'clean');

if postProcLevel >= 1
    mask = bwmorph(mask, 'erode');
    if postProcLevel == 2
        mask = bwmorph(mask, 'spur');
    end
    mask = bwmorph(mask, 'clean');
    mask = bwmorph(mask, 'thicken');
end



%     mask = imgMSP > 0;
%
%     % Clean mask
%     mask = bwmorph(mask, 'clean'); % remove isolated single pixels
%     mask = bwmorph(mask, 'fill'); % fill isolated 1-pixel holes
%     mask = imopen(mask, strel('square', 2));
%     mask = bwmorph(mask, 'thicken', 1);
%     %mask = imopen(mask, strel('disk', 1));



% rescale denoised image
%imgDenoised = imgDenoised * (maxI-minI) / max(imgDenoised(:));
imgDenoised = (imgDenoised-min(imgDenoised(:))) * (maxI-minI) / (max(imgDenoised(:))-min(imgDenoised(:)));


imgDenoised = mask.*imgDenoised;
localMax = locmax2d(imgDenoised, [9 9]);


%     rgbDisplay = zeros([size(mask) 3]);
%     c1 = scaleContrast(img);
%     c2 = c1;
%     c3 = c1;
%     c3(mask~=0 | detection1~=0) = 0;
%     xIdx = mask~=0 & detection1~=0;
%     oIdx = detection1~=0;
%     oIdx = xor(xIdx,oIdx);
%     nIdx = mask~=0;
%     nIdx = xor(xIdx,nIdx);
%     c1(nIdx) = 0;
%     c2(oIdx) = 0;
%
%     rgbDisplay(:,:,1) = c1;
%     rgbDisplay(:,:,2) = c2;
%     rgbDisplay(:,:,3) = c3;
%     figure;
%     ax(1) = subplot(1,2,1);
%     imagesc(img); colormap(gray(256)); axis image;
%     ax(2) = subplot(1,2,2);
%     imagesc(uint8(rgbDisplay)); colormap(gray(256)); axis image;
%     linkaxes(ax,'xy');
%     title('New vs. old detection');
%
%     maskRef = detection1;
%     maskRef(maskRef~=0) = 1;
%     figure; imagesc(mask-maskRef); colormap(gray(256)); axis image; colorbar;
%     sum(sum(abs(mask-maskRef)))

%===================================================
% Process connected components
%===================================================

[labels, nComp] = bwlabel(mask, 8);

area = zeros(nComp, 1);
totalInt = zeros(nComp, 1);
nMaxima = zeros(nComp, 1);
xmax = zeros(nComp, 1);
ymax = zeros(nComp, 1);
xcom = zeros(nComp, 1);
ycom = zeros(nComp, 1);
labelVect = zeros(nComp, 1);

xmax2 = cell(nComp, 1);
ymax2 = cell(nComp, 1);
area2 = cell(nComp, 1);
totalInt2 = cell(nComp, 1);
labelVect2 = cell(nComp, 1);

% Compute area and center of mass for each component
stats = regionprops(labels, imgDenoised, 'Area', 'WeightedCentroid', 'PixelIdxList');

% component labels of local maxima
maxLabels = labels .* (labels & localMax>0);
maxCoords(1:nComp) = struct('PixelIdxList', []);
mc = regionprops(maxLabels, 'PixelIdxList');
maxCoords(1:length(mc)) = deal(mc);


for n = 1:nComp
    %[yi,xi] = find(labels == n); % coordinates of nth component
    [yi,xi] = ind2sub([ny nx], stats(n).PixelIdxList);
    [ym,xm] = ind2sub([ny nx], maxCoords(n).PixelIdxList);
    area(n) = stats(n).Area;
    com = stats(n).WeightedCentroid;
    xcom(n) = com(1);
    ycom(n) = com(2);
    
    values = imgDenoised(stats(n).PixelIdxList);
    totalInt(n) = sum(values);
    
    nMaxima(n) = length(xm);
    %nMaxima(n)
    if nMaxima(n)==1
        xmax(n) = xm;
        ymax(n) = ym;
        nMaxima(n) = 1;
        labelVect(n) = labels(ym,xm);
    elseif nMaxima(n)==0 % no maximum was detected for this cluster
        maxValueIdx = find(values == max(values));
        xmax(n) = xi(maxValueIdx(1));
        ymax(n) = yi(maxValueIdx(1));
        nMaxima(n) = 1;
        labelVect(n) = labels(ymax(n), xmax(n));
    else % resolve multiple maxima cases
        maxValues = localMax(sub2ind(size(localMax), ym, xm)); % highest local max
        maxIdx = find(maxValues == max(maxValues));
        xmax(n) = xm(maxIdx(1));
        ymax(n) = ym(maxIdx(1));
        labelVect(n) = labels(ymax(n), xmax(n));
        
        % remove highest max from list
        xm(maxIdx(1)) = [];
        ym(maxIdx(1)) = [];
        
        % compute distance of secondary maxima to primary
        dist2max = sqrt((xmax(n)-xm).^2 + (ymax(n)-ym).^2);
        dist2com = sqrt((xcom(n)-xm).^2 + (ycom(n)-ym).^2);
        mindist = min(dist2max,dist2com);
        
        % retain secondary maxima where mindist > threshold
        idx2 = find(mindist > dthreshold);
        if ~isempty(idx2)
            xmax2{n} = xm(idx2);
            ymax2{n} = ym(idx2);
            nSecMax = length(idx2);
            nMaxima(n) = nSecMax+1;
            
            % split area
            area2{n} = area(n)*ones(nSecMax,1)/nMaxima(n);
            area(n) = area(n)/nMaxima(n);
            labelVect2{n} = labels(sub2ind(size(labels), ymax2{n}, xmax2{n}));
            
            %intensity values
            totalInt2{n} = totalInt(n)*ones(nSecMax,1)/nMaxima(n);
            totalInt(n) = totalInt(n)/nMaxima(n);
        end
    end
end

xmax2 =  vertcat(xmax2{:});
ymax2 = vertcat(ymax2{:});
totalInt2 = vertcat(totalInt2{:});
area2 = vertcat(area2{:});
labelVect2 = vertcat(labelVect2{:});

% assign
frameInfo.xmax = [xmax; xmax2(:)];
frameInfo.ymax = [ymax; ymax2(:)];
frameInfo.xcom = [xcom; xmax2(:)];
frameInfo.ycom = [ycom; ymax2(:)];
frameInfo.totalInt = [totalInt; totalInt2(:)];
frameInfo.area = [area; area2(:)];

frameInfo.nMaxima = nMaxima; % maxima per component
frameInfo.labels = [labelVect; labelVect2(:)];
frameInfo.nComp = nComp;



% prepare fields for tracker
nObj = length(frameInfo.xmax);
frameInfo.amp = zeros(nObj,2);
frameInfo.xCoord = zeros(nObj,2);
frameInfo.yCoord = zeros(nObj,2);

frameInfo.amp(:,1) = frameInfo.totalInt;
frameInfo.xCoord(:,1) = frameInfo.xcom;
frameInfo.yCoord(:,1) = frameInfo.ycom;

frameInfo.path = [];
frameInfo.maskPath = [];

%fprintf('Frame #%d, #maxima = %d, #components = %d\n', k, nMaxima, nComp);

% Save denoised and masked frame as compressed TIFF file
% imwrite(uint8(255*imgDenoised/max(imgDenoised(:))), frameInfo.maskPath, 'tif', 'compression' , 'lzw');
