function [frameInfo] = spotDetection(sourcePath, postProcLevel)
%
% Performs detection of local intensity clusters through a combination of 
% multiscale products and denoising by iterative filtering from
% significant coefficients:
% Olivo-Marin, "Extraction of spots in biological images using multiscale products," Pattern Recoginition 35, pp. 1989-1996, 2002.
% Starck et al., "Image Processing and Data Analysis," Section 2.3.4, p. 73
%
% INPUTS:   sourcePath               : path to TIF frame directory
%           postProcLevel (optional) : postprocessing level. 

% Parts of this function are based on code by Henry Jaqaman.
% Francois Aguet, March 2010

%===================================================
% Process inputs
%===================================================
if nargin<1 || isempty(sourcePath)
    sourcePath = [uigetdir('.tif', 'Select directory containing the movie frames') filesep];
end
if nargin<2 || isempty(postProcLevel)
   postProcLevel = 0; 
end

tifFiles = dir([sourcePath '*.tif*']);
nFrames = length(tifFiles);

detectionpath = [sourcePath 'Detection' filesep];
if ~(exist(detectionpath, 'dir')==7)
    mkdir(detectionpath);
end

maskPath = [sourcePath 'Detection' filesep 'Masks' filesep];
if ~(exist(maskPath, 'dir')==7)
    mkdir(maskPath);
end

dthreshold = 5;
S = 4; % max wavelet scale

%===================================================
% Process frames
%===================================================
frameInfo(1:nFrames) = struct('xmax', [], 'ymax', [], 'xcom', [], 'ycom', [],...
                              'totalInt', [], 'area', [], 'nMaxima', [], 'labels', [],...
                              'nComp', [], 'amp', [], 'xCoord', [], 'yCoord', []);

fprintf('Progress:     ');
for k = 1:nFrames
    %detection1 = double(imread([sourcePath 'images283' filesep 'cmask3001.jpg']));
    img = double(imread([sourcePath tifFiles(k).name]));
    
    maxI = max(img(:));
    minI = min(img(:));
    img = img*450/max(img(:));
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
    imgDenoised = imgDenoised - min(imgDenoised(:)); %%%%% This generates false detections, will be removed in next release.
        
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
    imgDenoised = imgDenoised * (maxI-minI) / max(imgDenoised(:));

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
    xmax2 = cell(nComp, 1);
    ymax2 = cell(nComp, 1);
    area2 = cell(nComp, 1);
    totalInt2 = cell(nComp, 1);
    labels2 = cell(nComp, 1);
        
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
        elseif nMaxima(n)==0 % no maximum was detected for this cluster
            maxValueIdx = find(values == max(values));
            xmax(n) = xi(maxValueIdx(1));
            ymax(n) = yi(maxValueIdx(1));
            nMaxima(n) = 1;
        else % resolve multiple maxima cases
            maxValues = localMax(sub2ind(size(localMax), ym, xm)); % highest local max
            maxIdx = find(maxValues == max(maxValues));
            xmax(n) = xm(maxIdx(1));
            ymax(n) = ym(maxIdx(1));
            
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
                labels2{n} = labels(n)*ones(nSecMax,1);
            
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
    
    % assign
    frameInfo(k).xmax = [xmax; xmax2(:)];
    frameInfo(k).ymax = [ymax; ymax2(:)];
    frameInfo(k).xcom = [xcom; xmax2(:)];
    frameInfo(k).ycom = [ycom; ymax2(:)];
    frameInfo(k).totalInt = [totalInt; totalInt2(:)];
    frameInfo(k).area = [area; area2(:)];
    
    frameInfo(k).nMaxima = nMaxima; % maxima per component
    %frameInfo(k).labels = [labels; [labels2{:}]]; % labels MxN -> N
    frameInfo(k).nComp = nComp;
    
    
    
    % prepare fields for tracker
    nObj = length(frameInfo(k).xmax);
    frameInfo(k).amp = zeros(nObj,2);
    frameInfo(k).xCoord = zeros(nObj,2);
    frameInfo(k).yCoord = zeros(nObj,2);
    
    frameInfo(k).amp(:,1) = frameInfo(k).totalInt;
    frameInfo(k).xCoord(:,1) = frameInfo(k).xcom;
    frameInfo(k).yCoord(:,1) = frameInfo(k).ycom;
    
    %fprintf('Frame #%d, #maxima = %d, #components = %d\n', k, nMaxima, nComp);
    
    % Save denoised and masked frame as compressed TIFF file
    imwrite(uint8(255*imgDenoised/max(imgDenoised(:))), [maskPath 'dmask_' num2str(k, ['%.' num2str(length(num2str(nFrames))) 'd']) '.tif'], 'tif', 'compression' , 'lzw');
    fprintf('\b\b\b\b%3d%%', round(100*k/nFrames));
end
fprintf('\n');
save([sourcePath 'Detection' filesep 'detectionResults.mat'], 'frameInfo');

% figure;
% xa = 1:max([frameInfo(:).area]);
% na = hist([frameInfo(:).area], xa);
% h = bar(xa, na);
% set(h, 'FaceColor', [0 0.8 0], 'EdgeColor', [0 0.4 0], 'BarWidth', 1);
% title('Area of detected spots (all frames)', 'FontName', 'Helvetica', 'FontSize', 14);

fprintf('Detection complete.\n');



function result = significantCoefficientDenoising(img, S)
% Former implementation
mask = zeros(size(img));
result = zeros(size(img));
W = awt(img, S);
for s = 1:S
    tmp = W(:,:,s);
    mask(abs(tmp) >= 3*std(tmp(:))) = 1;
    result = result + tmp.*mask;
end

% Correct version
% mask = zeros([size(img) S]);
% W = awt(img, S);
% for s = 1:S
%     tmp = W(:,:,s);
%     mask(:,:,s) = abs(tmp) >= 3*mad(tmp(:),1)/0.6745;
% end
% result = sum(mask.*W(:,:,1:S),3);

