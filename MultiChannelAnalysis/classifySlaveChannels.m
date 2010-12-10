% Francois Aguet, October 2010

function tracks = classifySlaveChannels(data, tracks)

% intensity-based classification of slave channels alone is unreliable
% measure: correlation (time-based max)
% compare with spline-based detection on Allen's data

% procedure:
% 1) Randomly selected in each frame are linked as phantom tracks


% Approach A.) independent of master channel
% 1) use track + mask ->
% decision frame by frame
% level above background: x * std(background)
% distrib. of x -> cutoff

% Significance thresholds
alpha = 0.05;
% sigmaT = icdf('normal', 1-alpha/2, 0, 1);
sigmaT = icdf('normal', 1-alpha, 0, 1); % weaker, single-tailed


%=================================
% Determine master/slave channels
%=================================
% detect number of channels (up to 4)
nChannels = length(data.channels);
% exclude master from list of channels
masterChannel = regexp(data.source, data.channels);
masterChannel = find([masterChannel{:}]);
slaveChannels = setdiff(1:nChannels, masterChannel);

sigma = getGaussianPSFsigma(data.NA, data.M, data.pixelSize, data.markers{masterChannel});
w = ceil(3*sigma);

%===============================================
% Classification
%===============================================


%-----------------------------------------------
% Utilization of a generic mask is not sufficient:
%-----------------------------------------------
% define mask for localization
bgMask = ones(2*w+1);
[xg,yg] = meshgrid(-w:w);
r = sqrt(xg.^2+yg.^2);
bgMask(r<=2*sigma) = 0;
% figure; imagesc(bgMask); colormap(gray(256)); axis image; colorbar;




%-----------------------------------------------
% 1. Generate mask for background points
%-----------------------------------------------
maskList = dir([data.source 'Detection' filesep 'Masks' filesep '*.tif']);

for k = 1:10:data.movieLength

    % load mask, dilate, -mask
    mask = double(imread([data.source 'Detection' filesep 'Masks' filesep maskList(k).name]));
    mask(mask~=0) = 1;
    %dmask = imdilate(mask, strel('disk', 2*w)) - imdilate(mask, strel('disk', w));
    dmask = imdilate(mask, strel('disk', 1*w));
    dmask(dmask~=0) = 1;
    %dmask = dmask - mask;
    %dmask(dmask<0) = 0;

    % generate background points
    % ny = data.imagesize(1);
    % nx = data.imagesize(2);
    
    
    % % in first frame, determine proportion of positive detections in background
    %
    % % Mask excluding CCPs in current frame
    % mask = double(imread([data.source 'Detection' filesep 'Masks' filesep maskList(10).name]));
    % mask(mask~=0) = 1;
    % mask = maskInt-imdilate(mask, strel('disk', w));
    % mask(mask<0) = 0;
    % % figure; imagesc(mask); colormap(gray(256)); axis image;
    %
    %
    % npRef = 3000;
    % % generate candidate points
    % x = (nx-2*w-1)*rand(1,npRef)+w+1;
    % y = (ny-2*w-1)*rand(1,npRef)+w+1;
    % xi = round(x);
    % yi = round(y);
    %
    % % remove points outside of mask or within border
    % idx = dmask(sub2ind([ny nx], yi,xi))==0 | xi<=w | yi<=w | xi>nx-w | yi>ny-w;
    % x(idx) = [];
    % y(idx) = [];
    % xi(idx) = [];
    % yi(idx) = [];
    % npRef = length(x);
    
    
    %
    % frameList = dir([data.channels{2} '*.tif*']);
    % frame = double(imread([data.channels{2} frameList(1).name]));
    % ctrlStatus = zeros(1, npRef);
    % for p = 1:npRef
    %     window = frame(yi(p)-w:yi(p)+w, xi(p)-w:xi(p)+w);
    %     c = mean(window(lmask==1));
    %     cStd = std(window(lmask==1));
    %     prm = fitGaussian2D(window, [x(p)-xi(p) y(p)-yi(p) max(window(:))-c sigma c], 'xyA');
    %     A = prm(3);
    %     ctrlStatus(p) = A > sigmaT*cStd;
    % end
    % fprintf('Ch. %d background detections: %d/%d (%.2f %%)\n', 2, sum(ctrlStatus), npRef, 100*sum(ctrlStatus)/npRef);
    
    
    %-----------------------------------------------
    % Pick random mask from detection
    %-----------------------------------------------
    load([data.source 'Detection' filesep 'detectionResults.mat']);
    
    frameList = dir([data.channels{2} '*.tif*']);
    frame = double(imread([data.channels{2} frameList(10).name]));
    %mask = double(imread([data.source 'Detection' filesep 'Masks' filesep maskList(k).name]));
    
    
    % % use mask and localized coordinates on background
    % xm = frameInfo(1).xloc;
    % ym = frameInfo(1).yloc;
    % xmi = round(xm);
    % ymi = round(ym);
    % idx = xmi<=w | ymi<=w | xmi>nx-w | ymi>ny-w | isnan(xmi);
    % xm(idx) = [];
    % ym(idx) = [];
    % xmi(idx) = [];
    % ymi(idx) = [];
    
    % ctrlStatus = zeros(1, npRef);
    N = sum(bgMask(:));
    
    % background at each pixel
    E = conv2(padarray(frame, [w w], 'replicate'), bgMask/N, 'valid');
    E2 = conv2(padarray((frame).^2, [w w], 'replicate'), bgMask/N, 'valid');
    cStdMap = sqrt(N/(N-1) * (E2-E.^2));
    cMap = E;
    
    % amplitude map
    gL2 = pi*sigma^2;
    g = exp(-(-w:w).^2/(2*sigma^2));
    aMap = conv2(g', g, padarray(frame, [w w], 'replicate'), 'valid');
    aMap = (aMap-2*pi*sigma^2*cMap) / gL2;
    
    % figure; imagesc(aMap); colormap(gray(256)); axis image; colorbar;
    testMap = aMap > sigmaT * cStdMap;
    % figure; imagesc(testMap); colormap(gray(256)); axis image; colorbar;
    % figure; imagesc(mask); colormap(gray(256)); axis image; colorbar;
    % figure; imagesc(ch2rgb(mask, testMap, dmask)); axis image;
    
    %falsePos = testMap-mask;
    %falsePos(falsePos<0) = 0;
    
    %p = sum(sum(dmask & testMap)) / sum(dmask(:))
    p(k) = sum(testMap(:)) / sum(dmask(:));
end
p = mean(p);

% for p = 1:npRef
%     window = frame(yi(p)-w:yi(p)+w, xi(p)-w:xi(p)+w);
%     % pick random mask
%     mi = unidrnd(length(xm));
%     maskWindow = ~mask(ymi(mi)-w:ymi(mi)+w, xmi(mi)-w:xmi(mi)+w);
%     
%     % background estimate
%     c = mean(window(maskWindow));
%     cStd = std(window(maskWindow));
% 
%     prm = fitGaussian2D(window, [0 0 max(window(:))-c sigma c], 'xyA');
%     %prm = fitGaussian2D(window, [x(p)-xi(p) y(p)-yi(p) max(window(:))-c sigma c], 'A');
%     %prm = fitGaussian2D(window, [xm(mi)-xmi(mi) ym(mi)-ymi(mi) max(window(:))-c sigma c], 'A');
%     A = prm(3);
%     ctrlStatus(p) = A > sigmaT*cStd;
%     
%     if A > sigmaT*cStd
%         prm
%         figure; imagesc(window); colormap(gray(256)); axis image;
%         hold on;
%         plot(xi(p) + prm(1), yi(p) + prm(2), 'rx');
%         figure; imagesc(maskWindow); colormap(gray(256)); axis image;
%         return
%     end
%     
%     
% end
% fprintf('Ch. %d background detections: %d/%d (%.2f %%)\n', 2, sum(ctrlStatus), npRef, 100*sum(ctrlStatus)/npRef);



% for c = slaveChannels
% 
%     % list of frames channel
%     
%     
% 
% bgMask = filterGauss2D(frame, 4, 'symmetric');
% bgMask = bgMask/max(bgMask(:));
% %bgMask = im2bw(bgMask, graythresh(bgMask(mask~=1)));
% 
% bgMask = im2bw(bgMask, graythresh(bgMask));
% sum(bgMask(:))
% 
% %bgMask = imclose(bgMask, strel('disk', 20));
% figure; imagesc(frame); colormap(gray(256)); axis image;
% figure; imagesc(bgMask); colormap(gray(256)); axis image;
% 
% % 
% 


%=================================
% Proportion of 'random' positives
%=================================

%=================================
% Classify tracks in slave channels 
%=================================

% steps for required to reject H_0: binomial
% nBinSteps = ceil(log(alpha)/log(0.5));

% Loops through all the tracks
for k = 1:length(tracks);
    cStatus = NaN(nChannels, 1);
    cStatusVector = NaN(nChannels, length(tracks(k).x)); 
    for c = 1:nChannels %slaveChannels
        binary = tracks(k).A(c,:) > sigmaT * tracks(k).cStd(c,:);
        % posLengths = find(diff([bd 0])==-1) - find(diff([0 bd 0])==1) + 1;
        
        %binoT = binoinv(0.95, length(tracks(k).x), sum(ctrlStatus)/npRef);
        binoT = binoinv(0.95, length(tracks(k).x), p);
        %if max(binary) > 0%nBinSteps
        if sum(binary) > binoT
            cStatus(c) = 1;
        else
            cStatus(c) = 0;
        end
        cStatusVector(c,:) = binary;
    end
    tracks(k).cStatus = cStatus;
    tracks(k).cStatusVector = cStatusVector;
end

for c = slaveChannels
    nPos = sum(arrayfun(@(x) x.cStatus(c), tracks));
    fprintf('Ch. %d positive tracks: %d/%d (%.2f %%)\n', c, nPos, length(tracks), 100*nPos/length(tracks));
end



% Cell mask: sum of first 10 frames
% frameInt = zeros(data.imagesize);
% for k = 1:10
%     frame = double(imread([data.source frameList(k).name]));
%     frameInt = frameInt + frame;
% end

% Sum all masks, generate cell mask
% maskInt = zeros(data.imagesize);
% for k = 1:data.movieLength
%     mask = double(imread([data.source 'Detection' filesep 'Masks' filesep maskList(k).name]));
%     maskInt = maskInt + mask;
% end
% maskInt(maskInt~=0) = 1;
% maskInt = imclose(maskInt, strel('disk', 2*w));
% maskInt = ~imclose(~maskInt, strel('disk', 2*w));
% 
% 
% 