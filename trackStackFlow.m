function [v,corLength,sigtVal] = trackStackFlow(stack,points,minCorL,varargin)
%trackStackFlow: Calculate the flow velocity from a stack of movie images.
%
% SYNOPSIS :
%    v = trackStackFlow(stack,points,minCorL)
%    [v,corLen] = trackStackFlow(stack,points,minCorL,maxCorL,varargin)
%    [v,corLen,sigtVal] = trackStackFlow(stack,points,minCorL,maxCorL,varargin)
%
% INPUT :
%    stack : An image stack (i.e. of dimensions n x m x l) to be correlated
% 
%    points : A set of points (size nP x 2) expressed in the xy 
%             coordinate system where the velocity is calculated.
%
%    minCorL : The minimum side length of an image block (or band)
%            that is to be cross correlated over frames to detect flow
%            velocity. Optimal block size will be searched in the range
%            [minCorL maxCorL] for the detection of coherent flow pattern
%            behind noisy data.
%
%    maxCorL : Optional - The maximum side length of an image block (or band)
%            that is to be cross correlated over frames to detect flow
%            velocity. Optimal block size will be searched in the range
%            [minCorL maxCorL] for the detection of coherent flow pattern
%            behind noisy data.
%            If not input, will be set to minCorL
%
%    The following optional parameters can be set as parameter/value pairs:
%
%    'bgAvgImg': A stack of stationary background frames to be substracted 
%                during image correlation. Default is zeros matrix.
%
%    'maxSpd'  : A numerical value that specifies the maximum speed that can
%                be detected (in pixels/frame). The default is 10.
%
%    'bgMask':   A stack of background masks which is used to remove
%                background pixels from being used in correlation.
%
%    'minFeatureSize': The minimum feature size in the image. This is 
%                      measured as the diameter of features.
%                      Default, 11 pixels (typical speckle size).
%
% OUTPUT :
%    v      : velocity vector of (size nP x2) expressed in the xy
%             coordinate system.
%
%    corLen : The optimal block length in the sense that it is the minimum
%             block length in the range [minCorLen, maxCorLe] that gives a
%             stable coherent flow.
%
%    sigtVal : The 1st, 2nd local maximum and the reference score for
%              significance test can also be output for use in
%              postprocessing.
%
% References:
% J. Li & G. Danuser, J. of microscopy, 220 150-167, 2005.

% Lin Ji, 2005
% Sebastien Besson, May 2011 (last modified Nov 2011)
% Adapted from imFlowTrack.m
% Sangyoon Han, October 2012 (last modified Mar 2013)

% Input check
ip= inputParser;
ip.addRequired('stack',@(x) isnumeric(x) && size(x,3)>=2);
ip.addRequired('points',@(x) isnumeric(x) && size(x,2)==2);
ip.addRequired('minCorL',@isscalar);
ip.addOptional('maxCorL',minCorL,@isscalar);
ip.addParamValue('maxSpd',10,@isscalar);
ip.addParamValue('bgMask',true(size(stack)),@(x) isequal(size(x),size(stack)));
ip.addParamValue('bgAvgImg', zeros(size(stack)),@isnumeric);
ip.addParamValue('minFeatureSize',11,@isscalar);
ip.parse(stack,points,minCorL,varargin{:});
maxCorL=ip.Results.maxCorL;
maxSpd=ip.Results.maxSpd;
minFeatureSize=ip.Results.minFeatureSize;
bgMask=ip.Results.bgMask;
bgAvgImg=ip.Results.bgAvgImg;

% SH: Poly-fit version

%We automatically update the speed search radius until a high limit is
% reached. If no significant maximum is detected, it means either the image
% quality is bad or the flow velocity is even higher than this limit.
maxSpdLimit = 2*maxSpd;

[imgW imgL numFrames] = size(stack);
x=points(:,1);
y=points(:,2);

%Initial maximum speed components in both direction.
initMaxFlowSpd = 10;
initMaxPerpSpd = 10;
closenessThreshold = 0.25;

%For isotropic correlation.
maxSpdLimit = max(maxSpdLimit,initMaxFlowSpd);

% Initialize output
nPoints = size(points,1);
v = zeros(nPoints,2);
corLength = minCorL*ones(nPoints,1);
sigtValues = NaN*ones(nPoints,3);

%We calculate a score for each sampling velocity. The score is an average of
% the normalized cross-correlation coefficient of an image block that moves
% over consecutive frames at the sampling velocity. We sample the velocity
% by sampling the components of the velocity in two orthogonal directions and
% in the range [-maxSpd maxSpd]. The two orthogonal directions are parallel to
% the two sides of the square block.
bandDir = [1 0];
perpDir = [-bandDir(2) bandDir(1)];

%We only use odd correlation lengths greater than 3 pixels
minCorL = max(3,minCorL+(1-mod(minCorL,2)));
maxCorL = max(minCorL,maxCorL+(1-mod(maxCorL,2)));

bAreaThreshold = min(0.95*minCorL^2,maxCorL^2*0.5);

%Options for optimization.
options = optimset('GradObj','on','Display','off');

% Creates the format string for the numerical indexes
L=length(num2str(nPoints));
strg=sprintf('%%.%dd',L);
backSpc =repmat('\b',1,L);

%Calculate the correlation coefficient for each sampling velocity at
% each point.
startTime = cputime;
fprintf(1,['   Start tracking (total: ' strg ' points): '],nPoints);
parfor k = 1:nPoints
    fprintf(1,[strg ' ...'],k);
    
    sigtVal = [NaN NaN NaN];
    
    %Always get back the initial max speed for the new point.
    maxFlowSpd = initMaxFlowSpd;
    maxPerpSpd = initMaxPerpSpd;
    
    %Alway start with 'minCorL' for each new point.
    corL = minCorL;
    
    pass = 0;
    while pass == 0 && corL <= maxCorL
        
        %Create kymograph around each point. 'bandDir' is used as the direction
        % of the kymographed line.
        xI = round(x(k));
        yI = round(y(k));
        if xI < 1 || xI > imgL || yI < 1 || yI > imgW
            %The point is outside the image. Mark it untrackable.
            pass = 0;
            corL = 2*maxCorL;
            continue;
        end
        %Always get back the initial max speed for new 'corL'.
        % We devide the max speed by 2 due the use of the while loop below.
        maxFlowSpd = initMaxFlowSpd/2;
        maxPerpSpd = initMaxPerpSpd/2;
        
        %Flag that indicates the quality of the score.
        pass = 0;
        while pass == 0 && maxFlowSpd < maxSpdLimit && maxPerpSpd < maxSpdLimit
            %If the quality of the score function is not good enough (pass == 0),
            % we increase the max sampling speed until the limit is reached.
            maxFlowSpd = min(maxSpdLimit,maxFlowSpd*2);
            maxPerpSpd = min(maxSpdLimit,maxPerpSpd*2);
            
            %Get sampling speed. Make sure it will not shift the template (block) outside of
            % the image area. We also use bigger stepsize for large speed.
            minSpdF = -min(floor(maxFlowSpd),max(0,xI-corL-1));
            maxSpdF = min(floor(maxFlowSpd),imgL-min(xI+corL,imgL));
%             dv = max(1,ceil(abs(minSpdF)/10));
%             vF = minSpdF:dv:min(0,maxSpdF);
%             if vF(end) ~= 0
%                 vF(end+1) = 0;
%             end
%             dv = max(1,ceil(abs(maxSpdF)/10));
%             vF = [vF vF(end)+dv:dv:maxSpdF];
%             
            minSpdP = -min(floor(maxPerpSpd),max(0,yI-corL-1));
            maxSpdP = min(floor(maxPerpSpd),imgW-min(yI+corL,imgW));
%             dv = max(1,ceil(abs(minSpdP)/10));
%             vP = minSpdP:dv:min(0,maxSpdP);
%             if vP(end) ~= 0
%                 vP(end+1) = 0;
%             end
%             dv = max(1,ceil(abs(maxSpdP)/10));
%             vP = [vP vP(end)+dv:dv:maxSpdP];
%             
            %Debugging
            vF = minSpdF:maxSpdF;
            vP = minSpdP:maxSpdP;
            %End of debugging
            
            hCLL    = min(xI-1,corL)+max(-vF(1),0);
            hCLR    = min(imgL-xI,corL)+max(vF(end),0);
            hCWL    = min(yI-1,corL)+max(-vP(1),0);
            hCWR    = min(imgW-yI,corL)+max(vP(end),0);
            cropL   = hCLL+hCLR+1;
            cropW   = hCWL+hCWR+1;
            kym     = zeros(cropW,cropL,numFrames);
            for k2 = 1:numFrames
                kym(:,:,k2) = double(stack(yI-hCWL:yI+hCWR,xI-hCLL:xI+hCLR,k2));
            end
            kymMask   = squeeze(bgMask(yI-hCWL:yI+hCWR,xI-hCLL:xI+hCLR,:));
            kymAvgImg = bgAvgImg(yI-hCWL:yI+hCWR,xI-hCLL:xI+hCLR,:);
            
            %The index of zero velocity.
            zeroI = [find(vP==0) find(vF==0)];
            
            %The index of the center of image block in the kymographed or cropped
            % image.
            centerI = [hCWL+1,hCLL+1];
            
            [score,blockIsTooSmall] = calScore(kym,centerI,corL,vP,vF, ...
                'bAreaThreshold',bAreaThreshold,'kymMask',kymMask,'kymAvgImg',kymAvgImg);
            
            if blockIsTooSmall
                %Tell the program to increase block size.
                pass = 0;
                maxFlowSpd = Inf;
                maxPerpSpd = Inf;
            else
                %Test the quality of the score function and find the index of the
                % maximum score.
                [pass,locMaxI,sigtVal] = findMaxScoreI(score,zeroI,minFeatureSize);
                if pass == 0 || corL < maxCorL
                    %Increase the block length and width by a factor of 5/4 to see if
                    % the ambiguity can be resovled. Also by comparing the two
                    % velocities returned from two block sizes, we identify the
                    % optimal block size that gives a coherent flow.
                    [score2] = calScore(kym,centerI,ceil(1.25*corL),...
                        vP,vF,'bAreaThreshold',bAreaThreshold,...
                        'kymMask',kymMask,'kymAvgImg',kymAvgImg);
                    if max(length(vF),length(vP))>160 %applying more generous threshold for higher velocity
                        [pass2,maxI2] = findMaxScoreI(score2,zeroI,minFeatureSize,0.8);
                    elseif max(length(vF),length(vP))>80 %applying more generous threshold for higher velocity
                        [pass2,maxI2] = findMaxScoreI(score2,zeroI,minFeatureSize,0.65);
                    elseif max(length(vF),length(vP))>40 %applying more generous threshold for higher velocity
                        [pass2,maxI2] = findMaxScoreI(score2,zeroI,minFeatureSize,0.58);
                    else
                        [pass2,maxI2] = findMaxScoreI(score2,zeroI,minFeatureSize,0.5);
                    end
                    if pass2 == 1
                        maxV2 = maxInterpfromScore(maxI2,score2,vP,vF);
                        %locMaxV = zeros(size(locMaxI,1),2);
                        locMaxV = [vP(locMaxI(:,1)).' vF(locMaxI(:,2)).'];

                        for j = 1:size(locMaxI,1)
                            maxIc = locMaxI(j,:);
                            maxV = maxInterpfromScore(maxIc,score,vP,vF);
                            locMaxV(j,:) = maxV;
                        end
                        
                        distToMaxV2 = sqrt(sum((locMaxV- ...
                            ones(size(locMaxV,1),1)*maxV2).^2,2));
                        [minD,ind] = min(distToMaxV2);
                        maxV = locMaxV(ind,:);
                        maxVNorm = max(norm(maxV2),norm(maxV));
                        if maxVNorm == 0 || ...
                                (pass == 1 && minD < 2*closenessThreshold*maxVNorm) || ...
                                (pass == 1 && maxVNorm < 0.5) || ...
                                (pass == 0 && minD < closenessThreshold*maxVNorm)
                            pass = 2;
                        else
                            pass = 0;
                        end
                    else
                        pass = 0;
                    end
                else
                    maxI = locMaxI;
                end
            end
        end
        
        if pass == 0
            if corL == maxCorL
                corL = Inf;
            else
                corL = min(maxCorL,floor(corL*3/2));
            end
        end
    end
    
    if pass == 0
        maxV = [NaN NaN];
        sigtVal = [NaN NaN NaN];
    elseif pass == 1
        maxV = maxInterpfromScore(maxI,score,vP,vF);
    end
    
    if ~isnan(maxV(1)) && ~isnan(maxV(2))
        rotv= maxV*[perpDir;bandDir];
        v(k,:) = [rotv(1) rotv(2)];
    else
        v(k,:) = [NaN NaN];
    end
    corLength(k) = corL;
    sigtValues(k,:) = sigtVal;
    
    fprintf(1,[backSpc '\b\b\b\b']);
end
nanInd = find(isnan(v(:,1)));
endTime = cputime;
fprintf(1,[strg '.\n'],nPoints);
fprintf(1,'   Tracking is done in %f sec (%f sec per point).\n', ...
    endTime-startTime,(endTime-startTime)/nPoints);
fprintf(1,'   Total tracked points: %d (out of %d).\n', ...
    nPoints-length(nanInd),nPoints);




function [score,blockIsTooSmall] = calScore(kym,centerI,corL,vP,vF,varargin)
% centerI : The coordinate index of the image block center in the 'kym' image.
% SH: this local function was updated for bead tracking instead of speckle
% tracking. The first slice of kym contains the reference frame and the
% second one has the current bead image.

%For boundary points, the cutoff of background can make the effective image area too small for
% stable correlation. We check this and report it back as output.
blockIsTooSmall = 0;

if mod(corL,2) == 0, corL = corL+1; end

numFrames = size(kym,3);
kymLen    = size(kym,2);
kymWidth  = size(kym,1);

%Check additional parameters
ip =inputParser;
ip.addParamValue('bAreaThreshold',0.5*corL^2,@isscalar);
ip.addParamValue('kymMask',[],@islogical)
ip.addParamValue('kymAvgImg',zeros(size(kym)),@isnumeric)
ip.parse(varargin{:});
bAreaThreshold=ip.Results.bAreaThreshold;
kymMask=ip.Results.kymMask;
kymAvgImg=ip.Results.kymAvgImg;

%The index of the correlating image block in the big cropped image.
bI1 = centerI(1)-(corL-1)/2:centerI(1)+(corL-1)/2;
bI2 = centerI(2)-(corL-1)/2:centerI(2)+(corL-1)/2;

%Find the part of the image block that is outside the cropped image and cut it off from the template.
bI1(bI1<1 | bI1>kymWidth) = [];
bI2(bI2<1 | bI2>kymLen) = [];

score = zeros(length(vP),length(vF));
if min(min(kymMask(:,:,1))) == 1
    %Background intensities are set to be zero. So, if there is no zero intensities, there is no
    % background.
    kymP2 = kym.*kym;
    
    %The norm of the kymographed image band at each frame.
    bNorm1 = sqrt(sum(sum(sum(kymP2(bI1,bI2,1:numFrames-1)))));
    
    for j1 = 1:length(vP)
        v1 = vP(j1);
        for j2 = 1:length(vF)
            v2 = vF(j2);
            corrM = kym(bI1,bI2,1:numFrames-1).* ...
                kym(bI1+v1,bI2+v2,2:numFrames);
            
            %The norm of the shifted image band at each frame.
            bNorm2 = sqrt(sum(sum(sum(kymP2(bI1+v1,bI2+v2,2:numFrames)))));
            
            %Normalize the correlation coefficients.
            score(j1,j2) = sum(corrM(:))/bNorm1/bNorm2;
        end
    end
else
    kym       = reshape(kym,kymLen*kymWidth,numFrames);
    kymMask   = reshape(kymMask,kymLen*kymWidth,numFrames);
    
    numAvgImgs    = size(kymAvgImg,3);
    kymAvgImg     = reshape(kymAvgImg,kymLen*kymWidth,numAvgImgs);
    lastKymAvgImg = kymAvgImg(:,numAvgImgs);
    
    %Extend 'kymAvgImg' to have 'numFrames' columns by repeating 'lastKymAvgImg'.
    kymAvgImg = [kymAvgImg lastKymAvgImg*ones(1,numFrames-numAvgImgs)];
    %kymAvgImg = kymAvgImg(:)*ones(1,numFrames);
    
    %%%%% Debugging %%%%%%%%%%%%%%%
    %kymAvgImg(:) = 0;
    %%%%% Debugging %%%%%%%%%%%%%%%
    
    kym0   = (kym-kymAvgImg).*kymMask;
    %kym0P2 = kym0.*kym0;
    
    %We only consider frames whose texture area (after cutting off background
    % area) is bigger than 'bAreaThreshold'. Note: Background pixel values are
    % zero.
    %First, Get the linear index of the image block in the big cropped image.
    [BI1,BI2] = ndgrid(bI1,bI2);
    bI = (BI2(:)-1)*kymWidth+BI1(:);  
    allFrames = 1:numFrames-1;
    nzInd = arrayfun(@(x)sum(kymMask(bI,x)),allFrames);    
    validFrames = allFrames(nzInd >= bAreaThreshold);
    
    %If the number of valid frames is less than half of the number of
    % correlating frames, we reject the tracking for this point with the default
    % zero score.
    if length(validFrames) < numFrames/2;
        blockIsTooSmall = 1;
        return;
    end
    
    %We consider template shift in both the positive (to next frame) flow direction and
    % negative (from previous frame) flow direction.
    kymValid   = kym0(bI,validFrames);
    kymValidP2 = kymValid.*kymValid;
    bNorm      = sqrt(sum(kymValidP2,1));
    
    % SB:old score calculation function from imFlowTrack
    %     for j1 = 1:length(vP)
    %       v1 = vP(j1);
    %       for j2 = 1:length(vF)
    %          v2 = vF(j2);
    %          v = v2*kymWidth+v1;
    %
    %          kymShift = kym(bI+v,validFrames+1)-kymAvgImg(bI,validFrames);
    %          bNormS   = sqrt(sum(kymShift.^2.*kymMask(bI,validFrames),1));
    %
    %          corrM   = -ones(1,length(validFrames));
    %          nzInd   = find(bNorm~=0);
    %          nzInd   = nzInd(bNormS(nzInd)~=0);
    %          zeroInd = find(bNorm==0);
    %          zeroInd = zeroInd(bNormS(zeroInd)==0);
    %
    %          corrM(zeroInd) = 1;
    %          corrM(nzInd)   = sum(kymValid(:,nzInd).*kymShift(:,nzInd),1);
    %
    %          zeroInd = find(bNorm==0);
    %          if ~isempty(zeroInd)
    %             bNorm(zeroInd) = 1;
    %          end
    %          zeroInd = find(bNormS==0);
    %          if ~isempty(zeroInd)
    %             bNormS(zeroInd) = 1;
    %          end
    %
    %          score(j1,j2) = mean(corrM./bNorm./bNormS);
    %       end
    %     end
    
    % SB: beginning of vectorized score calculation
    [v1,v2]=ndgrid(vP,vF);
    v = v2*kymWidth+v1;    
    [bI2,v2]=ndgrid(bI,v(:));
    allbI=bI2+v2;
    
    % Create a matrix of size (size(bI,1)xsize(validFrames)xsize(v))
    kymShiftMatrix= reshape(kym(allbI(:),validFrames+1),...
        length(bI),numel(v),length(validFrames))-...
        permute(repmat(kymAvgImg(bI,validFrames),[1 1 numel(v)]),[1 3 2]);
    kymShiftMatrix= permute(kymShiftMatrix,[2 3 1]);
    bNormS=squeeze(sqrt(sum(kymShiftMatrix.^2.*...
        permute(repmat(kymMask(bI,validFrames),[1 1 numel(v)]),[3 2 1]),3)));
    validCorrM = squeeze(sum(kymShiftMatrix.* ...
        permute(repmat(kym0(bI,validFrames),[1 1 numel(v)]),[3 2 1]),3));
    
    % Initialize the correlation matrix
    corrM = -ones(numel(v),length(validFrames));
    % Set the correlation value of zero correlation elements to 1
    corrM(repmat(bNorm==0,numel(v),1) & bNormS==0)=1;
    nzInd = repmat(bNorm~=0,numel(v),1) & bNormS~=0;
    corrM(nzInd) = validCorrM(nzInd);
    
    % Set bNorm and bNormS null components to 1 (to avoid division by zero)
    bNorm(bNorm==0)=1;
    bNorm=repmat(bNorm,numel(v),1);
    bNormS(bNormS==0)=1;
    score = mean(corrM./bNorm./bNormS,2);
    score = reshape(score,size(v));
    
    minusOnesI = find(score(:)==-1);
    nMOnesI    = (score(:)~=-1);
    [minScore,minScoreI] = min(score(nMOnesI));
    ind = 1:length(score(:));
    ind([minusOnesI minScoreI]) = [];
    score([minusOnesI minScoreI]) = min(score(ind));
end



function [pass,locMaxI,sigtVal] = findMaxScoreI(score,zeroI,minFeatureSize,sigThreshold)
%
% INPUT:
%    score : The cross-correlation score function.
%    zeroI : The index of 'score' that corresponds to zero velocity.
% OUTPUT:
%    pass  : If an unambiguous global maximum is found, pass = 1 is returned.
%            Otherwise, pass = 0 is returned indicating that the quality of the
%            score function is not good.
%    locMaxI : Index of local maximum whose scores pass the significant test.
%    sigtVal : A 1x3 vector that contains the scores of the 1st, 2nd local
%              maximum and the reference score.
if nargin < 4
    sigThreshold = 0.5; %0.5;
end

[m,n] = size(score);
numSamples = m*n;

avg = sum(score(:))/numSamples;
% dev = sum(abs(score(:)-avg))/numSamples;

%Some threshold used for quality test.
minFeatureRadius   = max(1,ceil((minFeatureSize-1)/2)); %Unit: pixel.
closenessThreshold = 0.2;
maxNbDist          = 0.5;

%We divide the score into 5x5 pixels blocks centered at 'zeroI' which is
% the index for zero velocity and find the local maximum in each block
% and then do a significant test on them in the sense that they have to be
% significantly bigger than the average. If there are more than one
% significant local maximums, return 0. To return 1, the unique significant
% local maximum also needs to be close the center of the sample region.

ind1 = 1:5:m;
if ind1(end) < m
    if ind1(end) == m-1
        ind1(end) = m;
    else
        ind1(end) = floor((ind1(end-1)+m)/2);
        ind1(end+1) = m;
    end
end
ind2 = 1:5:n;
if ind2(end) < n
    if ind2(end) == n-1
        ind2(end) = n;
    else
        ind2(end) = floor((ind2(end-1)+n)/2);
        ind2(end+1) = n;
    end
end

locMaxI    = [];
locMaxS    = [];
locAvgMinS = [];
locCount   = 0;
for k = 1:length(ind1)-1
    for j = 1:length(ind2)-1
        locScore = score(ind1(k):ind1(k+1),ind2(j):ind2(j+1));
        [tmp,index] = max(locScore,[],1);
        [maxS,i2]   = max(tmp);
        i1 = index(i2);
        %maxI = [ind1(k,i1) ind2(j,i2)];
        maxI = [ind1(k)+i1-1 ind2(j)+i2-1];
        
        %Check if it is true local maximum in the sense that it is bigger than
        % its surrounding pixels or it is a boundary maximum.
        if maxI(1) == 1
            yOffset = maxI(1):maxI(1)+2;
        elseif maxI(1) == m
            yOffset = maxI(1)-2:maxI(1);
        else
            yOffset = maxI(1)-1:maxI(1)+1;
        end
        if maxI(2) == 1
            xOffset = maxI(2):maxI(2)+2;
        elseif maxI(2) == n
            xOffset = maxI(2)-2:maxI(2);
        else
            xOffset = maxI(2)-1:maxI(2)+1;
        end
        if maxS >= max(max(score(yOffset,xOffset)))
            %maxS = sum(sum(score(yOffset,xOffset)))/9; %This caused flow
            %underestimation. We should use a single maximum value at the
            %maximum velocity position rather than averaging with neiboring
            %points. This can prevent a value at the border of the score
            %from not being captured as a miximum, which will lead to
            %expansion of correlation length. BTW, what was the reason of
            %averaging maximum score with neighboring scores? To prevent
            %very narrow peak from being true maximum? I don't think
            %that'll happen. - Sangyoon 3/2/2013
            
            %The following 'avgMinS' is used in the calibration of the
            % 'baseS' below. It is the averge of scores around the local
            % maximum (excluding the maximum).
            avgMinS = (sum(sum(score(yOffset,xOffset)))-score(maxI(1),maxI(2)))/8;
            
            %Further check if it is close to any previouse
            % local maximum.
            %Distance to previouse local maximum.
            dist = zeros(locCount,1);
            for jj = 1:locCount
                dist(jj) = norm(maxI-locMaxI(jj,:));
            end
            [minD,minDI] = min(dist);
            if locCount == 0 || minD >= ...
                    max(2,norm(locMaxI(minDI,:)-zeroI)*closenessThreshold)
                locCount = locCount+1;
                
                locMaxI(locCount,:)  = maxI;
                locMaxS(locCount)    = maxS;
                locAvgMinS(locCount) = avgMinS;
            elseif maxS > locMaxS(minDI)
                locMaxI(minDI,:)  = maxI;
                locMaxS(minDI)    = maxS;
                locAvgMinS(minDI) = avgMinS;
            end
        end
    end
end

[maxS,maxInd] = max(locMaxS);
maxI = locMaxI(maxInd,:);
[locMaxS,desI] = sort(locMaxS,'descend');
locMaxI    = locMaxI(desI,:);
maxI       = locMaxI(1,:);
locAvgMinS = locAvgMinS(desI);
avgMinS    = locAvgMinS(1); %Note: this is not global 'minS'. It is the
% minimum score among the 9 elements around
% 'maxI'. It is used in the calibartion of
% 'baseS' below.
maxINorm = norm(maxI-zeroI);

if size(locMaxI,1) == 1
    sigtVal = [maxS 0 maxS];
    if maxI(1) < m/4 || maxI(1) > 3*m/4 || ...
            maxI(2) < n/4 || maxI(2) > 3*n/4
        pass = 0;
    else
        pass = 1;
    end
    return;
end

%Calculate the distance from all local maximum to the global maximum. It
%will be used to determine the neiborhood for calculating local average
%around the global maximum point.
dist = sqrt((locMaxI(2:end,1)-maxI(1)).^2+(locMaxI(2:end,2)-maxI(2)).^2);
minD = min(dist);

%Calculate the local averge around the maximun point. This local average
% will be used to calculate 'baseS', the reference score.
offset = max(minFeatureRadius,ceil(min(minD,maxINorm)*maxNbDist));
if maxI(1) > offset
    y0 = maxI(1)-offset;
else
    y0 = 1;
end
yNeighbor = y0:min(m,y0+2*offset);

if maxI(2) > offset
    x0 = maxI(2)-offset;
else
    x0 = 1;
end
xNeighbor = x0:min(n,x0+2*offset);
locAvg = sum(sum(score(yNeighbor,xNeighbor)))/ ...
    length(yNeighbor)/length(xNeighbor);

baseS  = locMaxS;
for k = 2:length(locMaxS)
    %Anisotropically adapted reference score:
    % Calculate local minimum along the line between the two local maximums
    % for the testing of local maximum significance.
    %The distance between the two local maximum.
    xD = locMaxI(k,2)-maxI(2);
    yD = locMaxI(k,1)-maxI(1);
    if abs(xD) > abs(yD)
        xShift = 0:sign(xD):xD;
        yShift = floor(xShift*yD/xD+0.5);
    else
        yShift = 0:sign(yD):yD;
        xShift = floor(yShift*xD/yD+0.5);
    end
    lineMinS = min(score(m*(maxI(2)+xShift-1)+maxI(1)+yShift));
    baseS(k)   = max(lineMinS,locAvg);
    
    %In very rare cases, 'baseS' can be even bigger than 'maxS' since our
    % 'maxS' is the average of the 9 elements around 'maxI'. So, we need
    % the following control to make sure 'baseS' is always less than
    % 'maxS'.
    if baseS(k) > maxS
        baseS(k) = min(avgMinS,baseS(k));
    end
end
if length(baseS) == 1
    sigtVal = [maxS 0 maxS];
elseif baseS(2) == maxS
    sigtVal = [maxS locMaxS(2) maxS/2];
else
    sigtVal = [maxS locMaxS(2) baseS(2)];
end
inSigI = find(maxS-baseS > 0 & locMaxS-baseS < (maxS-baseS)*sigThreshold);
locMaxI(inSigI,:) = [];
locMaxS(inSigI)   = [];

if isempty(locMaxS)
    pass = 0;
    return;
elseif length(locMaxS) > 1
    pass = 0;
    return;
end

if maxI(1) < m/4 || maxI(1) > 3*m/4 || ...
        maxI(2) < n/4 || maxI(2) > 3*n/4
    pass = 0;
    return;
end

pass = 1;

function maxV2 = maxInterpfromScore(maxI2,score,vP,vF)
% Sangyoon: I made a change for this refining process to
% use parabola approximation. Once parabola fit is
% too much apart from integer maxV (maxVorg), I
% started to use the fmincon again for more correct refining process.
% parabola approximation
% input:    maxI2       :index for maxV in score
%           score       :score
%           vP,vF       :velocity range
% output:   maxV2       :refined velocity

subv = 1; % radius of subgroup for subscore
maxVorg  = [vP(maxI2(1)) vF(maxI2(2))];

bPolyTracked = 0;
if (maxI2(1)-subv)>=1 && (maxI2(1)+subv)<=size(score,1)...
   && (maxI2(2)-subv)>=1 && (maxI2(2)+subv)<=size(score,2)
    subv = 1; % radius of subgroup for subscore
    sub_score = score(max(1,maxI2(1)-subv):min(size(score,1),maxI2(1)+subv),...
                        max(1,maxI2(2)-subv):min(size(score,2),maxI2(2)+subv));
    % my field of interest
    subvP = vP(max(1,maxI2(1)-subv):min(size(score,1),maxI2(1)+subv));
    subvF = vF(max(1,maxI2(2)-subv):min(size(score,2),maxI2(2)+subv));

    [subvFG,subvPG]=meshgrid(subvF,subvP);
    subvF1D = reshape(subvFG,[],1);
    subvP1D = reshape(subvPG,[],1);
    sub_score1D = reshape(sub_score,[],1);

    % starting point estimation SH based on discretized maxV (-b/2a =
    % maxVorg(2)) in quadratical expression to avoid the random starting point warning SH
    asp = -0.026; %decided empirically
    bsp = -2*asp*maxVorg(2);
    csp = asp;
    dsp = -2*csp*maxVorg(1);
    esp = -0.5; %arbitrary number
    s = fitoptions('Method','NonlinearLeastSquares','StartPoint', [asp,bsp,csp,dsp,esp]); 
    f = fittype('a*x^2+b*x+c*y^2+d*y+e','independent', {'x', 'y'}, 'dependent', 'z','option',s);
    sf = fit( [subvF1D, subvP1D], sub_score1D, f);

    px = [sf.a sf.b sf.e]; py = [sf.c sf.d sf.e];
    maxV2 = [roots(polyder(py)) roots(polyder(px)) ];
    bPolyTracked = 1;
end

if ~bPolyTracked||norm(maxVorg-maxV2,2)>1 %checking for proximity of the fitted point to the original discrete point
    subv = 4; % expanding region to fit
    sub_score = score(max(1,maxI2(1)-subv):min(size(score,1),maxI2(1)+subv),...
                        max(1,maxI2(2)-subv):min(size(score,2),maxI2(2)+subv));
    % my field of interest
    subvP = vP(max(1,maxI2(1)-subv):min(size(score,1),maxI2(1)+subv));
    subvF = vF(max(1,maxI2(2)-subv):min(size(score,2),maxI2(2)+subv));

    sp   = csape({subvP,subvF},sub_score);
    dsp1 = fnder(sp,[1,0]);
    dsp2 = fnder(sp,[0,1]);
    options = optimset('Algorithm','interior-point'); % this generates an warning
%             options = optimset('Algorithm','sqp');% this too

    maxV2 = fmincon(@vFun,maxVorg,[],[],[],[], ...
        [max(subvP(1),maxVorg(1)-2) max(subvF(1),maxVorg(2)-2)], ...
        [min(subvP(end),maxVorg(1)+2), min(subvF(end),maxVorg(2)+2)],[], ...
        options,sp,dsp1,dsp2);            
end

