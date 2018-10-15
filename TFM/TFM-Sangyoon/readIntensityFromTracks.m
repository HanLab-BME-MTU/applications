function tracksNA = readIntensityFromTracks(tracksNA,imgStack, attribute, varargin)
% tracksNA = readIntensityFromTracks(tracksNA,imgStack, attribute) records
% pixel intensity from imgStack for entire track, even after ANA state in the last x,y
% position, and store it to attribute in tracksNA (1 means ampTotal, 2
% means forceMag, 3 for fret, 4 for flowSpeed, 5 for ampTotal2 and 6 for
% ampTotal3) assumes that tracks reflects coordinates from
% stage-drift-corrected images. imgStack is also from SDC output.
% Set 'reTrack' to be true if you didn't retrack using paxillin image
% stack and want to read values from other stacks.
% Sangyoon Han, Jan 2015
% Last modified, Feb 2016
ip =inputParser;
ip.addParamValue('extraLength',300,@isscalar); % selcted track ids
ip.addParamValue('reTrack',true,@islogical); % selcted track ids
ip.addParamValue('trackOnlyDetected',false,@islogical); % selcted track ids
ip.addParamValue('movieData',[],@(x) isa(x,'MovieData') || isempty(x)); % moviedata for utrack
ip.parse(varargin{:});
extraLengthForced=ip.Results.extraLength;
reTrack=ip.Results.reTrack;
trackOnlyDetected =ip.Results.trackOnlyDetected;
MD =ip.Results.movieData;
extraLength = ip.Results.extraLength;
% get stack size
numFrames = size(imgStack,3);
% w4 = 8;
sigma = median(tracksNA(1).sigma);
numTracks = numel(tracksNA);
if isempty(MD)
    searchRadius = 1;
    maxR = 2;
    maxGap = 3;
    brownScaling = 1.01;
else
    iTrackingProc =MD.getProcessIndex('TrackingProcess');
    trackingProc = MD.getProcess(iTrackingProc);
    trackingParams = trackingProc.funParams_;
    minR=trackingParams.costMatrices(2).parameters.minSearchRadius;
    maxR=trackingParams.costMatrices(2).parameters.maxSearchRadius;
    searchRadius = (minR+maxR)/2;
%     searchRadiusDetected = maxR;
    brownScaling = trackingParams.costMatrices(2).parameters.brownScaling(1)+1;
    
    maxGap = trackingParams.gapCloseParam.timeWindow;
end

if attribute==2
    halfWidth=4;
    halfHeight=4;
else
    halfWidth=1;
    halfHeight=1;
end
%% %%%%%%%%%%%%%%%%%%%%%%5
% poolobj = gcp('nocreate'); % If no pool, do not create new one.
% myCluster = parcluster('local');
% if isempty(poolobj)
%     parpool(myCluster.NumWorkers);
% end
    

%% Field creation before running parfor
if ~isfield(tracksNA,'startingFrameExtra')
    tracksNA(end).startingFrameExtra=[];
end
if ~isfield(tracksNA,'startingFrameExtraExtra')
    tracksNA(end).startingFrameExtraExtra=[];
end
if ~isfield(tracksNA,'endingFrameExtra')
    tracksNA(end).endingFrameExtra=numFrames;
end
if ~isfield(tracksNA,'endingFrameExtraExtra')
    tracksNA(end).endingFrameExtraExtra=[];
end
if ~isfield(tracksNA,'ampTotal') && attribute==1
    tracksNA(end).ampTotal=tracksNA(end).amp;
end
if attribute==2 && ~isfield(tracksNA,'forceMag')
    tracksNA(end).forceMag=[];
elseif attribute==3 && ~isfield(tracksNA,'fret')
    tracksNA(end).fret=[];
elseif attribute==4 && ~isfield(tracksNA,'flowSpeed')
    tracksNA(end).flowSpeed=[];
elseif attribute==5 && ~isfield(tracksNA,'ampTotal2')
    tracksNA(end).ampTotal2=[];
elseif attribute==6 && ~isfield(tracksNA,'ampTotal3')
    tracksNA(end).ampTotal3=[];
end    
%% Change the old format
if attribute>1 && isfield(tracksNA,'state')
    tracksNA = changeTrackStateFormat( tracksNA );
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parfor_progress(numTracks);
progressText(0,'Re-reading and tracking individual tracks');
% progressbar
tic
if numTracks<100
  parforArg = 0;
else
  parforArg = Inf;
end

parfor (k=1:numTracks, parforArg)
% for k=1:numTracks
%     startFrame = max(1, min(arrayfun(@(x) x.startingFrame,tracksNA))-extraLength);
%     endFrame = min(numFrames, max(arrayfun(@(x) x.endingFrame,tracksNA))+extraLength);
%     startFrame = max(1, tracksNA(k).startingFrame-extraLength);
%     endFrame = min(numFrames,tracksNA(k).endingFrame+extraLength);
%     curStartingFrame = tracksNA(k).startingFrame;
%     curEndingFrame = tracksNA(k).endingFrame;
%     tracksNA(k).startingFrameExtra = startFrame;
%     tracksNA(k).endingFrameExtra = endFrame;
%     curRange = curStartingFrame:curEndingFrame;
%     if extraLength==0
%         if attribute==1
%             tracksNA(k).ampTotal(curRange) = arrayfun(@(x) imgStack(round(tracksNA(k).yCoord(x)),round(tracksNA(k).xCoord(x)),x),curRange);
%         elseif attribute==2
%             tracksNA(k).forceMag(curRange) = arrayfun(@(x) imgStack(round(tracksNA(k).yCoord(x)),round(tracksNA(k).xCoord(x)),x),curRange);
%         end
%     else
    % initialize amptotal to have it have the same dimension as .amp
    curTrack=tracksNA(k);
    if attribute==1
        if isempty(MD)
            searchRadiusDetected = 2;
        else
            searchRadiusDetected = maxR;
        end
        curTrack.ampTotal = curTrack.amp;
        try
            curStartingFrame = curTrack.startingFrameExtra;
            curEndingFrame = curTrack.endingFrameExtra;
            if isempty(curStartingFrame)
                curStartingFrame = curTrack.startingFrame;
            end
            if isempty(curEndingFrame)
                curEndingFrame = curTrack.endingFrame;
            end
        catch
            curStartingFrame = curTrack.startingFrame;
            curEndingFrame = curTrack.endingFrame;
        end
        if ~trackOnlyDetected
            % for the earlier time-points - going backward
            startFrame = max(1, curTrack.startingFrame-extraLength);
            endFrame = min(numFrames,curTrack.endingFrame+extraLength);
            ii=curStartingFrame;
            x = curTrack.xCoord(ii);
            y = curTrack.yCoord(ii);
            A = curTrack.amp(curStartingFrame);
            c = curTrack.bkgAmp(curStartingFrame); 
            curTrack.startingFrameExtra = curStartingFrame;
            curTrack.endingFrameExtra = curEndingFrame;

            gapClosed=0;
            for ii=curStartingFrame:-1:startFrame
                curImg = imgStack(:,:,ii); %#ok<*PFBNS>
                p=-1;
    %             curSigma = sigma;
    %             pitFound = false;
                while p<=15 %30
                    %oldP = p;
                    p=p+1;
                    pmP = -p; %ceil(p/2)*(-1)^oldP; % I removed the 'decreasing mode' because it also captures too much noise.
                    curSigma = sigma*(20-pmP)/20; % increasing sigma by 5 percent per each iteration
                    pitFound = false;
                    curAlpha = 0.05+p/100;
                    pstruct = fitGaussians2D(curImg, x, y, A, curSigma, c, 'xyac','Alpha',curAlpha);
                    if ~isnan(pstruct.x) && abs(pstruct.x-x)<searchRadius && abs(pstruct.y-y)<searchRadius && pstruct.A>0 && pstruct.A<2*A
                        x = pstruct.x;
                        y = pstruct.y;
                        A = pstruct.A;
                        c = pstruct.c; 
                        xi = round(x);
                        yi = round(y);
                        xRange = max(1,xi-halfWidth):min(xi+halfWidth,size(imgStack,2));
                        yRange = max(1,yi-halfHeight):min(yi+halfHeight,size(imgStack,1));
                        curAmpTotal = curImg(yRange,xRange);
                        curAmpTotal = mean(curAmpTotal(:));
                        curTrack.startingFrameExtra = ii;
                        curTrack.xCoord(ii) = x;
                        curTrack.yCoord(ii) = y;
                        curTrack.amp(ii) = A;
                        curTrack.bkgAmp(ii) = c;
                        curTrack.ampTotal(ii) =  curAmpTotal;
                        curTrack.presence(ii) =  true;
                        curTrack.sigma(ii) = curSigma;
                        if curTrack.state(ii)==1 || curTrack.state(ii)==5 %'BA','ANA')
                            curTrack.state(ii) = 2; %'NA';
                        end
                        pitFound = true;
                        if gapClosed>0
                            % we have to fill in the gap 
                            for kk=1:gapClosed
                                curTrack.xCoord(ii+kk) = ((gapClosed+1-kk)*curTrack.xCoord(ii)+kk*curTrack.xCoord(ii+gapClosed+1))/(gapClosed+1);
                                curTrack.yCoord(ii+kk) = ((gapClosed+1-kk)*curTrack.yCoord(ii)+kk*curTrack.yCoord(ii+gapClosed+1))/(gapClosed+1);
                                curTrack.amp(ii+kk) = ((gapClosed+1-kk)*curTrack.amp(ii)+kk*curTrack.amp(ii+gapClosed+1))/(gapClosed+1);
                                curTrack.bkgAmp(ii+kk) = ((gapClosed+1-kk)*curTrack.bkgAmp(ii)+kk*curTrack.bkgAmp(ii+gapClosed+1))/(gapClosed+1);

                                x = curTrack.xCoord(ii+kk);
                                y = curTrack.yCoord(ii+kk);
                                xi = round(x);
                                yi = round(y);
                                xRange = max(1,xi-halfWidth):min(xi+halfWidth,size(imgStack,2));
                                yRange = max(1,yi-halfHeight):min(yi+halfHeight,size(imgStack,1));
                                curAmpTotal = curImg(yRange,xRange);
                                curAmpTotal = mean(curAmpTotal(:));
                                curTrack.ampTotal(ii+kk) =  curAmpTotal;
                            end
                        end
                        gapClosed = 0;
                        break
                    end
                end
                if ~pitFound && gapClosed >= maxGap
                    curTrack.startingFrameExtra = ii+1+gapClosed;
                    break % We still have to read the values
                elseif ~pitFound && gapClosed < maxGap
                    gapClosed = gapClosed+1;
                end
%                 if ~pitFound % We read values regardless
%                     xi = round(x);
%                     yi = round(y);
%                     xRange = max(1,xi-halfWidth):min(xi+halfWidth,size(imgStack,2));
%                     yRange = max(1,yi-halfHeight):min(yi+halfHeight,size(imgStack,1));
%                     curAmpTotal = curImg(yRange,xRange);
%                     curAmpTotal = mean(curAmpTotal(:));
%                     curTrack.xCoord(ii) = x;
%                     curTrack.yCoord(ii) = y;
%                     curTrack.ampTotal(ii) =  curAmpTotal;
%                     curTrack.presence(ii) =  false;
%                 end                    
            end
        end
        %% for the present period - it is necessary for ampTotal
        ii=curStartingFrame+1;
        x = curTrack.xCoord(ii);
        y = curTrack.yCoord(ii);
        A = curTrack.amp(ii);
%                 c = curTrack.bkgAmp(ii); 
        trackingFromStartingFrame = true;
        
        gapClosed=0;
        for ii=curStartingFrame+1:curEndingFrame
            curImg = imgStack(:,:,ii);
            if ~reTrack
                xi = round(x);
                yi = round(y);
                xRange = max(1,xi-halfWidth):min(xi+halfWidth,size(imgStack,2));
                yRange = max(1,yi-halfHeight):min(yi+halfHeight,size(imgStack,1));
                curAmpTotal = curImg(yRange,xRange);
                curAmpTotal = mean(curAmpTotal(:));
                curTrack.ampTotal(ii) =  curAmpTotal;
            else
                p=-1; %This seems waisting the time. Now I am skipping...
                while p<=16
                    p=p+1;
                    pitFound = false;
                    pmP = -p; %ceil(p/2)*(-1)^oldP; % I removed the 'decreasing mode' because it also captures too much noise.
                    curSigma = sigma*(20-pmP)/20; % increasing sigma by 5 percent per each iteration
%                     curSigma = sigma*(20-p)/20; % increasing sigma by 5 percent per each iteration
%                     if ismember(curTrack.state(ii),[3,4]) % if this is FC or FC, we use anisotropic fitting
% %                         pstruct = fitAnisoGaussian2D(curImg, x, y, A, curSigma, c, 'xyAc');
%                         filterSigma = sigma*sqrt(2);
%                         [~,T] = steerableDetector(double(curImg),2,filterSigma);
%                         P = zeros(1, 7);
%                         P(:,1) = x;
%                         P(:,2) = y;
%                         P(:,3) = A;
%                         P(:,4) = 2*curSigma;       % sigmaX
%                         P(:,5) = curSigma;         % sigmaY
%                         P(:,6) = T(round(y),round(x))+pi/2;    % angle
%                         [params] = fitAnisoGaussian2D(double(curImg), P, 'xyArsct');
%                         
%                         pstruct.x=params(1);
%                         pstruct.y=params(2);
%                         pstruct.A=params(3);
%                         pstruct.c=params(4); 
%                     else
                        pstruct = fitGaussians2D(curImg, x, y, A, curSigma, c, 'xyAsc');
%                     end
%                     pstruct2 = fitGaussianMixtures2D(double(curImg), x, y, A, curSigma, c);
                    
                    if ~isnan(pstruct.x) && abs(pstruct.x-x)<searchRadiusDetected*2 && abs(pstruct.y-y)<searchRadiusDetected*2 && pstruct.A>0 
                        if trackingFromStartingFrame
                            trackingFromStartingFrame=false;
                            curTrack.startingFrameExtra = ii;
                        end
                        x = pstruct.x;
                        y = pstruct.y;
                        A = pstruct.A;
                        c = pstruct.c; 
                        xi = round(x);
                        yi = round(y);
                        xRange = max(1,xi-halfWidth):min(xi+halfWidth,size(imgStack,2));
                        yRange = max(1,yi-halfHeight):min(yi+halfHeight,size(imgStack,1));
                        curAmpTotal = curImg(yRange,xRange);
                        curAmpTotal = mean(curAmpTotal(:));
                        curTrack.xCoord(ii) = x;
                        curTrack.yCoord(ii) = y;
                        curTrack.amp(ii) = A;
                        curTrack.bkgAmp(ii) = c;
                        curTrack.ampTotal(ii) =  curAmpTotal;
                        curTrack.presence(ii) =  1;
                        curTrack.sigma(ii) = curSigma;
                        if curTrack.state(ii)==1 || curTrack.state(ii)==5 %'BA','ANA')
                            curTrack.state(ii) = 2; %'NA';
                        end
                        pitFound = true;
                        if gapClosed>0
                            % we have to fill in the gap 
                            for kk=1:gapClosed
                                curTrack.xCoord(ii-kk) = ((gapClosed+1-kk)*curTrack.xCoord(ii)+kk*curTrack.xCoord(ii-gapClosed-1))/(gapClosed+1);
                                curTrack.yCoord(ii-kk) = ((gapClosed+1-kk)*curTrack.yCoord(ii)+kk*curTrack.yCoord(ii-gapClosed-1))/(gapClosed+1);
                                curTrack.amp(ii-kk) = ((gapClosed+1-kk)*curTrack.amp(ii)+kk*curTrack.amp(ii-gapClosed-1))/(gapClosed+1);
                                curTrack.bkgAmp(ii-kk) = ((gapClosed+1-kk)*curTrack.bkgAmp(ii)+kk*curTrack.bkgAmp(ii-gapClosed-1))/(gapClosed+1);

                                x = curTrack.xCoord(ii-kk);
                                y = curTrack.yCoord(ii-kk);
                                xi = round(x);
                                yi = round(y);
                                xRange = max(1,xi-halfWidth):min(xi+halfWidth,size(imgStack,2));
                                yRange = max(1,yi-halfHeight):min(yi+halfHeight,size(imgStack,1));
                                curAmpTotal = curImg(yRange,xRange);
                                curAmpTotal = mean(curAmpTotal(:));
                                curTrack.ampTotal(ii-kk) =  curAmpTotal;
                            end
                        end
                        gapClosed = 0;
                        if isempty(MD)
                            searchRadiusDetected = 2;
                        else
                            searchRadiusDetected = maxR;
                        end
                        break
                    end
                end
                if ~pitFound && gapClosed >= maxGap
                    curTrack.endingFrameExtra = ii-gapClosed;
                    curEndingFrame = curTrack.endingFrameExtra;
                    if trackingFromStartingFrame % this case, we need to change the startingFrameExtra
                        curTrack.startingFrameExtra = ii;
                        gapClosed=0;
                    else
                        break
                    end
                elseif ~pitFound && gapClosed < maxGap
                    gapClosed = gapClosed+1;
                    searchRadiusDetected = searchRadiusDetected^brownScaling;
                    % If the search is over in this regime, read whatever
                    % in the same x-y position
                    
                end
                % It is a rare case, but it is possible that at some point
                % there is no significant point source detected during this
                % re-tracking. In this case, we change the curEndingFrame
                % to be the previous time point
                if ii==curEndingFrame %&& isnan(pstruct.x) && 
                    curEndingFrame=ii-1;
%                     curTrack.endingFrame = curEndingFrame;
                    curTrack.endingFrameExtra = curEndingFrame;
                    break
                end
                if ~pitFound % We read values regardless
                    xi = round(x);
                    yi = round(y);
                    xRange = max(1,xi-halfWidth):min(xi+halfWidth,size(imgStack,2));
                    yRange = max(1,yi-halfHeight):min(yi+halfHeight,size(imgStack,1));
                    curAmpTotal = curImg(yRange,xRange);
                    curAmpTotal = mean(curAmpTotal(:));
                    curTrack.xCoord(ii) = x;
                    curTrack.yCoord(ii) = y;
                    curTrack.ampTotal(ii) =  curAmpTotal;
                    curTrack.presence(ii) =  false;
                end
            end
        end
        % for the later time-points - going forward, x and y are already
        % set as a last point.
        if ~trackOnlyDetected
            x = curTrack.xCoord(curEndingFrame);
            y = curTrack.yCoord(curEndingFrame);
            A = curTrack.amp(curEndingFrame);
            c = curTrack.bkgAmp(curEndingFrame);
            gapClosed=0;

            for ii=(curEndingFrame+1):endFrame
                curImg = imgStack(:,:,ii);
                pitFoundEnd = false;
                p=-1;
                while p<=30
%                         oldP = p;
                    p=p+1;
                    pmP = -p; %ceil(p/2)*(-1)^oldP; % I removed the 'decreasing mode' because it also captures too much noise.
                    curSigma = sigma*(20-pmP)/20; % increasing sigma by 5 percent per each iteration
                    curAlpha = 0.05+p/100;
                    pstruct = fitGaussians2D(curImg, x, y, A, curSigma, c, 'xyac','Alpha',curAlpha);
                    if ~isnan(pstruct.x) && abs(pstruct.x-x)<searchRadiusDetected && abs(pstruct.y-y)<searchRadiusDetected && pstruct.A>0 && pstruct.A<2*A
                        x = pstruct.x;
                        y = pstruct.y;
                        A = pstruct.A;
                        c = pstruct.c; 
                        xi = round(x);
                        yi = round(y);
                        xRange = max(1,xi-halfWidth):min(xi+halfWidth,size(imgStack,2));
                        yRange = max(1,yi-halfHeight):min(yi+halfHeight,size(imgStack,1));
                        curAmpTotal = curImg(yRange,xRange);
                        curAmpTotal = mean(curAmpTotal(:));
                        curTrack.endingFrameExtra = ii;
                        curTrack.xCoord(ii) = x;
                        curTrack.yCoord(ii) = y;
                        curTrack.amp(ii) = A;
                        curTrack.bkgAmp(ii) = c;
                        curTrack.ampTotal(ii) =  curAmpTotal;
                        curTrack.presence(ii) =  true;
                        curTrack.sigma(ii) = curSigma;
                        if curTrack.state(ii)==1 || curTrack.state(ii)==5 %'BA','ANA')
                            curTrack.state(ii) = 2; %'NA';
                        end
%                             if strcmp(curTrack.state{ii},'BA') || strcmp(curTrack.state{ii},'ANA')
%                                 curTrack.state{ii} = 'NA';
%                             end
                        pitFoundEnd = true;
                        if gapClosed>0
                            % we have to fill in the gap 
                            for kk=1:gapClosed
                                curTrack.xCoord(ii-kk) = ((gapClosed+1-kk)*curTrack.xCoord(ii)+kk*curTrack.xCoord(ii-gapClosed-1))/(gapClosed+1);
                                curTrack.yCoord(ii-kk) = ((gapClosed+1-kk)*curTrack.yCoord(ii)+kk*curTrack.yCoord(ii-gapClosed-1))/(gapClosed+1);
                                curTrack.amp(ii-kk) = ((gapClosed+1-kk)*curTrack.amp(ii)+kk*curTrack.amp(ii-gapClosed-1))/(gapClosed+1);
                                curTrack.bkgAmp(ii-kk) = ((gapClosed+1-kk)*curTrack.bkgAmp(ii)+kk*curTrack.bkgAmp(ii-gapClosed-1))/(gapClosed+1);

                                x = curTrack.xCoord(ii-kk);
                                y = curTrack.yCoord(ii-kk);
                                xi = round(x);
                                yi = round(y);
                                xRange = max(1,xi-halfWidth):min(xi+halfWidth,size(imgStack,2));
                                yRange = max(1,yi-halfHeight):min(yi+halfHeight,size(imgStack,1));
                                curAmpTotal = curImg(yRange,xRange);
                                curAmpTotal = mean(curAmpTotal(:));
                                curTrack.ampTotal(ii-kk) =  curAmpTotal;
                            end
                        end
                        gapClosed = 0;
                        if isempty(MD)
                            searchRadiusDetected = 2;
                        else
                            searchRadiusDetected = maxR;
                        end
                        break
                    end
                end
                if ~pitFoundEnd && gapClosed >= maxGap
                    curTrack.endingFrameExtra = ii-gapClosed-1;
                    break
                elseif ~pitFoundEnd && gapClosed < maxGap
                    gapClosed = gapClosed+1;
                    searchRadiusDetected = searchRadiusDetected^brownScaling;
                    % If the search is over in this regime, read whatever
                    % in the same x-y position
                    
                end
                
%                 if ~pitFound % We read values regardless
%                     xi = round(x);
%                     yi = round(y);
%                     xRange = max(1,xi-halfWidth):min(xi+halfWidth,size(imgStack,2));
%                     yRange = max(1,yi-halfHeight):min(yi+halfHeight,size(imgStack,1));
%                     curAmpTotal = curImg(yRange,xRange);
%                     curAmpTotal = mean(curAmpTotal(:));
%                     curTrack.xCoord(ii) = x;
%                     curTrack.yCoord(ii) = y;
%                     curTrack.ampTotal(ii) =  curAmpTotal;
%                     curTrack.presence(ii) =  false;
%                 end                    
            end
            if startFrame==curStartingFrame
                curTrack.startingFrameExtra = curStartingFrame;
            end
            if endFrame==curEndingFrame
                curTrack.endingFrameExtra = curEndingFrame;
            end
        else
            for ii=curStartingFrame+1:curEndingFrame
                curImg = imgStack(:,:,ii);
                x = curTrack.xCoord(ii);
                y = curTrack.yCoord(ii);
                xi = round(x);
                yi = round(y);
                xRange = max(1,xi-halfWidth):min(xi+halfWidth,size(imgStack,2));
                yRange = max(1,yi-halfHeight):min(yi+halfHeight,size(imgStack,1));
                curAmpTotal = curImg(yRange,xRange);
                curAmpTotal = mean(curAmpTotal(:));
                curTrack.ampTotal(ii) =  curAmpTotal;
            end
        end
        if ~isempty(extraLengthForced) && abs(extraLengthForced)>0
            curTrack.startingFrameExtraExtra = curTrack.startingFrameExtra;
            curTrack.endingFrameExtraExtra = curTrack.endingFrameExtra;
            if curTrack.startingFrameExtra>1
                curTrack.startingFrameExtraExtra = max(1, curTrack.startingFrameExtra-extraLengthForced);
                x = curTrack.xCoord(curTrack.startingFrameExtra);
                y = curTrack.yCoord(curTrack.startingFrameExtra);
                xi = round(x);
                yi = round(y);
                xRange = max(1,xi-halfWidth):min(xi+halfWidth,size(imgStack,2));
                yRange = max(1,yi-halfHeight):min(yi+halfHeight,size(imgStack,1));
                for ii=curTrack.startingFrameExtraExtra:curTrack.startingFrameExtra
                    curImg = imgStack(:,:,ii);
                    curAmpTotal = curImg(yRange,xRange);
                    curAmpTotal = mean(curAmpTotal(:));
                    curTrack.xCoord(ii) = x;
                    curTrack.yCoord(ii) = y;
                    curTrack.ampTotal(ii) =  curAmpTotal;
    %                     curTrack.presence(ii) =  1;
                end
            end
            if curTrack.endingFrameExtra<numFrames
                curTrack.endingFrameExtraExtra = min(numFrames,curTrack.endingFrameExtra+extraLengthForced);            
                x = curTrack.xCoord(curTrack.endingFrameExtra);
                y = curTrack.yCoord(curTrack.endingFrameExtra);
                xi = round(x);
                yi = round(y);
                xRange = max(1,xi-halfWidth):min(xi+halfWidth,size(imgStack,2));
                yRange = max(1,yi-halfHeight):min(yi+halfHeight,size(imgStack,1));
                for ii=curTrack.endingFrameExtra:curTrack.endingFrameExtraExtra
                    curImg = imgStack(:,:,ii);
                    curAmpTotal = curImg(yRange,xRange);
                    curAmpTotal = mean(curAmpTotal(:));
                    curTrack.xCoord(ii) = x;
                    curTrack.yCoord(ii) = y;
                    curTrack.ampTotal(ii) =  curAmpTotal;
    %                     curTrack.presence(ii) =  1;
                end
            end
        end

%         curTrack.lifeTime = curTrack.endingFrameExtra - curTrack.startingFrameExtra;
    elseif attribute==2 || attribute==5 || attribute==6
        try
            startFrame = curTrack.startingFrameExtraExtra;
            endFrame = curTrack.endingFrameExtraExtra;
        catch
            try
                startFrame = curTrack.startingFrameExtra;
                endFrame = curTrack.endingFrameExtra;
            catch
                startFrame = max(1, curTrack.startingFrame-extraLengthForced);
                endFrame = min(numFrames,curTrack.endingFrame+extraLengthForced);
            end
        end
        if reTrack
            frameRange = startFrame:endFrame;
        else
            frameRange = [curTrack.startingFrameExtraExtra:curTrack.startingFrameExtra curTrack.endingFrameExtra:curTrack.endingFrameExtraExtra];
        end
        if attribute==2
            curTrack.forceMag = curTrack.amp;
        elseif attribute==5
            curTrack.ampTotal2 = curTrack.amp;
        elseif attribute==6
            curTrack.ampTotal3 = curTrack.amp;
        end
        
        for ii=frameRange
            curImg = imgStack(:,:,ii);
            x = curTrack.xCoord(ii);
            y = curTrack.yCoord(ii);
            xi = round(x);
            yi = round(y);
            xRange = max(1,xi-halfWidth):min(xi+halfWidth,size(curImg,2));
            yRange = max(1,yi-halfHeight):min(yi+halfHeight,size(curImg,1));
            curAmpTotal = curImg(yRange,xRange);
            if attribute==2
                curTrack.forceMag(ii) = mean(curAmpTotal(:));
            elseif attribute==5
                curTrack.ampTotal2(ii) = mean(curAmpTotal(:));
            elseif attribute==6
                curTrack.ampTotal3(ii) = mean(curAmpTotal(:));
            end
        end        
    elseif attribute==3 || attribute==4 %This time it uses FA area
        startFrame = curTrack.startingFrame;
        endFrame = curTrack.endingFrame;
        frameRange = startFrame:endFrame;
        for ii=frameRange
            curImg = imgStack(:,:,ii);
            if attribute==3
                curImg(curImg==0)=NaN; %assuming FA value
            end
%             if strcmp(curTrack.state(ii),'FA') || strcmp(curTrack.state(ii),'FC') % this is FA
            if curTrack.state(ii)==4 || curTrack.state(ii)==3 %'FA','FC') % this is FA
                pixelList=curTrack.FApixelList{ii};
                pixelIdxList = sub2ind(size(curImg),pixelList(:,2),pixelList(:,1));
                curAmpTotal = curImg(pixelIdxList);
            else
                x = curTrack.xCoord(ii);
                y = curTrack.yCoord(ii);
                xi = round(x);
                yi = round(y);
                xRange = max(1,xi-halfWidth):min(xi+halfWidth,size(curImg,2));
                yRange = max(1,yi-halfHeight):min(yi+halfHeight,size(curImg,1));
                curAmpTotal = curImg(yRange,xRange);
            end
            if attribute==3
                curTrack.fret(ii) = nanmean(curAmpTotal(:));
            elseif attribute==4
                curTrack.flowSpeed(ii) = mean(curAmpTotal(:));
            end
        end        
    elseif extrapolState
        disp('Please choose 1 or 2 for attribute.')
    end
    tracksNA(k)=curTrack;
    % progressbar(k/(numTracks), 0, 're-reading and tracking individual tracks');
    % progressbar(ii/(nFrames-1), 0, 'Matching with segmented adhesions:')
    % tk = toc;
    % waitbar(k/(numTracks), wtBar, sprintf([logMsg(0) timeMsg(tk*(numTracks/k)-tk)]));
%     parfor_progress;
    progressText(k/(numTracks));
end
end
%         for ii=startFrame:endFrame
%             curImg = imgStack(:,:,ii);
%             if ii<curStartingFrame
% %                 x = tracksNA(k).xCoord(curStartingFrame);
% %                 y = tracksNA(k).yCoord(curStartingFrame);
% %                 A = tracksNA(k).amp(curStartingFrame);
% %                 c = tracksNA(k).bkgAmp(curStartingFrame); 
%                 extrapolState=1;
%             elseif ii>curEndingFrame
% %                 x = tracksNA(k).xCoord(curEndingFrame);
% %                 y = tracksNA(k).yCoord(curEndingFrame);
%                 extrapolState=1;
%             else
%                 x = tracksNA(k).xCoord(ii);
%                 y = tracksNA(k).yCoord(ii);
%                 A = tracksNA(k).amp(curStartingFrame);
%                 c = tracksNA(k).bkgAmp(curStartingFrame); 
%                 extrapolState=0;
%             end
%             if attribute==1 && extrapolState %intensity
%                 xi = floor(x);
%                 xres = x-xi;
%                 yi = floor(y);
%                 yres = y-yi;
%                 curAmpTotal = curImg(yi,xi);
% %                 window = curImg(yi-w4:yi+w4, xi-w4:xi+w4);
% %                 [prmVect]=fitGaussianMixture2D(window,[x-xi+w4,y-yi+w4,A,sigma,c],'xyas');
% %                 [prmVect]=fitGaussianMixture2D(window,[x-xi+w4,y-yi+w4,curAmpTotal-min(window(:)),sigma,min(window(:))],'xyas');
%                 pstruct = fitGaussians2D(curImg, x, y, A, sigma, c, 'xyasc','WindowSize',8);
%                 if (abs(prmVect(1)-xres)<2 && abs(prmVect(2)-yres)<2 && prmVect(3)>0 && prmVect(3)<2*(curAmpTotal-min(window(:))))
%                     tracksNA(k).xCoord(ii) = xi + prmVect(1);
%                     tracksNA(k).yCoord(ii) = yi + prmVect(2);
%                     tracksNA(k).amp(ii) = prmVect(3);
%                     tracksNA(k).bkgAmp(ii) = prmVect(5);
%                     tracksNA(k).ampTotal(ii) =  prmVect(3)+prmVect(5);
%                 end
%                 tracksNA(k).ampTotal(ii) = curAmpTotal;
%             elseif attribute==1 && ~extrapolState %intensity
%                 xi = round(x);
%                 yi = round(y);
%                 xRange = xi-1:xi+1;
%                 yRange = yi-1:yi+1;
%                 curAmpTotal = curImg(yRange,xRange);
%                 tracksNA(k).ampTotal(ii) = mean(curAmpTotal(:));
% %                 curAmpTotal = curImg(yi,xi);
% %                 tracksNA(k).ampTotal(ii) = curAmpTotal;
%             elseif attribute==2 %forceMag
%                 xi = round(x);
%                 yi = round(y);
%                 xRange = xi-1:xi+1;
%                 yRange = yi-1:yi+1;
%                 curAmpTotal = curImg(yRange,xRange);
%                 tracksNA(k).forceMag(ii) = mean(curAmpTotal(:));
% %                 tracksNA(k).forceMag(ii) = curImg(round(y),round(x));
%             elseif extrapolState
%                 disp('Please choose 1 or 2 for attribute.')
%             end
%         end
%     end
% end
% end
