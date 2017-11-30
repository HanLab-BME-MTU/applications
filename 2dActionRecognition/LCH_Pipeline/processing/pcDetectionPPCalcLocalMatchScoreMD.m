
function [] = pcDetectionPPCalcLocalMatchScoreMD(params,dirs)

% pcDetectionPPCalcLocalMatchScoreMD
% Input: detected cells, processed cross-corrlation matching scores (include phase contrast)
% Output for each detection: 
%   matching score, 30% (parameter) matching score in ROI, var (matching score), 
% %   filtered detection in (-5) - 5 (parameter) time frames ahead 
% Statistics are going to be collectde across all movies to decide on
% strategy for further filtering.

% Assaf Zaritsky, June 2015

if params.always
    unix(sprintf('rm %s*_detectionStats.mat',[dirs.detectPPData filesep]));    
    fprintf('cell detection (always): clean output directories\n'); 
end


% what percentile to use as "background" in cross correlation matching score 
detectionPPVicinityPrctile = params.detectionPPVicinityPrctile; 
% % how many farmes to look back and forth to filter out a detection because
% % the same hit was filtered earlier
% detectionPPFilterFPTime = params.detectionPPFilterFPTime;

% vicinity (in pixels) around a detection to calculate statistics at
% detectionPPVicinityPacthes = ceil(params.detectionPPVicinityPixels / params.pixelSize);
detectionPPVicinityPacthes = params.detectionPPVicinityPatch;
    
if ~isfield(params,'sTime')
    params.sTime = 1;
end

for t = params.sTime : params.nTime - params.frameJump
    detectPPFname = [dirs.detectPPData sprintf('%03d',t) '_detectionStats.mat'];
    
    fprintf(sprintf('cell post processing detection frame %d\n',t));
    
    if exist(detectPPFname,'file') && ~params.always
        if params.deepDebug
            fprintf('detection PP: existing file, continue\n');
        end;
        
        continue;
    end
    
    detectFname = [dirs.detectData sprintf('%03d',t) '_detections.mat'];
    mfFname = [dirs.mfData sprintf('%03d',t) '_mf.mat'];                    
    
    if params.deepDebug
            fprintf('detection PP: before loading\n');
    end;
    
    load(detectFname); % 'detections','combinedImage','pstruct'
    
    %     % temporary
    %     if isfield(detections,'LowRes')
    %         copyfile(detectFname,detectPPFname);
    %         continue;
    %     end
        
    
    load(mfFname); % dxs, dys
    
    if params.deepDebug
            fprintf('detection PP: after loading\n');
    end;
    
    n = length(detections.xLowRes);
        
    detections.LowRes.matchScore = nan(1,n);
    detections.LowRes.backgroundScore = nan(1,n);
    detections.LowRes.varScore = nan(1,n);
    detections.LowRes.diffScore = nan(1,n);
    detections.LowRes.speed = nan(1,n);
    
    for d = 1 : n
        if params.deepDebug
            fprintf(sprintf('detection PP: processing detection %d\n',d));
        end;
        curX = round(detections.xLowRes(d));
        curY = round(detections.yLowRes(d));
        
        BB = combinedImage(max(1,curY-detectionPPVicinityPacthes):min(size(combinedImage,1),curY+detectionPPVicinityPacthes),...
            max(1,curX-detectionPPVicinityPacthes):min(size(combinedImage,2),curX+detectionPPVicinityPacthes));
        
        detections.LowRes.matchScore(d) = combinedImage(curY,curX);
        detections.LowRes.backgroundScore(d) = prctile(BB(:),detectionPPVicinityPrctile);
        detections.LowRes.varScore(d) = var(BB(:));
        detections.LowRes.diffScore(d) = detections.LowRes.matchScore(d) - detections.LowRes.backgroundScore(d);
        xx=round(detections.x(d));
        yy=round(detections.y(d));
        detections.LowRes.speed(d) = sqrt(dxs(yy,xx).^2 + dys(yy,xx).^2);                
    end
         
    % TODO: find filtered detections over time
    
    save(detectPPFname,'detections','combinedImage','pstruct');
    
    if params.deepDebug
        fprintf('detection PP: done\n');
    end;
end
end