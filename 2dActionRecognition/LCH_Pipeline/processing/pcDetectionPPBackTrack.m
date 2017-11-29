
function [] = pcDetectionPPBackTrack(MD,params,dirs)

% pcDetectionPPBackTrack
% Input: detected cells, processed cross-corrlation matching scores
% (include phase contrast) for forground and background
% Output: adding FN based on matching score forward- backward <= 0.1 
% Algorithm: backtrack each detection based on (1) motion field to previous
% frame (2) lap (linear assignment problem) (3) cells with no close match:
% find local maxima & check difference with background
%
% After this function all detections will be in detections.y and in
% detections.yBacktrack (marked in 'c' for debugging)

% Assaf Zaritsky, July 2015
if params.always
    unix(sprintf('rm %s filesep*_backtrack.mat',[dirs.detectPPData filesep]));
    unix(sprintf('rm %s filesep*_backtrack.eps',[dirs.detectPPVis filesep]));
    fprintf('detection backtracking (always): clean output directories\n'); 
end


% params.always = true;

diffDetectionBackgroundScoreTH = 0.1-0.05; % based on the data from Analysis/metaAnalysis
detectionPPVicinityPrctile = params.detectionPPVicinityPrctile;
detectionPPVicinityPacthes = params.detectionPPVicinityPatch;
   
if ~isfield(params,'sTime')
    params.sTime = 1;
end

for t = params.nTime-1 : -1 : params.sTime
    PPBacktrackFname0 = [dirs.detectPPData sprintf('%03d',t) '_backtrack.mat'];
    PPBacktrackisFname = [dirs.detectPPVis sprintf('%03d',t) '_backtrack.eps'];
    
    fprintf(sprintf('detection backtracking frame %d\n',t));  
    
    if exist(PPBacktrackFname0,'file') && ~params.always        
        if params.deepDebug
            fprintf('detection backtracking: file exist, continue\n');
        end;
        
        continue;
    end
    
    mfFname0 = [dirs.mfData sprintf('%03d',t) '_mf.mat'];
    load(mfFname0);
    dxs = -(dxs); dys = -(dys);
    
    PPFilteredFname0 = [dirs.detectPPData sprintf('%03d',t) '_filtered.mat'];
    PPBacktrackFname1 = [dirs.detectPPData sprintf('%03d',t+1) '_backtrack.mat'];
    %     PPFilteredFname1 = [dirs.detectPPData sprintf('%03d',t+1) '_filtered.mat'];
       
    if ~exist(PPBacktrackFname1,'file')
        % just for the first frame   
        load(PPFilteredFname0);
        detections.yBacktrack = [];
        detections.xBacktrack = [];
        detections.yBacktrackLowRes = [];
        detections.xBacktrackLowRes = [];
        save(PPBacktrackFname0,'detections','combinedImage','pstruct');
        continue;
    else    
        load(PPBacktrackFname1);
    end
    
    if ~isfield(detections,'yBacktrack')
        assert(PPBacktrackFname1); % why do we need this condition?
        detections.yBacktrack = [];
        detections.xBacktrack = [];
    end
    
    detections1 = detections;
    detections1.y = [detections1.y, detections1.yBacktrack];
    detections1.x = [detections1.x, detections1.xBacktrack];
    detections1.yLowRes = [detections1.yLowRes, detections1.yBacktrackLowRes];
    detections1.xLowRes = [detections1.xLowRes, detections1.xBacktrackLowRes];
        
    load(PPFilteredFname0);% detections,combinedImage,pstruct    
    
    if ~isfield(detections,'yBacktrack')
        detections.yBacktrack = [];
        detections.xBacktrack = [];
    else
        assert(false); % do we ever get here?
    end
          
    if params.deepDebug
            fprintf('detection backtracking: file exist, continue\n');
    end;
    
    detections.yBacktrackLowRes = [];
    detections.xBacktrackLowRes = [];
    detections.yBacktrackLowResFiltered = [];
    detections.xBacktrackLowResFiltered = [];
           
    % TODO: fix detections1 by dxs, dys        
    
    %% linear assignment problem
    matchingThresholdInPixels = 20/params.pixelSize; % 10 um
    %     X0 = [[detections.y, detections.yBacktrack]',[detections.x, detections.yBacktrack]'];
    X0 = [detections.y',detections.x'];
    X1 = [detections1.y',detections1.x'];
    n1 = size(X1,1);
    n0 = size(X0,1);
            
    if n1 == 0
        save(PPBacktrackFname0,'detections','combinedImage','pstruct');
        continue;
    end
    
    if n0 == 0
        warning(['no detections frame ' num2str(t)]);
        for d = 1 : n1
            curYLowRes = round(detections1.yLowRes(d));
            curXLowRes = round(detections1.xLowRes(d));
            if isBacktrack(curXLowRes,curYLowRes,combinedImage,diffDetectionBackgroundScoreTH,detectionPPVicinityPrctile,detectionPPVicinityPacthes)
                detections.yBacktrackLowRes = [detections.yBacktrackLowRes curYLowRes];
                detections.xBacktrackLowRes = [detections.xBacktrackLowRes curXLowRes];
            else
                detections.yBacktrackLowResFiltered = [detections.yBacktrackLowResFiltered curYLowRes];
                detections.xBacktrackLowResFiltered = [detections.xBacktrackLowResFiltered curXLowRes];
            end
        end
        
        % no detections for backtracking
        if n0 == 0
            save(PPBacktrackFname0,'detections','combinedImage','pstruct');
            continue;
        end
    else
        if params.deepDebug
            fprintf('detection backtracking: before LAP\n');
        end;
        
        % get all distances
        allDistances = pdist2([detections1.y',detections1.x'],[detections.y', detections.x']);
        
        D = createSparseDistanceMatrix(X1, X0, matchingThresholdInPixels);
        [link10, ~] = lap(D, [], [], 1);
        
        link10 = link10(1:n1);
        matchIdx = link10<=n0;
        %     idx0 = find(matchIdx);
        %     idx1 = double(link10(matchIdx));
        
        % Find detections in X1 but NOT in X0
        cellsInds1 = find(~matchIdx);%[idx1(~ismember(idx1,idx0)) []];
        
        if params.deepDebug
            fprintf('detection backtracking: after LAP\n');
        end;
        
        for d = 1 : length(cellsInds1)
            curInd = cellsInds1(d);
            coordYX = X1(curInd,:);
            
            curY = round(detections1.y(curInd));
            curX = round(detections1.x(curInd));
            
            % A match exists (Cell splitting event: 1 --> 2)
            if(min(allDistances(curInd,:)) < matchingThresholdInPixels)
                %
                %         allDists = sqrt((detections.y - curY).^2 + (detections.x - curX).^2) < matchingThresholdInPixels;
                %         if sum(allDists)
                save(PPBacktrackFname0,'detections','combinedImage','pstruct');
                continue;
            end
            
            
            curYLowRes = round(detections1.yLowRes(curInd));
            curXLowRes = round(detections1.xLowRes(curInd));
            
            %             BB = combinedImage(max(1,curYLowRes-detectionPPVicinityPacthes):min(size(combinedImage,1),curYLowRes+detectionPPVicinityPacthes),...
            %                 max(1,curXLowRes-detectionPPVicinityPacthes):min(size(combinedImage,2),curXLowRes+detectionPPVicinityPacthes));
            %
            %             matchScore = combinedImage(curYLowRes,curXLowRes);
            %             backgroundScore = prctile(BB(:),detectionPPVicinityPrctile);
            %             diffScore = matchScore - backgroundScore;
            
            if isBacktrack(curXLowRes,curYLowRes,combinedImage,diffDetectionBackgroundScoreTH,detectionPPVicinityPrctile,detectionPPVicinityPacthes) %diffScore > diffDetectionBackgroundScoreTH
                detections.yBacktrackLowRes = [detections.yBacktrackLowRes curYLowRes];
                detections.xBacktrackLowRes = [detections.xBacktrackLowRes curXLowRes];
            else
                detections.yBacktrackLowResFiltered = [detections.yBacktrackLowResFiltered curYLowRes];
                detections.xBacktrackLowResFiltered = [detections.xBacktrackLowResFiltered curXLowRes];
            end                        
        end
    end
    
    if params.deepDebug
            fprintf('detection backtracking: finished backtracking, before visualization\n');
    end;
    
    detections.yBacktrack = (detections.yBacktrackLowRes-0.5) .* params.patchSize;
    detections.xBacktrack = (detections.xBacktrackLowRes-0.5) .* params.patchSize;
    
    detections.yBacktrackFiltered = (detections.yBacktrackLowResFiltered-0.5) .* params.patchSize;
    detections.xBacktrackFiltered = (detections.yBacktrackLowResFiltered-0.5) .* params.patchSize;
    
   % visualization
    
    I = MD.getChannel(1).loadImage(t);
    
     % overlay detections on original image
    h = figure();
    imagesc(I); colormap(gray);    
    hold on
    plot(detections.x,detections.y,'g+','MarkerSize',5,'LineWidth',1.5);
    plot(detections.xBacktrack,detections.yBacktrack,'c+','MarkerSize',5,'LineWidth',1.5);    
    haxes = get(h,'CurrentAxes');
    set(haxes,'XTick',[]);
    set(haxes,'YTick',[]);
    set(haxes,'XTickLabel',[]);
    set(haxes,'YTickLabel',[]);
    set(h,'Color','w');
    axis equal;
    axis off;
    axis tight;
    hold off;
    export_fig_biohpc(PPBacktrackisFname);
    
    close all;
    
    
    save(PPBacktrackFname0,'detections','combinedImage','pstruct');
    
    if params.deepDebug
            fprintf('detection backtracking: done\n');
    end;
end
end

%%

function [res] = isBacktrack(curXLowRes,curYLowRes,combinedImage,diffDetectionBackgroundScoreTH,detectionPPVicinityPrctile,detectionPPVicinityPacthes)
BB = combinedImage(max(1,curYLowRes-detectionPPVicinityPacthes):min(size(combinedImage,1),curYLowRes+detectionPPVicinityPacthes),...
    max(1,curXLowRes-detectionPPVicinityPacthes):min(size(combinedImage,2),curXLowRes+detectionPPVicinityPacthes));

matchScore = combinedImage(curYLowRes,curXLowRes);
backgroundScore = prctile(BB(:),detectionPPVicinityPrctile);
diffScore = matchScore - backgroundScore;

res = diffScore > diffDetectionBackgroundScoreTH;
end