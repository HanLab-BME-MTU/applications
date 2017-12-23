%% pcCellRoiLeverMD - takes ROI masks from LEVER segmentation
% Nov. 2017
% mask taken from dirs.lever and saved to the local file

function [] = pcCellRoiLeverMD(MD,params,dirs)

load([dirs.tracking 'cellIdTYX.mat']);% cellTYX

nCells = length(cellTYX); %#ok<USENS>


%% ROI for cells
for icell = 1 : nCells
    roiLeverFname = [dirs.roiLever sprintf('%d',icell) '_roi.mat'];
    %     roiVisMovie = [dirs.roiVis sprintf('%d',icell) '_roi.avi'];
    
    if exist(roiLeverFname,'file') && ~params.always        
        if params.deepDebug
            fprintf(sprintf('cell %d roi LEVER: continue\n',icell)); 
        end
        
        continue;
    end   
    
    fprintf(sprintf('ROI LEVER cell %d/%d\n',icell,nCells));
    
    curCell = cellTYX{icell};    
    ntime = length(curCell.ts);
    
    curCellRoi = cell(1,ntime);
    
    %% movie init
    %     vwriter = VideoWriter(roiVisMovie,'Uncompressed AVI'); %#ok<TNMLP>
    %     vwriter.FrameRate = 10;
    %     open(vwriter);
    %     W = nan; H = nan;
    %%
    
    t0 = curCell.ts(1);
    for t = curCell.ts
            curT = t - t0 + 1;
            curX = round(curCell.xs(curT));            
            curY = round(curCell.ys(curT));
                        
            I = MD.getChannel(1).loadImage(curT);
            
            inLeverFname = [dirs.lever filesep reduceZeros4Lever(dirs.expname) '.LEVER' filesep num2str(t) '.mat'];            
            
            if ~exist(inLeverFname,'file')
                continue;
            end
            
            load(inLeverFname); % ROI_LEVER
            
            bby0 = max(1,curY - params.FOVRadius);
            bby1 = min(size(ROI_LEVER,1),curY + params.RoiRadius); % RoiRadius = 120 um
            bbx0 = max(1,curX - params.RoiRadius);
            bbx1 = min(size(ROI_LEVER,2),curX + params.RoiRadius);                                    
                        
            %
            %             MASK = false(size(ROI_LEVER));
            curI = I(bby0:bby1,bbx0:bbx1);
            MASK = ROI_LEVER(bby0:bby1,bbx0:bbx1);
            %             MASK = ROI_LEVER(bby0:bby1,bbx0:bbx1);
            
%             combinedImageBlur = imgaussfilt(combinedImage,2);
%             combinedImageBlurSized = imresize(combinedImageBlur,size(scores));
%             
%             curScores = combinedImageBlurSized(bby0:bby1,bbx0:bbx1);
%              
%             curYbb = curY - bby0;
%             curXbb = curX - bbx0;
%             
%             %             curYbbPatch = round(curYbb./params.patchSize);
%             %             curXbbPatch = round(curXbb./params.patchSize);
%             
%             % optimization: imresize by 1/patchSize
%             %             curScoresPatch = imresize(curScores,1.0/params.patchSize); % go to patch resolution
%             prc95 = prctile(curScores(:),95);
%             backgroundGL = prctile(curScores(:),params.detectionPPVicinityPrctile);
%             curScoresNorm = curScores - backgroundGL;
%             %             curScoresPatchBlur = imgaussfilt(curScoresPatchNorm,2);
%             
%             MASK = curScoresNorm > max(prc95,params.diffDetectionBackgroundScoreTH);            
%             centerMask = false(size(MASK)); centerMask(curYbb,curXbb) = true;
%             centerMask = imdilate(centerMask,strel('disk',round(15.0/params.pixelSize))); % minimum of 10 x 10 um mask
%             MASK = MASK | centerMask;
%             
%             %             %             imresize(curScoresPatch,size(curScores));
%             %
%             %             % Exclude:             MASK = otsuSeg(curScores, 2*params.patchSize);
%             %
%             %
%             %             % connected components
%             %             %             CC = bwconncomp(MASK);
%             %             %             ROI = true(size(MASK)); % default - take the whole FOV
%             %             %             [L,nL] = bwlabel(MASK);
%             %             %             for iL = 1 : nL
%             %             %                 curL = (L == iL);
%             %             %                 if curL(curYbb,curXbb)
%             %             %                     ROI = curL;
%             %             %                     break;
%             %             %                 end
%             %             %             end
            
            ROI = imdilate(MASK,strel('square',round(3.0/params.pixelSize)));
            
            curCellRoi{curT}.bby0 = bby0;
            curCellRoi{curT}.bby1 = bby1;
            curCellRoi{curT}.bbx0 = bbx0;
            curCellRoi{curT}.bbx1 = bbx1;
            
            curCellRoi{curT}.roi = ROI; 
            curCellRoi{curT}.img = curI; 
                        
            curCellRoi{curT}.x = curX;
            curCellRoi{curT}.y = curY;
            curCellRoi{curT}.t = t;
            
            %% movie frame
            %             curPerim = bwperim(ROI);
            %             curPerim = imdilate(curPerim,strel('square',5));
            %             Iprc99 = prctile(curI(:),99.99);
            %             curI(curPerim) = Iprc99;%uint16(65535);%intmax('uint16(65535)');
            %
            %             h = figure('visible','off'); imagesc(curI); colormap(gray);
            %             text(size(I,1)-300,size(I,2)-500,sprintf('%d minutes',round(t*params.timePerFrame)),'color','w','FontSize',15);
            %             haxes = get(h,'CurrentAxes');
            %             set(haxes,'XTick',[]);
            %             set(haxes,'XTickLabel',[]);
            %             set(haxes,'YTick',[]);
            %             set(haxes,'YTickLabel',[]);
            %
            %             drawnow; pause(0.01);
            %             movieFrame = getframe(h);
            %
            %             if isnan(W)
            %                 [H,W,~] = size(movieFrame.cdata);
            %                 minH = H;
            %                 maxH = H;
            %                 minW = W;
            %                 maxW = W;
            %             end
            %
            %             if H ~= size(movieFrame.cdata,1) || W ~= size(movieFrame.cdata,2)
            %                 minH = min(H,size(movieFrame.cdata,1));
            %                 maxH = max(H,size(movieFrame.cdata,2));
            %                 minW = min(W,size(movieFrame.cdata,1));
            %                 maxW = max(W,size(movieFrame.cdata,2));
            %             end
            %
            %             movieFrameResized = uint8(zeros(H,W,3));
            %             movieFrameResized(:,:,1) = imresize(movieFrame.cdata(:,:,1),[H,W]);
            %             movieFrameResized(:,:,2) = imresize(movieFrame.cdata(:,:,2),[H,W]);
            %             movieFrameResized(:,:,3) = imresize(movieFrame.cdata(:,:,3),[H,W]);
            %             movieFrame.cdata = movieFrameResized;
            %
            %             writeVideo(vwriter,movieFrame);
            %             close all;
            %             fprintf(sprintf('cell %d roi movie frame %d\n',icell,curT));
    end
    %     close(vwriter);
    %     clear vwriter;
    save(roiLeverFname,'curCellRoi');        
    
end
end

function reducedStr = reduceZeros4Lever(origStr)
if strcmp(origStr(end-1),'0')
    reducedStr = [origStr(1:end-2) origStr(end)];
else
    reducedStr = origStr;
end
end


