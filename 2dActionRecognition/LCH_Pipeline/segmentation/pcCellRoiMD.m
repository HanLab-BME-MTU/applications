%% pcCellRoiMD - calculates rough ROI mask for single cells
% Using Andrew's code
% Nov. 2017

function [] = pcCellRoiMD(MD,params,dirs)

load([dirs.tracking 'cellIdTYX.mat']);% cellTYX

nCells = length(cellTYX); %#ok<USENS>



%% ROI for cells
for icell = 1 : nCells
    roiFname = [dirs.roiData sprintf('%d',icell) '_roi.mat'];
    roiVisMovie = [dirs.roiVis sprintf('%d',icell) '_roi.avi'];
    
    if exist(roiFname,'file') && ~params.always        
        if params.deepDebug
            fprintf(sprintf('cell %d roi: continue\n',icell)); 
        end
        
        continue;
    end   
    
    fprintf(sprintf('ROI cell %d/%d\n',icell,nCells));
    
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
            
            detectPPFname = [dirs.detectPPData sprintf('%03d',t) '_detectionStats.mat'];
            load(detectPPFname); % combinedImage                        
            
            mfFname = [dirs.mfData sprintf('%03d',t) '_mf.mat'];% scores
            load(mfFname); % scores
            
            I = MD.getChannel(1).loadImage(curT);
            
            bby0 = max(1,curY - params.RoiRadius);
            bby1 = min(size(scores,1),curY + params.RoiRadius); % RoiRadius = 120 um
            bbx0 = max(1,curX - params.RoiRadius);
            bbx1 = min(size(scores,2),curX + params.RoiRadius);                                    
            
            curI = I(bby0:bby1,bbx0:bbx1);  
            [imgOut, mask, vI, fvI, bfvI] = segCellLCH(curI);
            
            MASK = mask;
            
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
    save(roiFname,'curCellRoi');
    
end
end


%%
function [imgOut, mask, vI, fvI, bfvI] = segCellLCH(img, varargin)
% simple script to pre-process/segment LCH cells for Deep learning.
% use example: [mask vI fvI bfvI]=segCellLCH(imread('./14-May-2017_atcc_s06_t120_x998_y1586_t130_f6.png'));
%
%  INPUT: image
%           'align' : orient the mask along major axis
%           'preview': plot results for debugging
%             'avgBG':    default is false for using zeros, 
%                         otherwise fills the non-mask with the avg of the 
%                         original image background
%             'varFilterOut': outputs the variance filter of the segmentated 
%                         image region, note variance filter is done first on entire 
%                         image the mask is cut from this
%             'centerMass": place object in the center of the image and pad
%             with zero/background
%             
%  OUTPUT: binary mask
%    (optional outputs include image processing steps for debugging)
%
%   
%
% by Andrew R. Jamieson, Oct 2017

ip = inputParser;
ip.addRequired('img', @isnumeric);
ip.addOptional('align',false, @islogical);
ip.addOptional('preview', false, @islogical);
ip.addOptional('avgBG', false, @islogical);
ip.addOptional('varFilterOut', false, @islogical);
ip.addOptional('centerMass', false, @islogical);
ip.parse(img,varargin{:});
p = ip.Results;

% Core image processing steps
%%%%%%%%%%%%%%%%%%%%%%%
I = mat2gray(img);
gI = imfilter(I, fspecial('gaussian', 5,.25));
vI = stdfilt(gI);
fvI = imfilter(vI, fspecial('gaussian', 5,1));
bfvI = imbinarize(fvI, .02);
bfvI = imdilate(bfvI, strel('disk', 5));
maskAll = imclose(bfvI, strel('disk', 5));
% maskAll = imerode(maskAll, strel('disk', 4));
% maskAll = imerode(maskAll, strel('disk', 2));
maskAll = bwfill(maskAll,'holes');
%%%%%%%%%%%%%%%%%%%%%%%

% I = mat2gray(img);
% gI = imfilter(I, fspecial('gaussian', 5,1));
% vI = stdfilt(gI);
% fvI = imfilter(vI, fspecial('gaussian', 7,2));
% bfvI = imbinarize(fvI, .02);
% bfvI = imdilate(bfvI, strel('disk', 2));
% maskAll = imclose(bfvI, strel('disk', 7));
% maskAll = bwfill(maskAll,'holes');

% See how many objects
CC = bwconncomp(maskAll);

% check if in center
% if yes, keep just this object
% if multiple, select one that overlaps with center
if CC.NumObjects > 1
    centerPts = round(size(maskAll)./2);
    mask = bwselect(maskAll,centerPts,centerPts);

    % check if mask present at center.
    if isempty(find(mask,1))
        % if not, just take the larget
        mask = bwareafilt(maskAll,1);
    end
else
    mask = maskAll;
end


% check if 97% covering image, then
% rp = regionprops(mask);
% if rp.Area/size(maskAll,1)^2  >= .9
%     mask = 0;
% end
if p.varFilterOut
    Ivar = rangefilt(I);
    imgFG = mask.*Ivar;    
%     varI = stdfilt(I);
else
    imgFG = mask.*I;    
end


if p.avgBG && ~p.varFilterOut % (note, just use zeros for var filter outputs
    rp = regionprops(mask);
    if rp.Area ~= size(maskAll,1)^2
        imgBG = ~mask.*I;
        imgBG = ~mask .* mean(imgBG(~mask));
        imgOut = imgFG + imgBG;
    else
        imgOut = imgFG;
    end
    
else
    imgOut = imgFG;
end

if p.align
    disp('re-orienting mask')
    rp = regionprops(mask,'orientation');
    rp.Orientation
    imgOut = imrotate(imgOut,-1*rp.Orientation);
end

if p.centerMass
    rp = regionprops(mask);
    
    [r c] = size(mask);
    rShift = round(r/2 - rp.Centroid(2))
    cShift = round(c/2 - rp.Centroid(1))

    % Call circshift to move region to the center.
    imgOut = circshift(imgOut, [rShift cShift]);
    mask = circshift(mask, [rShift cShift]);
end


if p.align
    disp('re-orienting mask')
    rp = regionprops(mask,'orientation');
    rp.Orientation
    imgOut = imrotate(imgOut,-1*rp.Orientation);
end


if p.preview 
    % preview results
    figure; imshow(imgOut);
    figure;
    subplot(3,2,1);imshow(bfvI,[]); title('after thre & dilation');
    subplot(3,2,2);imshow(vI,[]);title('variance filter');
    subplot(3,2,3);imshow(fvI,[]); title('smoothing');
    subplot(3,2,4);imshow(imgOut,[]);title('final');
    subplot(3,2,5);imshow(maskAll,[]);title('all object mask');
    subplot(3,2,6);imshow(I,[]); title('original');
end

% Works roughly for 128x128 downsampling
% vI = stdfilt(I);
% fvI = imfilter(vI, fspecial('gaussian', 7,3));
% bfvI = imbinarize(fvI, .02);
% bfvI = imdilate(bfvI, strel('disk', 3));
% mask = imclose(bfvI, strel('disk', 7));
end