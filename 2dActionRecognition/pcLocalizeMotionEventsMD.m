%%

function [] = pcLocalizeMotionEventsMD(MD,params,dirs)

windowSize = round(35 / params.pixelSize); % 40 um
lbpMapping = getmapping(8,'riu2');

for t = 2 : params.nTime - params.frameJump - 1
    
    fprintf(sprintf('localize motion events frame %d\n',t));
    
    detectFname = [dirs.detectData pad(t,3) '_detections.mat'];
    cellsFname = [dirs.detectData pad(t,3) '_cells.mat'];
    
    if exist(cellsFname,'file') && ~params.always
        continue;
    end
    
    load(detectFname); % 'I','detections','combinedScores','BIN_SCORES','socresLocalMaxima'
    
    %     I = imread([dirs.images pad(t,3) '.tif']);
    I = MD.getChannel(1).loadImage(t);
    
    % Assumes one channel or 3 channels with same value
    if size(I,3) > 1
        I = I(:,:,1);
    end
    
    [sizeY,sizeX] = size(I);
    
    nDetections = length(detections.xs);
    ncells = 0;
    
    Idisp = uint16(zeros(size(I,1),size(I,2),3));
    Idisp(:,:,1) = I;
    Idisp(:,:,2) = I;
    Idisp(:,:,3) = I;
    
    BIN = false(size(I));
    
    cells.xs = [];
    cells.ys = [];
    for icell = 1 : nDetections
        x = detections.xs(icell);
        y = detections.ys(icell);
        
        bbxs = max(1,x-windowSize):min(sizeX,x+windowSize);
        bbys = max(1,y-windowSize):min(sizeY,y+windowSize);
        
        % Ugly patch because detections were found out side of the image!?
        if isempty(bbxs) || isempty(bbys)
            continue;
        end
        
        scoresCrop = combinedScores(bbys,bbxs);
        %         % TODO: maybe replace this code by code that segments the motion
        %         % maps?
        %         [BIN0] = segmentPhaseContrastLBPKmeans(Icrop,10,lbpMapping);
        %         BIN0 = OtsuBinarization(scoresCrop);
        TRI = otsuMultTH(scoresCrop,3);
        BIN0 = TRI == 1;
        BIN1 = localizeCellSegmentation(BIN0,params);
        BIN(bbys,bbxs) = BIN1;
        if sum(BIN1(:)) > 0
            ncells = ncells + 1;
            [cogX,cogY] = cog(BIN1);
            cogX = cogX + max(1,x-windowSize);
            cogY = cogY + max(1,y-windowSize);
            cells.xs = [cells.xs cogX];
            cells.ys = [cells.ys cogX];
            cells.seg{ncells} = BIN1;
            cells.bb{ncells}.bbxs = bbxs;
            cells.bb{ncells}.bbys = bbys;
        end
    end
    Idisp(:,:,1) = uint16(min(uint32(Idisp(:,:,1)) +  uint32(BIN .* 10000),uint32(intmax('uint16'))));
    
    h = figure('visible','off');
    imagesc(Idisp);
    hold on;
    plot(detections.xs,detections.ys,'sw','MarkerFaceColor','w','MarkerSize',3);
    plot(cells.xs,cells.ys,'sg','MarkerFaceColor','g','MarkerSize',2);
    
    text(size(I,1)-150,size(I,2)-500,sprintf('%d minutes',round(t*params.timePerFrame)),'color','w','FontSize',15);
    haxes = get(h,'CurrentAxes');
    set(haxes,'XTick',[]);
    set(haxes,'XTickLabel',[]);
    set(haxes,'YTick',[]);
    set(haxes,'YTickLabel',[]);
    hold off;        
    saveas(h,[dirs.detectVis pad(t,3) '_singleCells.jpg']);    
    close(h);   
    
    save(cellsFname,'I','cells');
end
end
%%

function [BIN1] = localizeCellSegmentation(BIN0,params)
y = round(size(BIN0,1) ./ 2);
x = round(size(BIN0,2) ./ 2);

BIN1 = false(size(BIN0));

BIN0 = morphOpen(BIN0,round(5./params.pixelSize));
BIN0 = morphClose(BIN0,round(5./params.pixelSize));
BIN0 = imfill(BIN0,'holes');

[L,nL] = bwlabel(BIN0,8);
for i = 1 : nL
    curBin = L == i;
    if curBin(y,x)
        cellSize = sum(curBin(:));
        if isSingleCell(curBin,params)
            BIN1 = curBin;
        end
        break;
    end
end
end

%% true - single cell, false -- multiple cells
function [val] = isSingleCell(BIN,params)
val =  sum(BIN(:)) > params.cellMinArea && sum(BIN(:)) < params.cellMaxArea;
end