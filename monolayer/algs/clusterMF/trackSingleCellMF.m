
% Tracks all cells tracks using filtered motion fields
function [] = trackSingleCellMF(dirname,prefix,patchSize,t0, phase1, phase2, phase3, padSize, always)

motionDir = strcat(dirname,'/MF/');
clustersDir = sprintf('%s/clusters/',motionDir); % bilateral filtering results in there...
tracksDir = sprintf('%s/tracks/',motionDir);
imgDir = strcat(dirname,'/images/');
roiDir = strcat(dirname,'/ROI/');

outputDataFname = sprintf('%s/%s_trajectories.mat',tracksDir,prefix);

if (exist(outputDataFname,'file') ~= 0 && ~always)
    return;
end

if (exist(sprintf('%s',clustersDir),'dir') == 0)
    eval(sprintf('mkdir %s',clustersDir));
end

if (exist(sprintf('%s',tracksDir),'dir') == 0)
    eval(sprintf('mkdir %s',tracksDir));
end

I = imread(sprintf('%s/%s%s.tif',imgDir,prefix,pad(t0,padSize)));
[sizeY,sizeX] = size(I);
load(sprintf('%s/%s%s_ROI',roiDir,prefix,pad(t0,padSize))); % binNoExcesses
ROI = binNoExcesses;


trajectories = {};
i = 1;
for y = 1 : 4*patchSize : sizeY
    for x = 1 : 1*patchSize : sizeX
        if ~ROI(y,x)
            trajectories{i}.path.ys = y;
            trajectories{i}.path.xs = x;
            i = i+1;
        end
    end
end

for t = t0+1 : phase3-1
    clustersInfoFname = sprintf('%s/%s%s_clustersInfo_%d.mat',clustersDir,prefix,pad(t,padSize),patchSize);
    load(clustersInfoFname);
    for i = 1 : length(trajectories)
        curY = round(trajectories{i}.path.ys(t - t0));
        curX = round(trajectories{i}.path.xs(t - t0));
        
        nextY = curY + outImgDy1(curY,curX);%/2.0;
        nextX = curX + outImgDx1(curY,curX);%/2.0;
        % borders
        nextY = min(sizeY,max(nextY,1));
        nextX = min(sizeX,max(nextX,1));
        
        trajectories{i}.path.ys(t - t0 + 1) = nextY;
        trajectories{i}.path.xs(t - t0 + 1) = nextX;
    end
end

phase1 = phase1 - t0 - 1;
phase2 = phase2 - t0 - 1;
phase3 = phase3 - t0 - 1;

I = imread(sprintf('%s/%s%s.tif',imgDir,prefix,pad(t0,padSize)));
Itrack = I;
h = figure; imagesc(Itrack); colormap(gray);
hold on;
for i = 1 : length(trajectories)
    plot(trajectories{i}.path.xs(1:phase1),trajectories{i}.path.ys(1:phase1),'g');
    plot(trajectories{i}.path.xs(phase1+1:phase2),trajectories{i}.path.ys(phase1+1:phase2),'r');
    plot(trajectories{i}.path.xs(phase2+1:end),trajectories{i}.path.ys(phase2+1:end),'b');
end
hold off;

outputImgFname = sprintf('%s/%s_trajectories.jpg',tracksDir,prefix);
eval(sprintf('print -djpeg %s', outputImgFname));
close(h);
save(outputDataFname,'trajectories', 'I');

%% AVI of tracks

aviFilename = sprintf('%s/%s.avi',tracksDir,prefix);
aviobj = avifile(aviFilename,'fps',3,'compression','None');

for f = t0+1 : phase3
    curFrame = f - t0;
    I = imread(sprintf('%s/%s%s.tif',imgDir,prefix,pad(f,padSize)));
    try
        h = figure; colormap(gray); imagesc(I);
        caxis([0 255]);
        hold on;
        for i = 1 : length(trajectories)
            plot(trajectories{i}.path.xs(1:phase1),trajectories{i}.path.ys(1:phase1),'g');
            plot(trajectories{i}.path.xs(phase1+1:phase2),trajectories{i}.path.ys(phase1+1:phase2),'r');
            plot(trajectories{i}.path.xs(phase2+1:end),trajectories{i}.path.ys(phase2+1:end),'b');
        end
        hold off;
        drawnow;
        movieFrame = getFrame(h);
        aviobj = addframe(aviobj, movieFrame);
        close all;
    catch
        fprintf(sprintf('No image %d',f));
        break;
    end
end

aviobj = close(aviobj);
close all;

end