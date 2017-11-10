function [] = whTrajecotires(params,dirs)

time = 1 : params.nTime - params.frameJump - 1; % -1 because segmentation is based on motion

%% Bilateral filtering
for t = time   
    accelerationFname = [dirs.acceleration pad(t,3) '_acceleration.mat'];
    
    if exist(accelerationFname,'file') && ~params.always
        continue;
    end
    
    mfFname = [dirs.mfData pad(t,3) '_mf.mat']; % dxs, dys
        
    load(mfFname);
    
    bilateralFname = [bilateralDir prefix int2str(t) '_bilateral.mat'];
    
    if exist(bilateralFname,'file') && ~params.always
        load(bilateralFname);
    else
        [filterMfsDxs] = bilateralFiltering(dxs,patchSize,nBilateralIter);
        [filterMfsDys] = bilateralFiltering(dys,patchSize,nBilateralIter);
        save(bilateralFname,'filterMfsDxs','filterMfsDys');
    end
end

%% Trajectories for visualization

load([dirs.roiData pad(t,3) '_roi.mat']); % ROI
[sizeY,sizeX] = size(ROI);
DIST = bwdist(ROI);
I0 = imread([dirs.images pad(t,3) '.tif']);


trajectoriesVis = {};
i = 1;
for xx = 1 : (params.patchSize/params.pixelSize) : sizeX
    for yy = 1 : 1 : sizeY/70
        x = xx;
        y = randi(sizeY);
        if ROI(y,x)
            trajectoriesVis{i}.ys = y;
            trajectoriesVis{i}.xs = x;           
            i = i+1;
        end
    end
end
for t = time
    bilateralFname = [bilateralDir prefix int2str(t) '_bilateral.mat'];
    load(bilateralFname); % filterMfsDxs, filterMfsDys
    for i = 1 : length(trajectoriesVis)
        curY = round(trajectoriesVis{i}.ys(t));
        curX = round(trajectoriesVis{i}.xs(t));
        
        filterMfsDys1 = imresize(filterMfsDys,size(ROI),'nearest');
        filterMfsDxs1 = imresize(filterMfsDxs,size(ROI),'nearest');
        
        nextY = curY + filterMfsDys1(curY,curX);
        nextX = curX + filterMfsDxs1(curY,curX);
        % borders
        nextY = min(sizeY,max(nextY,1));
        nextX = min(sizeX,max(nextX,1));
        
        trajectoriesVis{i}.ys = [trajectoriesVis{i}.ys nextY];
        trajectoriesVis{i}.xs = [trajectoriesVis{i}.xs nextX];
    end
end
  
% Visualize cells' trajectories
h = figure('visible','off'); imagesc(I0); colormap(gray);
hold on;
for j = 1 : length(trajectoriesVis)
    path = trajectoriesVis{j}.path;
    plot(path.xs(1),path.ys(1),'or','MarkerFaceColor','r');
    plot(path.xs(2:end),path.ys(2:end),'-r');
end
haxes = get(h,'CurrentAxes');
set(haxes,'XTick',[]);
set(haxes,'XTickLabel',[]);
set(haxes,'YTick',[]);
set(haxes,'YTickLabel',[]);
hold off;
trajVisFname = [dirs.trajectories dirs.expname '_trajectories.bmp'];
saveas(h,trajVisFname);
close(h);


%% Dense trajectories
allTrajectories = {};
i = 1;
for xx = 1 : params.patchSize : sizeX
    for yy = 1 : params.patchSize : sizeY
        x = xx;
        y = yy;
        if ROI(y,x)
            allTrajectories{i}.ys = y;
            allTrajectories{i}.xs = x;
            dist = DIST(round(trajectories{i}.ys),round(trajectories{i}.xs));
            allTrajectories{i}.dist = dist;
            i = i+1;
        end
    end
end
for t = time
    bilateralFname = [bilateralDir prefix int2str(t) '_bilateral.mat'];
    load(bilateralFname); % filterMfsDxs, filterMfsDys
    for i = 1 : length(allTrajectories)
        curY = round(allTrajectories{i}.ys(t));
        curX = round(allTrajectories{i}.xs(t));
        
        filterMfsDys1 = imresize(filterMfsDys,size(ROI),'nearest');
        filterMfsDxs1 = imresize(filterMfsDxs,size(ROI),'nearest');
        
        nextY = curY + filterMfsDys1(curY,curX);
        nextX = curX + filterMfsDxs1(curY,curX);
        % borders
        nextY = min(sizeY,max(nextY,1));
        nextX = min(sizeX,max(nextX,1));
        
        allTrajectories{i}.ys = [allTrajectories{i}.ys nextY];
        allTrajectories{i}.xs = [allTrajectories{i}.xs nextX];
    end
end

ncells = length(allTrajectories);
    PERS.persi = zeros(1,ncells);
    PERS.dists = zeros(1,ncells);
    for i = 1 : ncells
        xs = allTrajectories{i}.xs;
        ys = allTrajectories{i}.ys;
        allTrajectories{i}.persi = chemotacticIndex(xs(1:end),ys(1:end));
        PERS.persi(i) = allTrajectories{i}.persi;
        PERS.dists(i) = allTrajectories{i}.dist;
    end

save(trajFname,'allTrajectories','trajectoriesVis','I0','ROI','DIST','PERS','ncells');
end