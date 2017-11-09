function [] = whTrajectories(params,dirs)

time = 1 : params.nTime - params.frameJump - 1; % -1 because segmentation is based on motion

fprintf('starting trajectories\n');

%% Bilateral filtering
for t = time
    bilateralFname = [dirs.mfBilateral pad(t,3) '_bilateral.mat'];
    
    if exist(bilateralFname,'file') && ~params.always
        return;
    end
    
    mfFname = [dirs.mfData pad(t,3) '_mf.mat']; % dxs, dys
    
    load(mfFname);
    
    [filterMfsDxs] = bilateralFiltering(dxs,params.patchSize,params.nBilateralIter);
    [filterMfsDys] = bilateralFiltering(dys,params.patchSize,params.nBilateralIter);
    save(bilateralFname,'filterMfsDxs','filterMfsDys');
end

%% Trajectories for visualization

trajVisFname = [dirs.trajectories dirs.expname '_trajectories.jpg'];

if ~exist(trajVisFname,'file') || params.always
    
    load([dirs.roiData pad(t,3) '_roi.mat']); % ROI
    [sizeY,sizeX] = size(ROI);
    DIST = bwdist(ROI);
    I0 = imread([dirs.images pad(1,3) '.tif']);
    
    
    trajectoriesVis = {};
    i = 1;
    for xx = 1 : (params.patchSize/params.pixelSize) : sizeX
        for yy = 1 : (params.patchSize/params.pixelSize) : sizeY
            x = round(xx);
            y = randi(sizeY);
            if ROI(y,x)
                trajectoriesVis{i}.ys = y;
                trajectoriesVis{i}.xs = x;
                i = i+1;
            end
        end
    end
    for t = time
        bilateralFname = [dirs.mfBilateral pad(t,3) '_bilateral.mat'];
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
    
    fprintf('display trajectories visualization\n');
    
    % Visualize cells' trajectories
    % h = figure('visible','off'); imagesc(I0); colormap(gray);
    h = figure; imagesc(I0); colormap(gray);
    hold on;
    for j = 1 : length(trajectoriesVis)
        path = trajectoriesVis{j};
        plot(path.xs(1),path.ys(1),'or','MarkerFaceColor','r');
        plot(path.xs(2:end),path.ys(2:end),'-r');
    end
    haxes = get(h,'CurrentAxes');
    set(haxes,'XTick',[]);
    set(haxes,'XTickLabel',[]);
    set(haxes,'YTick',[]);
    set(haxes,'YTickLabel',[]);
    hold off;
    saveas(h,trajVisFname);
    close(h);
end

%% Dense trajectories

fprintf('before dense trajectories\n');
trajFname = [dirs.trajectories dirs.expname '_trajectories.mat'];

if ~exist(trajFname,'file') || params.always
    
    allTrajectories = {};
    i = 1;
    for xx = 1 : params.patchSize : sizeX
        for yy = 1 : params.patchSize : sizeY
            x = xx;
            y = yy;
            if ROI(y,x)
                allTrajectories{i}.ys = y;
                allTrajectories{i}.xs = x;
                dist = DIST(round(allTrajectories{i}.ys),round(allTrajectories{i}.xs));
                allTrajectories{i}.dist = dist;
                i = i+1;
            end
        end
    end
    for t = time
        bilateralFname = [dirs.mfBilateral pad(t,3) '_bilateral.mat'];
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
    
    fprintf('before dense trajectories\n');
end
end

%% Utility functions
function [filterMfs] = bilateralFiltering(mfs,patchSize,nIterations)
maxVelocityOrig = max(abs(mfs(:)));
mfs =  mfs + maxVelocityOrig;
maxVelocityNew = max(mfs(:));
mfs =  mfs./maxVelocityNew;
lowResMfs = imresize(mfs,1.0/patchSize);
lowResMfs = max(lowResMfs,0);
lowResMfs = min(lowResMfs,1);

% bilateral filtering
filterMfs = lowResMfs;
for i = 1 : nIterations
    filterMfs = bfilter2(filterMfs,2);
end

filterMfs = filterMfs .* maxVelocityNew;
filterMfs = filterMfs - maxVelocityOrig;

% In case we want to go back to orinigal resolution
% newMfs = imresize(mfs,size(mfs));

end