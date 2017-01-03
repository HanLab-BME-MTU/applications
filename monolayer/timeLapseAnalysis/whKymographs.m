function [] = whKymographs(params,dirs)

params.always = true;

params.xStepMinutes = 60;%240;
params.yStepUm = 50;
params.yMaxUm = params.kymoResolution.maxDistMu;%200;
params.fontsize = 24;

fprintf('start kymographs\n');
close all;
generateSpeedKymograph(params,dirs);
generateDirectionalMigrationKymograph(params,dirs);
generateCoordinatedMigrationKymograph(params,dirs);
% generateStrainRateKymograph(params,dirs);
% generateAccelerationKymograph(params,dirs);

close all;
end

function [] = generateSpeedKymograph(params,dirs)

speedKymographFname = [dirs.speedKymograph dirs.expname '_speedKymograph.mat'];

if exist(speedKymographFname,'file') && ~params.always
    return;
end

speedKymograph = nan(params.nstrips,params.nTime-params.frameJump);
speedKymographX = nan(params.nstrips,params.nTime-params.frameJump);
speedKymographY = nan(params.nstrips,params.nTime-params.frameJump);
for t = 1 : params.nTime - params.frameJump
    roiFname = [dirs.roiData pad(t,3) '_roi.mat']; % ROI
    mfFname = [dirs.mfData pad(t,3) '_mf.mat']; % dxs, dys
    
    load(roiFname);
    load(mfFname);
    
    speed = sqrt(dxs.^2 + dys.^2);
    DIST = bwdist(~ROI);
    
    for d = 1 : params.nstrips
        inDist = ...
            (DIST > (params.strips(d)-params.kymoResolution.stripSize)) & ...
            (DIST < params.strips(d)) & ...
            ~isnan(speed);
        speedInStrip = speed(inDist);
        speedKymograph(d,t) = mean(speedInStrip);
        % For directional migration
        speedInStripX = dxs(inDist);
        speedKymographX(d,t) = mean(abs(speedInStripX));
        speedInStripY = dys(inDist);
        speedKymographY(d,t) = mean(abs(speedInStripY));
    end
end

% Translate to mu per hour
speedKymograph = speedKymograph .* params.toMuPerHour;

save(speedKymographFname,'speedKymograph','speedKymographX','speedKymographY');

metaData.fname = [dirs.speedKymograph dirs.expname '_speedKymograph.eps'];
metaData.fnameFig = [dirs.speedKymograph dirs.expname '_speedKymograph.fig'];
metaData.caxis = [0 60];
% metaData.caxis = [8 25]; % Zhuo

% plotKymograph(speedKymograph,metaData,params);
params.caxis = metaData.caxis;
params.fname = metaData.fname;

plotKymograph(speedKymograph,params);

end

function [] = generateDirectionalMigrationKymograph(params,dirs)

directionalityKymographFname = [dirs.directionalityKymograph dirs.expname '_directionalityKymograph.mat'];

if exist(directionalityKymographFname,'file') && ~params.always
    return;
end


load([dirs.speedKymograph dirs.expname '_speedKymograph.mat']); % 'speedKymograph','speedKymographX','speedKymographY';
if params.isDx
    directionalityKymograph = speedKymographX ./ speedKymographY;
else
    directionalityKymograph = speedKymographY ./ speedKymographX;
end

save(directionalityKymographFname,'directionalityKymograph');

metaData.fname = [dirs.directionalityKymograph dirs.expname '_directionalityKymograph.eps'];
metaData.fnameFig = [dirs.directionalityKymograph dirs.expname '_directionalityKymograph.fig'];
metaData.caxis = [0 8];
% metaData.caxis = [0.9 1.4]; % Zhuo
metaData.caxis = [0 10]; % Georgio

% plotKymograph(directionalityKymograph,metaData,params);

params.caxis = metaData.caxis;
params.fname = metaData.fname;

plotKymograph(directionalityKymograph,params);
end


%%
function [] = generateCoordinatedMigrationKymograph(params,dirs)
coordinationKymographFname = [dirs.coordinationKymograph dirs.expname '_coordinationKymograph.mat'];

if exist(coordinationKymographFname,'file') && ~params.always
    return;
end

coordinationKymograph = nan(params.nstrips,params.nTime-params.frameJump);
for t = 1 : params.nTime - params.frameJump
    roiFname = [dirs.roiData pad(t,3) '_roi.mat']; % ROI
    coordinationFname = [dirs.coordination pad(t,3) '_coordination.mat']; % ROIclusters
    
    load(roiFname);
    load(coordinationFname);
        
    DIST = bwdist(~ROI);
    
    for d = 1 : params.nstrips        
        inDist = ((DIST > (params.strips(d)-params.kymoResolution.stripSize)) & (DIST < params.strips(d))) & ~isnan(ROIclusters);
                
        coordinationInStrip = ROIclusters(inDist);
        coordinationKymograph(d,t) = sum(coordinationInStrip)/length(coordinationInStrip);
    end
end

save(coordinationKymographFname,'coordinationKymograph');

metaData.fname = [dirs.coordinationKymograph dirs.expname '_coordinationKymograph.eps'];
metaData.fnameFig = [dirs.coordinationKymograph dirs.expname '_coordinationKymograph.fig'];
metaData.caxis = [0 1];

% plotKymograph(coordinationKymograph,metaData,params);

params.caxis = metaData.caxis;
params.fname = metaData.fname;

plotKymograph(coordinationKymograph,params);
end

%%
function [] = generateStrainRateKymograph(params,dirs)

strainRateKymographFname = [dirs.strainRateKymograph dirs.expname '_strainRateKymograph.mat'];

if exist(strainRateKymographFname,'file') && ~params.always
    return;
end

strainRateKymograph = nan(params.nstrips,params.nTime-params.frameJump);
for t = 1 : params.nTime - params.frameJump
    roiFname = [dirs.roiData pad(t,3) '_roi.mat']; % ROI
    strainRateFname = [dirs.strainRate pad(t,3) '_strainRate.mat']; % sr
    
    load(roiFname);
    load(strainRateFname);
        
    DIST = bwdist(~ROI);
    
    for d = 1 : params.nstrips        
        inDist = ((DIST > (params.strips(d)-params.kymoResolution.stripSize)) & (DIST < params.strips(d))) & ~isnan(sr);
                
        strainRateInStrip = sr(inDist);
        strainRateKymograph(d,t) = mean(strainRateInStrip);
    end
end

% Translate to mu per hour
strainRateKymograph = strainRateKymograph .* params.toMuPerHour ./ (2 * params.patchSize * params.pixelSize);

save(strainRateKymographFname,'strainRateKymograph');

metaData.fname = [dirs.strainRateKymograph dirs.expname '_strainRateKymograph.bmp'];
metaData.fnameFig = [dirs.strainRateKymograph dirs.expname '_strainRateKymograph.fig'];
metaData.caxis = [-0.2 0.2];

% plotKymograph(strainRateKymograph,metaData,params);

params.caxis = metaData.caxis;
params.fname = metaData.fname;

plotKymograph(strainRateKymograph,params);
end

function [] = generateAccelerationKymograph(params,dirs)

accelerationKymographFname = [dirs.accelerationKymograph dirs.expname '_accelerationKymograph.mat'];

if exist(accelerationKymographFname,'file') && ~params.always
    return;
end

accelerationKymograph = nan(params.nstrips,params.nTime-params.frameJump);
for t = 2 : params.nTime - params.frameJump - 1
    roiFname = [dirs.roiData pad(t,3) '_roi.mat']; % ROI
    accelerationFname = [dirs.acceleration pad(t,3) '_acceleration.mat']; % sr
    
    load(roiFname);
    load(accelerationFname); % acc
        
    DIST = bwdist(~ROI);
    
    for d = 1 : params.nstrips        
        inDist = ((DIST > (params.strips(d)-params.kymoResolution.stripSize)) & (DIST < params.strips(d))) & ~isnan(acc);
                
        accelerationInStrip = acc(inDist);
        accelerationKymograph(d,t) = mean(accelerationInStrip);
    end
end

% Translate to mu per hour
accelerationKymograph = accelerationKymograph .* params.toMuPerHour;

save(accelerationKymographFname,'accelerationKymograph');

metaData.fname = [dirs.accelerationKymograph dirs.expname '_accelerationKymograph.bmp'];
metaData.fnameFig = [dirs.accelerationKymograph dirs.expname '_accelerationKymograph.fig'];
metaData.caxis = [-4 8];

% plotKymograph(accelerationKymograph,metaData,params);

params.caxis = metaData.caxis;
params.fname = metaData.fname;

plotKymograph(accelerationKymograph,params);
end

function [] = plotKymographBackup(kymograph,metaData,params)
fprintf('start plot kymographs\n');
ntime = params.nTime - params.frameJump;
maxTime = ntime * params.timePerFrame;
maxTimeMod = (maxTime - mod(maxTime,100));
maxDistMod = (params.kymoResolution.max - mod(params.kymoResolution.max,100));

xTick = 1:(200/params.timePerFrame):((maxTimeMod/params.timePerFrame)+1);
xTickLabel = 0:200:maxTimeMod;
yTick = 1:(100/(params.pixelSize*params.patchSize)):((params.kymoResolution.maxDistMu/(params.pixelSize*params.patchSize))+1);
yTickLabel = 0:100:maxDistMod;

h = figure;
imagescnan(kymograph);
hold on;
caxis(metaData.caxis); colorbar;
haxes = get(h,'CurrentAxes');
set(haxes,'XLim',[1,maxTime/params.timePerFrame]);
set(haxes,'XTick',xTick);
set(haxes,'XTickLabel',xTickLabel);
set(haxes,'YTick',yTick);
set(haxes,'YTickLabel',yTickLabel);
set(haxes,'FontSize',32);
xlabel('Time (minutes)','FontSize',32); ylabel('Distance from edge (\mum)','FontSize',32);
hold off;
export_fig_biohpc(metaData.fname);
% eval(sprintf('print -dbmp16m  %s', metaData.fname));
% savefig(h,metaData.fnameFig);
fprintf('finish plot kymographs\n');
end