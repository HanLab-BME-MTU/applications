function [] = whKymographsStd(params,dirs)

params.always = false;

fprintf('start kymographs\n');
close all;
generateSpeedKymographStd(params,dirs);
generateDirectionalMigrationKymographStd(params,dirs);

close all;
end

function [] = generateSpeedKymographStd(params,dirs)

speedKymographFname = [dirs.speedKymographStd dirs.expname '_speedKymographStd.mat'];

if exist(speedKymographFname,'file') && ~params.always
    return;
end

speedKymographStd = nan(params.nstrips,params.nTime-params.frameJump);
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
        speedKymographStd(d,t) = std(speedInStrip);        
    end
end

% Translate to mu per hour
speedKymographStd = speedKymographStd .* params.toMuPerHour;

save(speedKymographFname,'speedKymographStd');

metaData.fname = [dirs.speedKymographStd dirs.expname '_speedKymographStd.eps'];
metaData.fnameFig = [dirs.speedKymographStd dirs.expname '_speedKymographStd.fig'];
metaData.caxis = [0 20];

plotKymograph(speedKymographStd,metaData,params);

end

function [] = generateDirectionalMigrationKymographStd(params,dirs)

directionalityKymographFname = [dirs.directionalityKymographStd dirs.expname '_directionalityKymographStd.mat'];

if exist(directionalityKymographFname,'file') && ~params.always
    return;
end

directionalityKymographStd = nan(params.nstrips,params.nTime-params.frameJump);

for t = 1 : params.nTime - params.frameJump
    roiFname = [dirs.roiData pad(t,3) '_roi.mat']; % ROI
    mfFname = [dirs.mfData pad(t,3) '_mf.mat']; % dxs, dys
    
    load(roiFname);
    load(mfFname);
        
    DIST = bwdist(~ROI);
    
    for d = 1 : params.nstrips
        inDist = ...
            (DIST > (params.strips(d)-params.kymoResolution.stripSize)) & ...
            (DIST < params.strips(d)) & ...
            ~isnan(dxs) & ~isnan(dys);
        
        speedInStripX = dxs(inDist);        
        speedInStripY = dys(inDist);
        
        directInStrip = abs(speedInStripY ./ speedInStripX);
        directInStrip = directInStrip(~isnan(directInStrip));
        directInStrip = directInStrip(directInStrip <= 8);        
        
        directionalityKymographStd(d,t) = std(directInStrip);
    end
end

save(directionalityKymographFname,'directionalityKymographStd');

metaData.fname = [dirs.directionalityKymographStd dirs.expname '_directionalityKymographStd.eps'];
metaData.fnameFig = [dirs.directionalityKymographStd dirs.expname '_directionalityKymographStd.fig'];
metaData.caxis = [0 1];

plotKymograph(directionalityKymographStd,metaData,params);

end


function [] = plotKymograph(kymograph,metaData,params)
fprintf('start plot std kymographs\n');
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
colormap 'jet';
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