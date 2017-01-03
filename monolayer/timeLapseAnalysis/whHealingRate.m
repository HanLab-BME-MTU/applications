function [] = whHealingRate(params,dirs)

healingRateFname = [dirs.roiVis dirs.expname '_healingRate.eps'];
healingRateMetaFname = [dirs.healingRate dirs.expname '_healingRate.mat'];

% if exist(healingRateFname,'file') && exist(healingRateMetaFname,'file') && ~params.always
%     return;
% end

time = 1 : params.nTime - params.frameJump - 1; % -1 because segmentation is based on motion
ntime = length(time);

healingRate = nan(1,ntime);
averageHealingRate = nan(1,ntime);

load([dirs.roiData pad(1,3) '_roi.mat']); % ROI
sumInitROI = sum(ROI(:)); clear ROI;

fprintf('calculating healing rate\n');

for t = time
    load([dirs.roiData pad(t,3) '_roi.mat']); % ROI
    ROI0 = ROI; clear ROI;
    load([dirs.roiData pad(t+1,3) '_roi.mat']); % ROI
    ROI1 = ROI; clear ROI;
    
    nDiffPixels = sum(ROI1(:)) - sum(ROI0(:));
    nDiffPixelsMeta = sum(ROI1(:)) - sumInitROI;
    if params.isDx
        healingRate(t) = params.toMuPerHour * nDiffPixels / size(ROI0,1);
        averageHealingRate(t) = params.toMuPerHour * nDiffPixelsMeta / (size(ROI0,1) * t);
    else
        warning('currently supporting only dx');
        healingRate(t) = params.toMuPerHour * nDiffPixels / size(ROI0,2);
        averageHealingRate(t) = params.toMuPerHour * nDiffPixelsMeta / (size(ROI0,2) * t);
    end
end

maxTime = ntime * params.timePerFrame;
maxTimeMod = (maxTime - mod(maxTime,100));
maxSpeed = 80; % um / hour

h = figure; 
hold on;
plot(time*params.timePerFrame,min(healingRate,maxSpeed),'or','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',10);
xlabel('Time (minutes)','FontSize',32); 
ylabel('Healing rate (\mum hour{-1})','FontSize',32);
haxes = get(h,'CurrentAxes');
set(haxes,'XLim',[0,maxTime]);
set(haxes,'XTick',0:200:maxTimeMod);
set(haxes,'XTickLabel',0:200:maxTimeMod);
set(haxes,'YLim',[0,maxSpeed]);
set(haxes,'YTick',0:20:maxSpeed);
set(haxes,'YTickLabel',0:20:maxSpeed);
set(haxes,'FontSize',32);
hold off;
% eval(sprintf('print -dbmp16m  %s', healingRateFname));

% for building package
% if isunix
% export_fig_biohpc(healingRateFname); % change 
% else
print(healingRateFname, '-dpdf');
% end

save(healingRateMetaFname,'averageHealingRate','healingRate');
end