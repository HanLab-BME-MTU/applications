function [] = whSeedsOfPlithotaxis(params,dirs)

fprintf('Starting seeds of plithotaxis\n');

% plithotaxisFname = [dirs.plithotaxis dirs.expname '_plithotaxis.eps'];
plithotaxisFname = [dirs.plithotaxis dirs.expname '_plithotaxis'];

% time = 1 : params.nTime - params.frameJump - 1; % -1 because segmentation is based on motion
time = 1 : 91 - params.frameJump - 1; % -1 because segmentation is based on motion
ntime = length(time);

load([dirs.roiData pad(1,3) '_roi.mat']); % ROI
[ySize,xSize] = size(ROI);

if exist([plithotaxisFname '.mat'],'file') &&  ~params.always
    outputSeeds(plithotaxisFname,params,ntime,ySize);
    return;
end

seeds = false(ySize,ntime);
seedsLocation = nan(ySize,ntime);

fprintf('Seeds of plithotaxis: before loop\n');
for t = time
    load([dirs.roiData pad(t,3) '_roi.mat']); % ROI
    ROI0 = ROI; clear ROI;
    load([dirs.roiData pad(t+1,3) '_roi.mat']); % ROI
    ROI1 = ROI; clear ROI;
    
    DIFF = ROI1 & ~ROI0;
    [ys,xs] = find(DIFF);
    
    %     assert(length(ys)==length(unique(ys)));
    
    seeds(ys,t) = true;    
    seedsLocation(ys,t) = xs;        
end
fprintf('Seeds of plithotaxis: after loop\n');

% strainEvents, strainEventsOrigResolution, motionEvents, ncellsLeaders, ncellsEffected, numEventsYAxis, seedsVis
strainEventsOutput = getStrainEvents(seeds,seedsLocation,params.patchSize); % output.nStrainEventsPercentage - measure for Meta Analysis!


fprintf('Seeds of plithotaxis: after getStrainEvents\n');
%% evolution of contour over time
% intXSpace = 5;
% space = 20;
% contourEvolution = zeros(ySize,intXSpace + space * time(end) + xSize);
% % figure; imagescnan(contourEvolution);
% % hold on;
% for t = timecumsumFlow
%     load([dirs.roiData pad(t,3) '_roi.mat']); % ROI
%     boundaries = bwboundaries(ROI);
%     assert(length(boundaries) == 1);
%     bws = boundaries{1};
%     inds = bws(:,1) > 1 & bws(:,2) > 1 & bws(:,1) < size(ROI,1);
%     bwys = bws(inds,1)';
%     bwxs = bws(inds,2)';
%     %     bwxs = bwxs' - min(bwxs);
%     
%     contourEvolutionFlat = contourEvolution(:);
%     contourEvolutionFlat(bwys + (intXSpace + bwxs + space*t - 1)*size(contourEvolution,1)) = 50 + t;
%     contourEvolution = reshape(contourEvolutionFlat,size(contourEvolution));
%     
%     %     for i = 1 : length(bwxs)
%     %         contourEvolution(bwys(i),intXSpace + bwxs(i) + space*t) = 50 + t;
%     %     end
%     %     contourEvolution(bwys',intXSpace + bwxs' + space*t) = true;
%     %     plot(intXSpace + bwxs' + space*t,bwys','-r');
% end
% contourEvolution = imdilate(contourEvolution,strel('square',3));
% % contourEvolution(contourEvolution == 0) = nan;
% 
% yTick = 1:200:ySize;
% yTickLabel = (1:200:ySize)-1;
% h = figure; 
% hold on;
% imagesc(contourEvolution);
% haxes = get(h,'CurrentAxes');
% set(haxes,'XTickLabel',[]);
% set(haxes,'YTick',yTick);
% set(haxes,'YTickLabel',yTickLabel);
% set(haxes,'FontSize',32);
% hold off;
% contourEvolutionFname = [dirs.plithotaxis dirs.expname '_contourEvolution'];
% eval(sprintf('print -dbmp16m  %s', [contourEvolutionFname '.bmp']));
% 
% hold off;
%%

maxTime = ntime * params.timePerFrame;

% xTick = 1:(60/params.timePerFrame):((maxTime/params.timePerFrame)+1);
% xTickLabel = 0:60:maxTime;
xTick = 1:(200/params.timePerFrame):((maxTime/params.timePerFrame)+1);
xTickLabel = 0:200:maxTime;
yTick = 1:200:ySize;
yTickLabel = (1:200:ySize)-1;

% h = figure;
% imagescnan(seeds);
% hold on;
% haxes = get(h,'CurrentAxes');
% set(haxes,'XLim',[1,maxTime/params.timePerFrame]);
% set(haxes,'XTick',xTick);
% set(haxes,'XTickLabel',xTickLabel);
% set(haxes,'YTick',yTick);
% set(haxes,'YTickLabel',yTickLabel);
% set(haxes,'FontSize',32);
% xlabel('Time (minutes)','FontSize',32); ylabel('Y-axis','FontSize',32);
% set(h,'Color','none');
% hold off;
% export_fig_biohpc(plithotaxisFname);

fprintf('Seeds of plithotaxis: starting figures\n');

h = figure;
imagescnan(seedsLocation);
hold on;
haxes = get(h,'CurrentAxes');
set(haxes,'XLim',[1,maxTime/params.timePerFrame]);
set(haxes,'XTick',xTick);
set(haxes,'XTickLabel',xTickLabel);
set(haxes,'YTick',yTick);
set(haxes,'YTickLabel',yTickLabel);
set(haxes,'FontSize',32);
xlabel('Time (minutes)','FontSize',32); ylabel('Y-axis','FontSize',32);
set(h,'Color','none');
hold off;
fprintf('Seeds of plithotaxis: just before saving\n');
drawnow;
export_fig_biohpc([plithotaxisFname '.eps']);
% eval(sprintf('print -dbmp16m  %s', [plithotaxisFname '.bmp']));

fprintf('Seeds of plithotaxis: done figure 1\n');

% export_fig_biohpc(plithotaxisLocationFname);

% ROI
load([dirs.roiData pad(time(1),3) '_roi.mat']); % ROI
ROI0 = ROI;
load([dirs.roiData pad(time(end),3) '_roi.mat']); % ROI
ROI1 = ROI; clear ROI;

fprintf('Seeds of plithotaxis: before saving\n');

save([plithotaxisFname '.mat'],'seeds','seedsLocation','strainEventsOutput'); % 'xSize','ySize','ROI0','ROI1'

fprintf('Seeds of plithotaxis: before outputSeeds\n');

outputSeeds(plithotaxisFname,params,ntime,ySize);

fprintf('Seeds of plithotaxis: after outputSeeds\n');
% xTick = 1:200:xSize;
% xTickLabel = (1:200:xSize)-1;
% 
% h = figure;
% imagesc(ROI0);
% hold on;
% haxes = get(h,'CurrentAxes');
% set(haxes,'XTick',xTick);
% set(haxes,'XTickLabel',xTickLabel);
% set(haxes,'YTick',yTick);
% set(haxes,'YTickLabel',yTickLabel);
% set(haxes,'FontSize',32);
% xlabel('X','FontSize',32); ylabel('Y','FontSize',32);
% set(h,'Color','none');
% hold off;
% eval(sprintf('print -dbmp16m  %s', [dirs.plithotaxis dirs.expname '_roi0.bmp']));
% 
% h = figure;
% imagesc(ROI1);
% hold on;
% haxes = get(h,'CurrentAxes');
% set(haxes,'XTick',xTick);
% set(haxes,'XTickLabel',xTickLabel);
% set(haxes,'YTick',yTick);
% set(haxes,'YTickLabel',yTickLabel);
% set(haxes,'FontSize',32);
% xlabel('X','FontSize',32); ylabel('Y','FontSize',32);
% set(h,'Color','none');
% hold off;
% eval(sprintf('print -dbmp16m  %s', [dirs.plithotaxis dirs.expname '_roi1.bmp']));

%% Speed of leader cells
% allSpeeds = [];
% for t = time
%     load([dirs.mfData pad(t,3) '_mf.mat']); % ROI
%     speed = sqrt(dxs.^2 + dys.^2);
%     %     ys = find(~isnan(seedsLocation(:,t)));
%     %     xs = seedsLocation(ys,t);
%     %     allSpeeds = [allSpeeds speed(ys + (xs - 1)*size(dxs,1))'];
%     
%     load([dirs.roiData pad(t,3) '_roi.mat']); % ROI
%     boundaries = bwboundaries(ROI);    
%     bws = boundaries{1};
%     inds = bws(:,1) > 1 & bws(:,2) > 1 & bws(:,1) < size(ROI,1);
%     bwys = bws(inds,1)';
%     bwxs = bws(inds,2)';
%     allSpeeds = [allSpeeds speed(bwys + (bwxs - 1)*size(dxs,1))];
% end
% 
% meanSpeed = nanmean(allSpeeds);
% 
% seedsSpeed = nan(ySize,ntime);
% for t = time
%     load([dirs.mfData pad(t,3) '_mf.mat']); % ROI
%     speed = sqrt(dxs.^2 + dys.^2);
%     ys = find(~isnan(seedsLocation(:,t)));
%     xs = seedsLocation(ys,t);
%     speeds = speed(ys + (xs - 1)*size(dxs,1));
%     %     meanSpeed = nanmean(speeds);
%     seedsSpeedFlat = seedsSpeed(:);
%     seedsSpeedFlat(ys + (t - 1)*size(seedsSpeed,1)) = speeds./meanSpeed;
%     seedsSpeed = reshape(seedsSpeedFlat,size(seedsSpeed));
% end

%% Video of "events"
% createVideo(strainEventsOutput,dirs,time);

fprintf('Seeds of plithotaxis: finish\n');

end

%% Revision: getStrainEvents moved to a seperate file

% % %%
% % function [output] = getStrainEvents(seeds,seedsLocation,patchSize)
% % 
% % seeds1 = false(ceil(size(seeds,1)/patchSize),size(seeds,2));
% % seedsLocation1 = nan(size(seeds1));
% % for y = 1 : size(seeds1,1)
% %     for x = 1 : size(seeds1,2)
% %         yorig = (y*patchSize-patchSize+1);
% %         seeds1(y,x) = (sum(seeds(yorig:min(yorig+patchSize-1,size(seeds,1)),x)) > length(yorig:min(yorig+patchSize-1,size(seeds,1)))/3); % just arbitrary 1/3
% %         if seeds1(y,x)
% %             seedsLocation1(y,x) = min(seedsLocation(yorig:min(yorig+patchSize-1,size(seeds,1)),x)); % max??
% %         end
% %     end
% % end
% % 
% % [patterns] = getEventPatterns();
% % 
% % [ySize,ntime] = size(seeds1);
% % strainEvents = cell(1,ntime);
% % 
% % ncellsEffected = 0;
% % for t = 1 : ntime-3    
% %     nevents = 0;    
% %     events = {};
% %     
% %     % forward    
% %     for y = 1 : ySize - 3
% %         bb = seeds1(y:y+3,t:t+3);
% %         xs = seedsLocation1(y:y+3,t:t+3);
% %         curEventsFwd = findMatchingPattern(bb,patterns,xs,y,false,seedsLocation1,t);
% %         if ~isempty(curEventsFwd)
% %             ncellsEffected = ncellsEffected + 1;
% %             for e = 1 : length(curEventsFwd)
% %                 nevents = nevents + 1;
% %                 events{nevents} = curEventsFwd{e};
% %             end
% %         end
% %     end
% % 
% % 
% %     % backward
% %     for y = ySize : -1 : 4
% %         bb = seeds1(y:-1:y-3,t:t+3);
% %         %             bb = bb(4:-1:1,:);
% %         xs = seedsLocation1(y:-1:y-3,t:t+3);
% %         %             xs = xs(4:-1:1,:);
% %         curEventsBwd = findMatchingPattern(bb,patterns,xs,y,true,seedsLocation1,t);
% %         if ~isempty(curEventsBwd)
% %             ncellsEffected = ncellsEffected + 1;
% %             for e = 1 : length(curEventsBwd)
% %                 nevents = nevents + 1;
% %                 events{nevents} = curEventsBwd{e};
% %             end
% %         end
% %     end
% %     
% % %     if   ~isempty(curEventsFwd) && ~isempty(curEventsBwd)
% % %         ncellsLeaders = ncellsLeaders + 1;
% % %     end
% %     
% %     strainEvents{t} = events;
% % end
% % 
% % % Go back to original resolution
% % seedsVis = int8(seeds1);
% % strainEventsOrigResolution = strainEvents;
% % for t = 1 : ntime
% %     events = strainEvents{t};
% %     eventsOrigRes = strainEvents{t};
% %     for ee = 1 : length(events)
% %         eventsOrigRes{ee}.ys = round(eventsOrigRes{ee}.ys * patchSize - patchSize + 1 + patchSize/2);
% %         seedsVis(events{ee}.ys(1),t) = 2;
% %     end
% %     strainEventsOrigResolution{t} = eventsOrigRes;
% % end
% % 
% % ncellsLeaders = ncellsEffected - sum(seedsVis(:) == 2);
% % nMotionEvents = sum(seeds1,1);
% % nStrainEvents = sum(seedsVis==2,1);
% % 
% % numEventsYAxis = sum(seedsVis==2,2);
% % 
% % output.strainEvents = strainEvents; % < / \ events with 3 sets of coordinates
% % output.strainEventsOrigResolution = strainEventsOrigResolution; % events at original resolution
% % output.nMotionEvents = nMotionEvents; % # motion events for each time point
% % output.nStrainEvents = nStrainEvents; % # strain events for each time point
% % output.nStrainEventsPercentage = sum(nStrainEvents)/sum(nMotionEvents); % # strain events for each time point
% % output.ncellsLeaders = ncellsLeaders; % numnber of < events
% % output.ncellsEffected = ncellsEffected; % number of cells effected (1 per strain event detected)
% % output.numEventsYAxis = numEventsYAxis; % for each y-position, how many strain events detected at that position
% % output.seedsVis = seedsVis; % 1 - motion event, 2 - strain event
% % output.cumsumVis = cumsum(seedsVis==2,2);
% % 
% % end
% % %%
% % function [patterns] = getEventPatterns()
% % patterns = ...
% % {...
% % logical([[1,0,0,0];[0,1,0,0];[0,0,1,0];[0,0,0,0]]),...
% % logical([[1,0,0,0];[0,1,0,0];[0,0,0,0];[0,0,1,0]]),...
% % logical([[1,0,0,0];[0,1,0,0];[0,0,0,1];[0,0,0,0]]),...
% % logical([[1,0,0,0];[0,1,0,0];[0,0,0,0];[0,0,0,1]]),...
% % logical([[1,0,0,0];[0,0,1,0];[0,0,0,1];[0,0,0,0]]),...
% % logical([[1,0,0,0];[0,0,1,0];[0,0,0,0];[0,0,0,1]]),...
% % logical([[1,0,0,0];[0,0,0,0];[0,1,0,0];[0,0,1,0]]),...
% % logical([[1,0,0,0];[0,0,0,0];[0,1,0,0];[0,0,0,1]]),...
% % logical([[1,0,0,0];[0,0,0,0];[0,0,1,0];[0,0,0,1]])...
% %     };
% % end
% % 
% % %%
% % 
% % function [events] = findMatchingPattern(bb,patterns,xs,y,backward,seedsLocation,t)
% % events = {};
% % nevents = 0;
% % ys = repmat(1:4,4,1)' - 1;
% % 
% % if backward
% %     ys = (-1)*ys;
% % end
% % 
% % for p = 1 : length(patterns)
% %     if sum(sum(bb & patterns{p})) < 3
% %         continue;
% %     end
% %     % validate
% %     [patternYs,patternXs] = find(patterns{p});
% %     candidateXs = zeros(1,3);
% %     for i = 1 : 3
% %         candidateXs(i) = xs(patternYs(i),patternXs(i));
% %     end
% %     if sum(candidateXs == sort(candidateXs)) < 3 % so the leader cells is located before the follower cell in x-axis!
% %         continue;
% %     end
% %     
% %     nevents = nevents + 1;
% %     events{nevents}.xs = candidateXs; % no need?
% %     events{nevents}.ys = y + ys(patternYs);
% %     %     events{nevents}.xs = updateXs(events{nevents}.ys,seedsLocation,t); % go back and find location
% % end
% % end

%%

% %% TODO: go back and find correct X (do we need candidateXs for that?)
% function [newXs] = updateXs(ys,seedsLocation,t)
% newXs = nan(1,3);
% for x = 1 : 3
%     for it = t : -1 : 1
%         if seedsLocation(ys(x),)
%         end
%     end
% end
% end

%% Create video
function [] = createVideo(strainEventsOutput,dirs,time)
close all;

timePerFrame = 5;

videoFname = [dirs.plithotaxis dirs.expname '_plithotaxis.avi'];

aviobj = avifile(videoFname,'fps',3,'compression','None');

for t = time    
    %     load([dirs.roiData pad(t,3) '_roi.mat']); % ROI
    I = imread([dirs.images pad(t,3) '.tif']);
    %     perimWidth = round(max(size(I)) / 200);
    %     I(dilate(bwperim(ROI,8),perimWidth)) = max(255,max(I(:)));
    
    %     h = figure('visible','off'); imagesc(I); colormap(gray);
    h = figure; imagesc(I); colormap(gray);
    hold on;
    
    %     events = strainEventsOutput.strainEventsOrigResolution{t};
    %     for ee = 1 : length(events)
    %         %         plot(events{ee}.xs(1),events{ee}.ys(1),'or','MarkerFaceColor','r','MarkerSize',5);
    %         plot(events{ee}.xs(1),events{ee}.ys(1),'or','LineWidth',2,'MarkerSize',5);
    %     end
    
    plot([1,size(I,2)],[266,266],'-r','LineWidth',3);
    plot([1,size(I,2)],[445,445],'-r','LineWidth',3);

    text(size(I,1)-500,size(I,2)-100,sprintf('%d minutes',round(t*timePerFrame)),'color','w','FontSize',15);
    haxes = get(h,'CurrentAxes');
    set(haxes,'XTick',[]);
    set(haxes,'XTickLabel',[]);
    set(haxes,'YTick',[]);
    set(haxes,'YTickLabel',[]);
    %     drawnow;
    movieFrame = getframe(h);
    aviobj = addframe(aviobj, movieFrame);
    hold off;
    close all;
end
aviobj = close(aviobj);
end

function [] = outputSeeds(plithotaxisFname,params,ntime,ySize)
close all;
load([plithotaxisFname '.mat']);

maxTime = ntime * params.timePerFrame;

% xTick = 1:(60/params.timePerFrame):((maxTime/params.timePerFrame)+1);
% xTickLabel = 0:60:maxTime;
xTick = 1:(200/params.timePerFrame):((maxTime/params.timePerFrame)+1);
xTickLabel = 0:200:maxTime;
yTick = 1:200:ySize;
yTickLabel = (1:200:ySize)-1;

fprintf('outputSeeds: before plotting\n');

h = figure;
imagescnan(seedsLocation);
hold on;

strainEvents = strainEventsOutput.strainEventsOrigResolution;
for t = 1 : ntime
     curEvents = strainEvents{t};
     for e = 1: length(curEvents)
         curY = curEvents{e}.ys(1);
         plot(t,curY,'ok','MarkerFaceColor','k','MarkerSize',4,'LineWidth',2);         
     end
end
    
haxes = get(h,'CurrentAxes');
set(haxes,'XLim',[1,maxTime/params.timePerFrame]);
set(haxes,'XTick',xTick);
set(haxes,'XTickLabel',xTickLabel);
set(haxes,'YTick',yTick);
set(haxes,'YTickLabel',yTickLabel);
set(haxes,'FontSize',32);
xlabel('Time (minutes)','FontSize',32); ylabel('Y-axis','FontSize',32);
set(h,'Color','none');
hold off;
fprintf('outputSeeds: before saving\n');
export_fig_biohpc([plithotaxisFname '_seeds.eps']);
% eval(sprintf('print -dbmp16m  %s', [plithotaxisFname '_seeds.bmp']));
fprintf('outputSeeds: after saving\n');
end
