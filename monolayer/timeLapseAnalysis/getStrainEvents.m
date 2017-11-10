function [output] = getStrainEvents(seeds,seedsLocation,patchSize)

seeds1 = false(ceil(size(seeds,1)/patchSize),size(seeds,2));
seedsLocation1 = nan(size(seeds1));
for y = 1 : size(seeds1,1)
    for x = 1 : size(seeds1,2)
        yorig = (y*patchSize-patchSize+1);
        seeds1(y,x) = (sum(seeds(yorig:min(yorig+patchSize-1,size(seeds,1)),x)) > length(yorig:min(yorig+patchSize-1,size(seeds,1)))/3); % just arbitrary 1/3
        if seeds1(y,x)
            seedsLocation1(y,x) = min(seedsLocation(yorig:min(yorig+patchSize-1,size(seeds,1)),x)); % max??
        end
    end
end

seedsLocation2 = fillSeedsLocation(seedsLocation1);

[patterns] = getEventPatterns();

[ySize,ntime] = size(seeds1);
strainEvents = cell(1,ntime);

ncellsEffected = 0;
for t = 1 : ntime-3    
    nevents = 0;    
    events = {};
    
    % forward    
    for y = 1 : ySize - 3
        bb = seeds1(y:y+3,t:t+3);
        xs = seedsLocation2(y:y+3,t:t+3);
        curEventsFwd = findMatchingPattern(bb,patterns,xs,y,false,seedsLocation2,t);
        if ~isempty(curEventsFwd)
            ncellsEffected = ncellsEffected + 1;
            for e = 1 : length(curEventsFwd)
                nevents = nevents + 1;
                events{nevents} = curEventsFwd{e};
            end
        end
    end


    % backward
    for y = ySize : -1 : 4
        bb = seeds1(y:-1:y-3,t:t+3);
        %             bb = bb(4:-1:1,:);
        xs = seedsLocation2(y:-1:y-3,t:t+3);
        %             xs = xs(4:-1:1,:);
        curEventsBwd = findMatchingPattern(bb,patterns,xs,y,true,seedsLocation2,t);
        if ~isempty(curEventsBwd)
            ncellsEffected = ncellsEffected + 1;
            for e = 1 : length(curEventsBwd)
                nevents = nevents + 1;
                events{nevents} = curEventsBwd{e};
            end
        end
    end
    
%     if   ~isempty(curEventsFwd) && ~isempty(curEventsBwd)
%         ncellsLeaders = ncellsLeaders + 1;
%     end
    
    strainEvents{t} = events;
end

% Revision: filter consecutive reoccuring stretching events in same Y
% location in 3 time points ahead
for t = 2 : ntime
    sEvents0 = strainEvents{t-1};
    sEvents1 = strainEvents{t};
    % get y positions, if ys from 1 appear in 0 than exclude from 1
    n0 = length(sEvents0);
    ys0 = nan(1,n0);
    for e = 1 : n0
        ys0(e) = sEvents0{e}.ys(1);
    end
    n1 = length(sEvents1);
    ys1 = nan(1,n1);
    for e = 1 : n1
        ys1(e) = sEvents1{e}.ys(1);
    end
    
    indsToKeep = find(~ismember(ys1,ys0));
    strainEvents{t} = sEvents1(indsToKeep);
end

% Go back to original resolution
seedsVis = int8(seeds1);
strainEventsOrigResolution = strainEvents;
for t = 1 : ntime
    events = strainEvents{t};
    eventsOrigRes = strainEvents{t};
    for ee = 1 : length(events)
        eventsOrigRes{ee}.ys = round(eventsOrigRes{ee}.ys * patchSize - patchSize + 1 + patchSize/2);
        seedsVis(events{ee}.ys(1),t) = 2;
    end
    strainEventsOrigResolution{t} = eventsOrigRes;
end

ncellsLeaders = ncellsEffected - sum(seedsVis(:) == 2);
nMotionEvents = sum(seeds1,1);
nStrainEvents = sum(seedsVis==2,1);

numEventsYAxis = sum(seedsVis==2,2);

output.strainEvents = strainEvents; % < / \ events with 3 sets of coordinates
output.strainEventsOrigResolution = strainEventsOrigResolution; % events at original resolution
output.nMotionEvents = nMotionEvents; % # motion events for each time point
output.nStrainEvents = nStrainEvents; % # strain events for each time point
output.nStrainEventsPercentage = sum(nStrainEvents)/sum(nMotionEvents); % # strain events for each time point
output.ncellsLeaders = ncellsLeaders; % numnber of < events
output.ncellsEffected = ncellsEffected; % number of cells effected (1 per strain event detected)
output.numEventsYAxis = numEventsYAxis; % for each y-position, how many strain events detected at that position
output.seedsVis = seedsVis; % 1 - motion event, 2 - strain event
output.cumsumVis = cumsum(seedsVis==2,2);
end

%%
function [patterns] = getEventPatterns()
patterns = ...
{...
logical([[1,0,0,0];[0,1,0,0];[0,0,1,0];[0,0,0,0]]),...
logical([[1,0,0,0];[0,1,0,0];[0,0,0,0];[0,0,1,0]]),...
logical([[1,0,0,0];[0,1,0,0];[0,0,0,1];[0,0,0,0]]),...
logical([[1,0,0,0];[0,1,0,0];[0,0,0,0];[0,0,0,1]]),...
logical([[1,0,0,0];[0,0,1,0];[0,0,0,1];[0,0,0,0]]),...
logical([[1,0,0,0];[0,0,1,0];[0,0,0,0];[0,0,0,1]]),...
logical([[1,0,0,0];[0,0,0,0];[0,1,0,0];[0,0,1,0]]),...
logical([[1,0,0,0];[0,0,0,0];[0,1,0,0];[0,0,0,1]]),...
logical([[1,0,0,0];[0,0,0,0];[0,0,1,0];[0,0,0,1]])...
    };
end

%%

function [events] = findMatchingPattern(bb,patterns,xs,y,backward,seedsLocation,t)
events = {};
nevents = 0;
ys = repmat(1:4,4,1)' - 1;

if backward
    ys = (-1)*ys;
end

for p = 1 : length(patterns)
    if sum(sum(bb & patterns{p})) < 3
        continue;
    end
    % validate
    [patternYs,patternXs] = find(patterns{p});
    candidateXs = zeros(1,3);
    for i = 1 : 3
        candidateXs(i) = xs(patternYs(i),patternXs(i));
    end
    
    %% Revision: no acceptence of protrusion pattern [1 1 x x]    
    if sum(bb(1:2,1)) == 2
        continue;
    end
    
    %% Revision: first cell is leading the way for sure
    %     if sum(candidateXs == sort(candidateXs)) < 3 % so the leader cells is located before the follower cell in x-axis!
    %         continue;
    %     end
    
    if xs(1,1) < max(xs(:,1)) % so the leader cells is located ahead (not behind) the follower cell in x-axis!
        continue;
    end
    
    nevents = nevents + 1;
    events{nevents}.xs = candidateXs; % no need?
    events{nevents}.ys = y + ys(patternYs);
    %     events{nevents}.xs = updateXs(events{nevents}.ys,seedsLocation,t); % go back and find location
end
end


function [seedsLocation2] = fillSeedsLocation(seedsLocation1)
seedsLocation2 = seedsLocation1;
for i = 2:size(seedsLocation1,1)
    seedsLocation2(:,i) = max(seedsLocation2(:,i),seedsLocation2(:,i-1));
end
end

% function [xsFull] = fillXs(xs)
% xsFull = xs;
% for i = 2:size(xs,1)
%     xsFull(:,i) = max(xs(:,i),xs(:,i-1));
% end
% end

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
%%
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
eval(sprintf('print -dbmp16m  %s', [plithotaxisFname '_seeds.bmp']));
end
