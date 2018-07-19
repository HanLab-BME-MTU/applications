function [globalMeas] = GCAfindPausingInNeuriteOutgrowthTrajectory(neuriteOutgrowth,varargin)
%GCAfindPausingInNeuriteTrajectory
% This function fits a cubic spline to the neurite trajectory measurements-
% calculates the velocity from this spline and then breaks up the trajectory
% using a user supplied velocity threshold to detect growth cone 'pausing'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%  neuriteOutgrowth:    REQUIRED: nx1 double array: containing the neurite
%                                 length measurements (in um/min)
%                                 where n is the number of time points
%                                 measured
%
%  threshForPause:      PARAM:    scalar: velocity threshold in (um/min), below will be
%                                 considered neurite pausing
%                                 DEFAULT: 0.5 um/min
%
%  secPerFrame          PARAM:    scalar: time interval per frame in sec
%                                 DEFAULT: 5 sec
%
%  makePlot:            PARAM:    logical: flag for sanity plots
%                                 DEFAULT: true
%
%  splineParam:         PARAM:    scalar >=0 and <=1 or empty: smoothing factor for cubic spline
%                                 where 0 is the least-squares straight line
%                                 fit to the data and 1 is the 'natural
%                                 cubic spline interpolant. If empty will
%                                 allow csaps function to calculate
%                                 smoothing parameter
%                                 DEFAULT: 0.01;
%
%  outPath:             PARAM:    character or empty: if is empty do nothing
%                                 else save results (including any figures)
%                                 in directory specfied by outPath
%                                 DEFAULT: [];
%
% OUTPUT:
%  globalMeas:                  stucture with fields
%                                 .outgrowth.grouping : scalar of the grouping
%                                 variable for the time series
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check Input
ip = inputParser;
ip.addRequired('neuriteOutgrowth',@(x) isvector(x));

ip.addParameter('threshPause',0.5,@isscalar); % currently in um/min - default 0.5 um/min
ip.addParameter('secPerFrame',5, @isscalar);
ip.addParameter('makePlot',true,@islogical);

ip.addParameter('splineParam',0.01,@(x) isscalar(x) || isempty(x));
ip.addParameter('outPath' ,[],@(x) ischar(x) || isempty(x));
ip.addParameter('forTitle',[]);

ip.parse(neuriteOutgrowth,varargin{:});

threshPause = ip.Results.threshPause;
makePlot = ip.Results.makePlot;
secPerFrame = ip.Results.secPerFrame;
splineParam = ip.Results.splineParam;
outPath = ip.Results.outPath;

%% START
d= neuriteOutgrowth;

[nTime,~]=size(neuriteOutgrowth);
if d(1)~= 0
    % subtract distance from the first frame
    d = d-d(1);
end

% perform spline filter on the neurite distance
sd_spline= csaps(linspace(1,nTime,nTime),d,splineParam);

sd=ppval(sd_spline,linspace(1,nTime,nTime));

% calculate velocity from the filtered distance using analytic
% expression
sd_spline2= csaps(linspace(1,nTime,nTime),d,splineParam);

sv_spline=sd_spline2;
sv_spline.order=3;
sv_spline.coefs(:,1)=3*sd_spline2.coefs(:,1);
sv_spline.coefs(:,2)=2*sd_spline2.coefs(:,2);
sv_spline.coefs(:,3)=1*sd_spline2.coefs(:,3);
sv=ppval(sv_spline,linspace(1,nTime,nTime));

a =diff(sv);

% convert to to um per min
sv = sv*60/secPerFrame;

% calculate the transitions in the trajectory
% find transitions-
% pausing regions
trans = diff(sv<threshPause & sv>-threshPause);
trans = trans~=0 ;

[maxVel,locMaxVel] = findpeaks(sv);

[minVel,locMinVel] = findpeaks(-sv);

% remove any points that fall under the velocity threshold.
locMaxVel = locMaxVel(maxVel>threshPause | maxVel <-threshPause);
maxVel = maxVel(maxVel>threshPause | maxVel<-threshPause);
% are they part of the pausing if so remove

% find all locMinVel that are not within the pause region
locMinVel = locMinVel(-minVel>threshPause);
minVel = minVel(-minVel>threshPause);
trans = find(trans);

% add back just the growth trans : Keep if decide to keep growth portion
% transG = trans;

%% Get the groups with local maximum, acceleration, deceleration etc
trans = [trans locMaxVel locMinVel];
%trans = [trans locMaxVel];
trans = sort(trans);

% Sometimes there is overlap in the local max and the
% transition point out of pause for instance if have one point just above
% the threshold.
trans = unique(trans);

grp1 = ones(trans(1)-1+1,1);
trans = [ trans length(d)];

% add another transition at max velocity
grpCell = arrayfun(@(x) repmat(x,[trans(x)-trans(x-1),1]),2:length(trans),'uniformoutput',0);

grp = [grp1; vertcat(grpCell{:})];
nGroups = length(unique(grp));

% get average velocity per piece
meanNeuriteVel = arrayfun(@(i) mean(sv(grp==i)),1:nGroups);

% associate with a state "pause", "acc", "dec", "retract"
pauseFrames = find(sv<threshPause & sv>-threshPause);
% growthFrames = find(sv>threshPause);
retractFrames = find(sv<-threshPause);
aFrames = find(a>0 & sv(1:end-1)>threshPause);
dFrames = find(a <0& sv(1:end-1)>threshPause);
% take out the transition frame
aFrames = setdiff(aFrames,locMaxVel);
dFrames = setdiff(dFrames,locMaxVel);
%% create the grouping variables for pause/acc/dec/retraction
%  While remembering the conversion is a bit cumbersome : numbers are a bit easier to work with
%  than characters and there isn't that many states.
%  pause = 1
%  rectract = 2
%  acc = 3
%  dec = 4

% Grouping Variables By Block in time
if ~isempty(outPath)
    % initiate a global param grouping structure
    globalMeas.outgrowth.grouping = grp;
    frames = 1:length(grp); % make frames
    
    gpFrames = arrayfun(@(i) frames(grp==i) ,1:nGroups,'uniformoutput',0);
    
    stateGrp = zeros(numel(gpFrames),1);
    
    aGroups = cellfun(@(x) ~isempty(intersect(x,aFrames)),gpFrames);
    % make sure not in more than one group
    stateGrp(aGroups,1) = 3;
    
    dGroups = cellfun(@(x) ~isempty(intersect(x,dFrames)),gpFrames);
    stateGrp(dGroups,1) = 4;
    
    % Grouping Variable By Frame (repmat of above)
    
    % find those frames that overlap with pausing
    pauseGroups = cellfun(@(x) ~isempty(intersect(x,pauseFrames)), gpFrames)  ;
    stateGrp(pauseGroups,1) = 1;
    
    retractGroups = cellfun(@(x) ~isempty(intersect(x,retractFrames)), gpFrames)  ;
    stateGrp(retractGroups,1) = 2;
    
    statesAllFrames = arrayfun(@(x) repmat(stateGrp(x),length(gpFrames{x}),1),1:numel(stateGrp),'uniformoutput',0);
    statesAllFrames = vertcat(statesAllFrames{:});
    
    %% Store Info
    globalMeas.outgrowth.groupedDistOrig = arrayfun(@(i) d(grp==i) ,1:nGroups,'uniformoutput',0);
    globalMeas.outgrowth.groupedVelSmoothed = arrayfun(@(i) sv(grp==i) ,1:nGroups,'uniformoutput',0);
    globalMeas.outgrowth.groupedFrames = gpFrames;
    % globalMeas.outgrowth.colors = c; % in case want to maintain consistency ...
    %globalMeas.outgrowth.pauseFrames = pauseFrames;
    globalMeas.outgrowth.maxNeuriteVelFrames = locMaxVel;
    globalMeas.outgrowth.stateGrp  = stateGrp;
    globalMeas.outgrowth.stateAllFrames = statesAllFrames;
    
    %% Make Plot if user specifies
    
    % Plot with acc/dec
    if makePlot
        % plot the velocity in red and the transitions
        fsFigure(0.75,'visible','off');
        
        subplot(2,1,2)
        
        % plot original points
        inTime = (1:length(sd)).*secPerFrame;
        inTime = inTime-secPerFrame;
        
        hold on
        c(1) = 'c';
        c(2) = 'b';
        c(3) = 'r';
        c(4) = 'm';
        c(5) = 'k';
        
        forC = globalMeas.outgrowth.stateGrp';
        forC(forC==0) = 5;
        forC = num2cell(forC);
        
        velBinned = globalMeas.outgrowth.groupedVelSmoothed;
        distBinned = globalMeas.outgrowth.groupedDistOrig;
        frameBinned = globalMeas.outgrowth.groupedFrames;
        
        
        cellfun(@(x,y,z) plot(x.*secPerFrame-secPerFrame,y,'LineWidth',2,'color',c(z)),frameBinned,velBinned,forC);
        % plot(inTime,sv,'r');
        hold on
        
        line([0,inTime(end)],[0,0],'color','k')
        arrayfun(@(x) line([trans(x)*secPerFrame-secPerFrame,trans(x)*secPerFrame-secPerFrame],[-4,10],'color','k'),1:length(trans));
        line([0,0],[-4 10],'color','k');
        axis([0 inTime(end) -4 10]);
        xlabel('Time (s)');
        ylabel('Neurite Outgrowth Velocity (um/min)');
        line([0 inTime(end)],[ip.Results.threshPause,ip.Results.threshPause],'color','k','lineStyle','--');
        line([0 inTime(end)],[-ip.Results.threshPause,-ip.Results.threshPause],'color','k','lineStyle','--');
        
        scatter(inTime(locMaxVel),maxVel,10,'k','filled');
        scatter(inTime(locMinVel),-minVel,10,'k','filled');
        % get the average of each piece for plotting the text
        placeInTime = arrayfun(@(i) mean(inTime(grp==i)),1:length(grp));
        
        arrayfun(@(i) text(placeInTime(i),10,num2str(meanNeuriteVel(i),2),'FontSize',7,'FontName','Arial'),1:nGroups);
        
        subplot(2,1,1);
        
        hold on
        ylabel('Neurite Length (um)');
        %arrayfun(@(x) plot(inTime(trans(x):trans(x+1)),sd(trans(x):trans(x+1)),'color',c(x+1,:),'Linewidth',2),1:length(trans(1:end-1)));
        
        hold on
        arrayfun(@(x) line([trans(x)*secPerFrame-secPerFrame,trans(x)*secPerFrame-secPerFrame],[-10,25],'color','k'),1:length(trans));
        % plot the original
        cellfun(@(x,y,z) scatter(x.*secPerFrame-secPerFrame,y,40,c(z),'filled'),frameBinned,distBinned,forC);
        
        plot(inTime,sd,'color','k');
        
        axis([0 inTime(end) -10 25]);
        
        if ~isempty(ip.Results.forTitle)
            title(ip.Results.forTitle,'Color','k','FontName','Arial','FontSize',14);
        end
        
        if ~isempty(outPath)
            saveas(gcf,[outPath filesep 'NeuriteTrajectorySubUnits.fig']);
            saveas(gcf,[outPath filesep 'NeuriteTrajectorySubUnits.png']);
            saveas(gcf,[outPath filesep 'NeuriteTrajectorySubUnits.eps'],'psc2');
        end
    end % if makePlot
    %% just do the same for transG
    %  clear grp1 grp
    %  grp1 = ones(trans(1)-1+1,1);
    %  transG = [ transG length(d)];
    %
    %  grpCell = arrayfun(@(x) repmat(x,[transG(x)-transG(x-1),1]),2:length(transG),'uniformoutput',0);
    %  hold on
    %  grp = [grp1; vertcat(grpCell{:})];
    %  nGroups = length(unique(grp));
    %
    % % get average velocity per piece
    % meanNeuriteVel = arrayfun(@(i) mean(sv(grp==i)),1:nGroups);
    %
    % pauseFrames = find(sv<threshPause & sv>-threshPause);
    % growthFrames = find(sv>threshPause);
    % retractFrames = find(sv<-threshPause);
    %
    %
    % % Grouping Variables By Block in time
    % if ~isempty(outPath)
    %     % initiate a global param grouping structure
    %     globalMeas.outgrowth.groupingG = grp;
    %     frames = 1:length(grp); % make frames
    %
    %     gpFrames = arrayfun(@(i) frames(grp==i) ,1:nGroups,'uniformoutput',0);
    %     %state = cell(length(gpFrames),1);
    %
    %     stateGrp = zeros(numel(gpFrames),1);
    %
    %     % Grouping Variable By Frame (repmat of above)
    %
    %     growthGroups = cellfun(@(x) ~isempty(intersect(x,growthFrames)),gpFrames);
    %     stateGrp(growthGroups,1) = 5; % generic growth is now 5
    %
    %     % find those frames that overlap with pausing
    %     pauseGroups = cellfun(@(x) ~isempty(intersect(x,pauseFrames)), gpFrames)  ;
    %     stateGrp(pauseGroups,1) = 1;
    %
    %     retractGroups = cellfun(@(x) ~isempty(intersect(x,retractFrames)), gpFrames)  ;
    %     stateGrp(retractGroups,1) = 2;
    %
    %     statesAllFrames = arrayfun(@(x) repmat(stateGrp(x),length(gpFrames{x}),1),1:numel(stateGrp),'uniformoutput',0);
    %     statesAllFrames = vertcat(statesAllFrames{:});
    %
    %     %% Store Info
    %     globalMeas.outgrowth.groupedDistOrigG = arrayfun(@(i) d(grp==i) ,1:nGroups,'uniformoutput',0);
    %     globalMeas.outgrowth.groupedVelSmoothedG = arrayfun(@(i) sv(grp==i) ,1:nGroups,'uniformoutput',0);
    %     globalMeas.outgrowth.groupedFramesG = gpFrames;
    %     % globalMeas.outgrowth.colors = c; % in case want to maintain consistency ...
    %     %globalMeas.outgrowth.pauseFrames = pauseFrames;
    %     globalMeas.outgrowth.maxNeuriteVelFramesG = locMaxVel;
    %     globalMeas.outgrowth.stateGrpG  = stateGrp;
    %     globalMeas.outgrowth.stateAllFramesG = statesAllFrames;
    %
    %
    %     %% Make Plot if user specifies
    %
    %     % Plot with acc/dec
    %     if makePlot
    %         % plot the velocity in red and the transitions
    %         fsFigure(0.75,'visible','off');
    %
    %         subplot(2,1,2)
    %
    %         % plot original points
    %         inTime = (1:length(sd)).*secPerFrame;
    %         inTime = inTime-secPerFrame;
    %
    %         hold on
    %         c(1) = 'c';
    %         c(2) = 'b';
    %         c(3) = 'g';
    % %         c(4) = 'm';
    %         c(5) = 'k';
    %
    %         forC = globalMeas.outgrowth.stateGrpG';
    %         forC(forC==0) = 5;
    %         forC = num2cell(forC);
    %
    %         velBinned = globalMeas.outgrowth.groupedVelSmoothedG;
    %         distBinned = globalMeas.outgrowth.groupedDistOrigG;
    %         frameBinned = globalMeas.outgrowth.groupedFramesG;
    %
    %
    %         cellfun(@(x,y,z) plot(x.*ip.Results.secPerFrame-ip.Results.secPerFrame,y,'LineWidth',2,'color',c(z)),frameBinned,velBinned,forC);
    %         % plot(inTime,sv,'r');
    %         hold on
    %
    %         line([0,inTime(end)],[0,0],'color','k')
    %         arrayfun(@(x) line([trans(x)*ip.Results.secPerFrame-ip.Results.secPerFrame,trans(x)*ip.Results.secPerFrame-ip.Results.secPerFrame],[-4,10],'color','k'),1:length(trans));
    %         line([0,0],[-4 10],'color','k');
    %         axis([0 inTime(end) -4 10]);
    %         xlabel('Time (s)');
    %         ylabel('Neurite Elongation Velocity (um/min)');
    %         line([0 inTime(end)],[ip.Results.threshPause,ip.Results.threshPause],'color','k','lineStyle','--');
    %         line([0 inTime(end)],[-ip.Results.threshPause,-ip.Results.threshPause],'color','k','lineStyle','--');
    %
    %         scatter(inTime(locMaxVel),maxVel,10,'k','filled');
    %         scatter(inTime(locMinVel),-minVel,10,'k','filled');
    %         % get the average of each piece for plotting the text
    %         placeInTime = arrayfun(@(i) mean(inTime(grp==i)),1:length(grp));
    %
    %         arrayfun(@(i) text(placeInTime(i),10,num2str(meanNeuriteVel(i),2),'FontSize',7,'FontName','Arial'),1:nGroups);
    %
    %
    %         subplot(2,1,1);
    %
    %
    %         hold on
    %         ylabel('Neurite Length (um)');
    %         %arrayfun(@(x) plot(inTime(trans(x):trans(x+1)),sd(trans(x):trans(x+1)),'color',c(x+1,:),'Linewidth',2),1:length(trans(1:end-1)));
    %
    %         hold on
    %         arrayfun(@(x) line([trans(x)*ip.Results.secPerFrame-ip.Results.secPerFrame,trans(x)*ip.Results.secPerFrame-ip.Results.secPerFrame],[-10,25],'color','k'),1:length(trans));
    %         % plot the original
    %         cellfun(@(x,y,z) scatter(x.*secPerFrame-secPerFrame,y,40,c(z),'filled'),frameBinned,distBinned,forC);
    %
    %         plot(inTime,sd,'color','k');
    %
    %
    %         axis([0 inTime(end) -10 25]);
    %
    %         if ~isempty(ip.Results.forTitle)
    %             title(ip.Results.forTitle,'Color','k','FontName','Arial','FontSize',14);
    %         end
    %
    %
    %         if ~isempty(outPath)
    %
    %             saveas(gcf,[outPath filesep 'NeuriteTrajectorySubUnitsG.fig']);
    %             saveas(gcf,[outPath filesep 'NeuriteTrajectorySubUnitsG.png']);
    %             saveas(gcf,[outPath filesep 'NeuriteTrajectorySubUnitsG.eps'],'psc2');
    %         end
    %
    %     end % if makePlot
    %%
    % add the parameters used
    globalMeas.outgrowth.partitionParams.threshPause = ip.Results.threshPause;
    globalMeas.outgrowth.partitionParams.splineParam = ip.Results.splineParam;
    save([outPath filesep 'globalMeas'],'globalMeas');
end % isempty
end % function