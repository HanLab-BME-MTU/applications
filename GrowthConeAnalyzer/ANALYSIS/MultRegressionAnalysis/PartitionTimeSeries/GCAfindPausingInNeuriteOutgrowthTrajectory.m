function [globalParams] = GCAfindPausingInNeuriteOutgrowthTrajectory(neuriteOutgrowth,varargin)
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
%  globalParams:                  stucture with fields
%                                 .outgrowth.grouping : scalar of the grouping
%                                 variable for the time series 
%  
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%Check Input
ip = inputParser;
ip.addRequired('neuriteOutgrowth',@(x) isvector(x));

ip.addParamValue('threshPause',0.5,@isscalar); % currently in um/min - default 0.5 um/min
ip.addParamValue('secPerFrame',5, @isscalar); 
ip.addParamValue('makePlot',true,@islogical); 

ip.addParamValue('splineParam',0.01,@(x) isscalar(x) || isempty(x)); 
ip.addParamValue('outPath' ,[],@(x) ischar(x) || isempty(x)); 
ip.addParamValue('forTitle',[]); 
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
sd_spline= csaps(linspace(1,nTime,nTime),d,0.01);

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
sv = sv*60/secPerFrame; % convert to um per frame*60 

% calculate the transitions in the trajectory 
% find transitions 
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
 trans = [trans locMaxVel locMinVel]; 
 %trans = [trans locMaxVel]; 
 trans = sort(trans); 

grp1 = ones(trans(1)-1+1,1); 
trans = [ trans length(d)]; 

% add another transition at max velocity 


grpCell = arrayfun(@(x) repmat(x,[trans(x)-trans(x-1),1]),2:length(trans),'uniformoutput',0);
hold on 
grp = [grp1; vertcat(grpCell{:})]; 
nGroups = length(unique(grp)); 

% get average velocity per piece 
meanNeuriteVel = arrayfun(@(i) mean(sv(grp==i)),1:nGroups);


% associate with a state "pause", "acc", "dec", "retract" 

pauseFrames = find(sv<threshPause & sv>-threshPause);
growthFrames = find(sv>threshPause); 
retractFrames = find(sv<-threshPause);
aFrames = find(a>0 & sv(1:end-1)>threshPause);  
dFrames = find(a <0& sv(1:end-1)>threshPause); 
% take out the transition frame 
aFrames = setdiff(aFrames,locMaxVel);
dFrames = setdiff(dFrames,locMaxVel); 
%% create the grouping variables for pause/acc/dec/retraction 
% I know this is a bit cumbersome but numbers are a bit easier to work with
% than characters and there isn't that many states. 
%  pause = 1
%  rectract = 2
%  acc = 3
%  dec = 4

% Grouping Variables By Block in time 
if ~isempty(outPath)
    % initiate a global param grouping structure 
    globalParams.outgrowth.grouping = grp;
    frames = 1:length(grp); % make frames
    
    gpFrames = arrayfun(@(i) frames(grp==i) ,1:nGroups,'uniformoutput',0);
    %state = cell(length(gpFrames),1);
    
    stateGrp = zeros(numel(gpFrames),1); 
    
    
   
    
    aGroups = cellfun(@(x) ~isempty(intersect(x,aFrames)),gpFrames);
    % make sure not in more than one group 
    stateGrp(aGroups,1) = 3; 
%     for i = 1:length(aGroups)
%         state{aGroups(i),1} = 'acc';
%     end
   
    dGroups = cellfun(@(x) ~isempty(intersect(x,dFrames)),gpFrames);
    stateGrp(dGroups,1) = 4; 
    
%     for i = 1:length(dGroups)
%         state{dGroups(i),1} = 'dec';
%     end

% Grouping Variable By Frame (repmat of above) 

 % find those frames that overlap with pausing
    pauseGroups = cellfun(@(x) ~isempty(intersect(x,pauseFrames)), gpFrames)  ;
    stateGrp(pauseGroups,1) = 1; 
%     for i = 1:length(pauseGroups)
%         state{pauseGroups(i),1} = 'pause';
%     end
    
    retractGroups = cellfun(@(x) ~isempty(intersect(x,retractFrames)), gpFrames)  ;
    stateGrp(retractGroups,1) = 2; 
%     for i = 1:length(retractGroups)
%         state{retractGroups(i),1} = 'retract';
%     end
    %acGroups =



   statesAllFrames = arrayfun(@(x) repmat(stateGrp(x),length(gpFrames{x}),1),1:numel(stateGrp),'uniformoutput',0); 
   statesAllFrames = vertcat(statesAllFrames{:}); 

%% Store Info 
    globalParams.outgrowth.groupedDistOrig = arrayfun(@(i) d(grp==i) ,1:nGroups,'uniformoutput',0);  
    globalParams.outgrowth.groupedVelSmoothed = arrayfun(@(i) sv(grp==i) ,1:nGroups,'uniformoutput',0); 
    globalParams.outgrowth.groupedFrames = gpFrames; 
    % globalParams.outgrowth.colors = c; % in case want to maintain consistency ... 
    %globalParams.outgrowth.pauseFrames = pauseFrames; 
    globalParams.outgrowth.maxNeuriteVelFrames = locMaxVel; 
    globalParams.outgrowth.stateGrp  = stateGrp; 
    globalParams.outgrowth.stateAllFrames = statesAllFrames;    
   
  
   %% make plot if user specifies
   if makePlot == true;
       % plot the velocity in red and the transitions
       fsFigure(0.75,'visible','off');
       
       
       subplot(2,1,1)
       
       % plot original points
       inTime = (1:length(sd)).*secPerFrame;
       inTime = inTime-secPerFrame;
       
       hold on
       c(1) = 'c'; 
       c(2) = 'b';
       c(3) = 'r';
       c(4) = 'm';
       c(5) = 'k';
       
       forC = globalParams.outgrowth.stateGrp';
       forC(forC==0) = 5; 
       forC = num2cell(forC); 
       
       velBinned = globalParams.outgrowth.groupedVelSmoothed;
       distBinned = globalParams.outgrowth.groupedDistOrig;
       frameBinned = globalParams.outgrowth.groupedFrames;
       
       
       cellfun(@(x,y,z) plot(x.*5-5,y,'LineWidth',2,'color',c(z)),frameBinned,velBinned,forC);
      % plot(inTime,sv,'r');
      hold on 
      
       line([0,inTime(end)],[0,0],'color','k')
       arrayfun(@(x) line([trans(x)*5-5,trans(x)*5-5],[-4,10],'color','k'),1:length(trans));
       line([0,0],[-4 10],'color','k');
       axis([0 inTime(end) -4 10]);
       xlabel('Time (s)');
       ylabel('Neurite Velocity (um/min)');
       line([0 600],[0.5,0.5],'color','k','lineStyle','--');
       line([0 600],[-0.5,-0.5],'color','k','lineStyle','--');
       
       scatter(inTime(locMaxVel),maxVel,10,'k','filled');
       scatter(inTime(locMinVel),-minVel,10,'k','filled'); 
       % get the average of each piece for plotting the text
       placeInTime = arrayfun(@(i) mean(inTime(grp==i)),1:length(grp));
       
       arrayfun(@(i) text(placeInTime(i),10,num2str(meanNeuriteVel(i),2),'FontSize',7,'FontName','Arial'),1:nGroups);
       % Assign heat map colors
       %  cMapLength=128; cMap=jet(cMapLength);
       % %
       % %  %data = dwell(idxPer);
       %   dataRange = 2;
       %   meanNeuriteVel(meanNeuriteVel>dataRange)=dataRange;
       %                          mapper=linspace(0,dataRange,cMapLength)';
       % %
       % %                         % get closest colormap index for each feature
       %                          D=createDistanceMatrix(meanNeuriteVel,mapper);
       %                          [sD,idxCMap]=sort(abs(D),2);
       
       %  for k = 1:cMapLength
       %      plot(inTime(trans(x):trans(x+1)),sd(trans(x):trans(x+1)))
       %  end
       
       subplot(2,1,2);
       
       % assign heat map colors
       
       
       %c = linspecer(length(unique(grp)));
       
      % plot(inTime(1:trans(1)),sd(1:trans(1)),'color',c(1,:),'LineWidth',2);
      
       hold on
       ylabel('Neurite Outgrowth (um)');
       %arrayfun(@(x) plot(inTime(trans(x):trans(x+1)),sd(trans(x):trans(x+1)),'color',c(x+1,:),'Linewidth',2),1:length(trans(1:end-1)));
       
       hold on
       arrayfun(@(x) line([trans(x)*5-5,trans(x)*5-5],[-10,25],'color','k'),1:length(trans));
       % plot the original
       cellfun(@(x,y,z) scatter(x.*5-5,y,40,c(z),'filled'),frameBinned,distBinned,forC); 
       
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
   



   
   
save([outPath filesep 'globalParams'],'globalParams'); 
end 
end

