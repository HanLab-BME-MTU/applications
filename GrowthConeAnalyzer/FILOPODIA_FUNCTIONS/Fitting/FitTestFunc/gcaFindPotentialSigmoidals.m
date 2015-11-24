function [slopeMaxNeg,slopeMaxPos] = gcaFindPotentialSigmoidals(yData,varargin)
%testLocalSmoothMax
% Small helper function: fits a cubic spline to the yData 
% calculates the smoothed gradient from this spline 
% and detects the peaks. (something similar is likely buried in the repo) 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT: 
%  yData:    REQUIRED: nx1 double array: containing the neurite 
%                                 length measurements (in um/min)
%                                 where n is the number of time points
%                                 measured 
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
% slopeMaxPos:  sites in the yData that have a local maximum positive slope
% (potential sigmoidal sites)
% 
% slopeMinPos: sites in the yData that have a local minimum positive slope.
% 
%  
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%Check Input
ip = inputParser;
ip.addRequired('yData',@(x) isvector(x));



ip.addParamValue('makePlot',true,@islogical); 

ip.addParamValue('splineParam',0.1,@(x) isscalar(x) || isempty(x)); 
ip.addParamValue('outPath' ,[],@(x) ischar(x) || isempty(x)); 
ip.addParamValue('forTitle',[]); 
ip.addParamValue('backEst',[]); 
ip.addParamValue('mask',[]); % noise mask 
ip.addParamValue('img',[]); % img 
ip.parse(yData,varargin{:});

% threshPause = ip.Results.threshPause;
makePlot = ip.Results.makePlot; 
% secPerFrame = ip.Results.secPerFrame; 
splineParam = ip.Results.splineParam; 
outPath = ip.Results.outPath;

%% START 
d= yData; 

[nTime,~]=size(yData);

% perform spline filter on the yData
sd_spline= csaps(linspace(1,nTime,nTime),d,splineParam);

sd=ppval(sd_spline,linspace(1,nTime,nTime));

% calculate gradient of the yData
sd_spline2= csaps(linspace(1,nTime,nTime),d,splineParam);

sv_spline=sd_spline2;
sv_spline.order=3;
sv_spline.coefs(:,1)=3*sd_spline2.coefs(:,1);
sv_spline.coefs(:,2)=2*sd_spline2.coefs(:,2);
sv_spline.coefs(:,3)=1*sd_spline2.coefs(:,3);
sv=ppval(sv_spline,linspace(1,nTime,nTime));

%% find places where slope is maximally increasing/decreasing 
% these sites will be potential target sites for the fitting. 
% ploblem is techically one of finding multiple significant gaussians then 
% testing which is the first one significantly above the noise 
[~,slopeMaxPos] = findpeaks(sv);
slopeMaxPos = slopeMaxPos(sv(slopeMaxPos)>0); % find increasing portions of the curve

[~,slopeMaxNeg] = findpeaks(-sv); 

%% 
if ~isempty(ip.Results.img)
    nPlots = 3;
else
    nPlots = 2;
end
   %% make plot if user specifies
   if makePlot == true;
       % plot the velocity in red and the transitions
       fsFigure(0.75,'visible','on');
       
       
       subplot(nPlots,1,1)
       
       % plot original points
       %inTime = (1:length(sd)).*secPerFrame;
       inTime = 1:length(sd); 
      % inTime = inTime-secPerFrame;
       
       hold on
       c(1) = 'c'; 
       c(2) = 'b';
       c(3) = 'r';
       c(4) = 'm';
       c(5) = 'k';
       
%        forC = globalParams.outgrowth.stateGrp';
%        forC(forC==0) = 5; 
%        forC = num2cell(forC); 
%        
%        velBinned = globalParams.outgrowth.groupedVelSmoothed;
%        distBinned = globalParams.outgrowth.groupedDistOrig;
%        frameBinned = globalParams.outgrowth.groupedFrames;
       
       plot(inTime,sv,'Linewidth',2,'Color','k'); 
%        cellfun(@(x,y,z) plot(x,y,'LineWidth',2,'color',c(z)),frameBinned,velBinned,forC);
      % plot(inTime,sv,'r');
      hold on 
      
       line([0,length(sd)],[0,0],'color','k')
%        arrayfun(@(x) line([trans(x),trans(x)],[-4,10],'color','k'),1:length(trans));
      % line([0,0],[-4 10],'color','k');
      % axis([0 inTime(end) -4 10]);
%        xlabel('Time (s)');
%        ylabel('Neurite Velocity (um/min)');
%        line([0 600],[0.5,0.5],'color','k','lineStyle','--');
%        line([0 600],[-0.5,-0.5],'color','k','lineStyle','--');
       
      % scatter(inTime(locMaxVel),maxVel,10,'k','filled');
       scatter(inTime(slopeMaxNeg),sv(slopeMaxNeg),50,'r','filled'); 
       % get the average of each piece for plotting the text
%        placeInTime = arrayfun(@(i) mean(inTime(grp==i)),1:length(grp));
       
%        arrayfun(@(i) text(placeInTime(i),10,num2str(meanNeuriteVel(i),2),'FontSize',7,'FontName','Arial'),1:nGroups);
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
       
       subplot(nPlots,1,2);
       
       % assign heat map colors
       
       
       %c = linspecer(length(unique(grp)));
       
      % plot(inTime(1:trans(1)),sd(1:trans(1)),'color',c(1,:),'LineWidth',2);
      
       hold on
       %ylabel('Neurite Outgrowth (um)');
       %arrayfun(@(x) plot(inTime(trans(x):trans(x+1)),sd(trans(x):trans(x+1)),'color',c(x+1,:),'Linewidth',2),1:length(trans(1:end-1)));
       
%        hold on
%        arrayfun(@(x) line([trans(x),trans(x)],[-10,25],'color','k'),1:length(trans));
%        % plot the original
%        cellfun(@(x,y,z) scatter(x,y,40,c(z),'filled'),frameBinned,distBinned,forC); 
       
       plot(inTime,sd,'color','k');
       scatter(inTime,yData,10,'k','filled'); 
          scatter(inTime(slopeMaxNeg),yData(slopeMaxNeg),50,'r','filled'); 
          scatter(inTime(slopeMaxPos),yData(slopeMaxPos),50,'b','filled'); 
       
       %axis([0 inTime(end) -10 25]);
       
       if ~isempty(ip.Results.forTitle)
           title(ip.Results.forTitle,'Color','k','FontName','Arial','FontSize',14);          
       end
       
       if ~isempty(ip.Results.backEst) 
           line([inTime(1),inTime(end)],[ip.Results.backEst,ip.Results.backEst],'Color','k','Linestyle','--'); 
       end 
       
       if nPlots == 3 
           subplot(nPlots,1,3);
           imshow(-img,[]); 
           hold on 
           spy(mask,'r');           
       end 
       
       
       if ~isempty(outPath)
           
           saveas(gcf,[outPath filesep  ip.Results.forTitle{2} '.fig']);
           saveas(gcf,[outPath filesep  ip.Results.forTitle{2} '.png']);
           %saveas(gcf,[outPath filesep 'NeuriteTrajectorySubUnits.eps'],'psc2'); 
       end
       
      
       
   end % if makePlot
   



   
   
% save([outPath filesep 'globalParams'],'globalParams'); 
end 
%end

