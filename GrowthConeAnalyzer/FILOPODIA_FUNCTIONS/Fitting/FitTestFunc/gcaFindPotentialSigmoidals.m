function [slopeMaxNeg,slopeMaxPos,valuesNeg,valuesPos,sd,sv] = gcaFindPotentialSigmoidals(yData,varargin)
% gcaFindPotentialSigmoidals
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

ip.addParamValue('makePlot',false,@islogical); 

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
% problem is techically one of finding multiple significant gaussians then 
% testing which is the first one significantly above the noise 
[~,slopeMaxPos] = findpeaks(sv);
slopeMaxPos = slopeMaxPos(sv(slopeMaxPos)>0); % find increasing portions of the curve

[~,slopeMaxNeg] = findpeaks(-sv); 
slopeMaxNeg = slopeMaxNeg(sv(slopeMaxNeg)<0); % find decreasing portions of the curve

valuesNeg = sv(slopeMaxNeg); 
valuesPos = sv(slopeMaxPos); 
%% 
if ~isempty(ip.Results.img)
    nPlots = 3;
else
    nPlots = 2;
end
%% make plot if user specifies
if makePlot == true;
    % plot the velocity in red and the transitions
   hFig = fsFigure(0.75,'visible','off');
    
    %% Plot the original 
    subplot(nPlots,1,1);
    
    hold on
    inTime = 1:length(sd); 
    plot(inTime,sd,'color','k');
    h1 = scatter(inTime,yData,10,'k','filled');
     scatter(inTime(slopeMaxNeg),yData(slopeMaxNeg),50,'r','filled');
     scatter(inTime(slopeMaxPos),yData(slopeMaxPos),50,'g','filled');
  
    if ~isempty(ip.Results.forTitle)
        title(strrep(ip.Results.forTitle,'_',' '),'Color','k','FontName','Arial','FontSize',12);
    end
    
    if ~isempty(ip.Results.backEst)
       h2 =  line([inTime(1),inTime(end)],[ip.Results.backEst,ip.Results.backEst],'Color','k','Linestyle','--');
    end
    
   
    ylabel('Fluorescence Intensity (AU)','FontSize',12,'FontName','Arial');
    
    
    if nPlots == 3
        subplot(nPlots,1,3);
        imshow(-img,[]);
        hold on
        spy(mask,'r');
    end
    
    %% Plot the smoothed gradient 
    subplot(nPlots,1,2)
    
    % plot original points
    %inTime = (1:length(sd)).*secPerFrame;
    inTime = 1:length(sd);
    
    plot(inTime,sv,'Linewidth',2,'Color','k');
    
    hold on
    
    % line at zero 
    line([0,length(sd)],[0,0],'color','k')

    % red is potential minimum (maximally decreasing)
    h1 = scatter(inTime(slopeMaxNeg),sv(slopeMaxNeg),50,'r','filled');
    
    % green is potential positive 
    h2 = scatter(inTime(slopeMaxPos),sv(slopeMaxPos),50,'g','filled');
    
    ylabel({'Smoothed Derivative' ; ['SplineParam= ' num2str(splineParam)]},'FontSize',12,'FontName','Arial'); 
    xlabel({'Distance Along Filopodia' ; '(Pixel Number Before Convert)'},'FontSize',12,'FontName','Arial'); 
    box off
    
    
    if ip.Results.makePlot
        saveas(gcf,[outPath filesep  ip.Results.forTitle '.fig']);
        saveas(gcf,[outPath filesep  ip.Results.forTitle '.png']);
    end
       
end % if makePlot

end 


