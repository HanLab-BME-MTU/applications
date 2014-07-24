function [glbMaxPosMean,glbMaxValMean,glbMaxPosSEM,glbMaxValSEM,glbMaxPosCI95,glbMaxValCI95,glbMaxPosList,glbMaxValList]=perfJackKnife(corrSets,maxLag,normVar)
% This is rather a bootstrap

% skip always one data set
numSets=length(corrSets);
wholeList=1:numSets;

for NAindex=1:100 %wholeList
    clear edgF;
    
%     % for JackKnife:
%     idxList=setdiff(wholeList,NAindex);
%     
%     % for JackKnife with n>1 (random resampling, not all possible permutations!):
%     % take out only e.g. 0.66% of the data:
%     cutOff=max(1,floor(0.66*numSets));
%     perm=randperm(numSets);
%     idxList=sort(perm(1:cutOff));
    
    % for Bootstrap (random resampling):
    idxList=randi(numSets,1,numSets);
     
    %**************************************************************************
    % 3) Extract the time course of forces and intensities                    *
    %**************************************************************************
    % run through all clusters and collect the data.
    inIdx=1;
    for idx=idxList
        edgF(inIdx).observations = corrSets(idx).fcorr;
        edgI(inIdx).observations = corrSets(idx).Icorr;
        inIdx=inIdx+1;
    end
    
    % Calculate the autocorrelation:
    out1=crossCorr(edgF,edgF,maxLag,normVar); % in the ideal case these entries are all -1!
    cFF=out1(:,1);
    cFF_std=out1(:,2);
    
    out1=crossCorr(edgI,edgI,maxLag,normVar); % in the ideal case these entries are all -1!
    cII=out1(:,1);
    cII_std=out1(:,2);
    
    % Calculate the cross correlation:
    out1=crossCorr(edgF,edgI,maxLag,normVar); % in the ideal case these entries are all -1!
    cFI=out1(:,1);
    cFI_std=out1(:,2);
    
    %display(['cFI(dt=0): ',num2str(cFI(maxLag+1)),'+-(',num2str(cFI_std(maxLag+1)),')']);
    
    % generate the smoothing spline:
    x=-maxLag:1:maxLag;
    y=cFI;
    p=0.2;%0.2; %take 0.2 for 8kPa and 0.5 for 35kPa, the 35kPa is more spiked and therefore should be smoothed less.
    w=min(cFI_std.^2)./(cFI_std.^2); % in this way 0<=w<=1
    % generate the spline form:
    sp=csaps(x,y,p,[],w);
    
    % calculate the first derivative:
    sp_p    = fnder(sp,1);
    sp_pp   = fnder(sp,2);
    
    % find the extrema:
    zeroVals = fnzeros(sp_p);
    zeroVals=zeroVals(1,:);
    
    % check that second derivative is >0
    zeroVals_pp=fnval(sp_pp,zeroVals);
    
    % the max pos are:
    maxPos=zeroVals(zeroVals_pp<=0);
    
    glbMaxVal=NaN;
    glbMaxPos=NaN;
    if ~isempty(maxPos)
        maxVals=fnval(sp,maxPos);
        % the global maximum:
        [glbMaxVal,id]=max(maxVals);
        glbMaxPos=maxPos(id);
    end
    glbMaxPosList(NAindex)=glbMaxPos;
    glbMaxValList(NAindex)=glbMaxVal;
    
    
    figure(111)
    title('The cross correlation for cFI')
    errorbar(-maxLag:1:maxLag,cFI,cFI_std,'r')
    hold on
    fnplt(sp);
    plot(-maxLag:1:maxLag,cFI,'k')
    % ylim([-1 1])
    xlim([-maxLag maxLag])
    title('F / I')
    xlabel('dframes')
    ylabel('corr')
    hold off
    
    display(['run: ',num2str(NAindex),'cFI: max x-corr val: ',num2str(glbMaxVal), ' at time-lag [frames]: ',num2str(glbMaxPos)]);
end

% Check if there are some unreasonable values:
badMax=[];
badMax=(glbMaxPosList<-5 | glbMaxPosList>5 | isnan(glbMaxPosList));
if ~isempty(badMax)
    display(['rejected ',num2str(sum(badMax)),' bad maxima: ',num2str(glbMaxPosList(badMax))])
end

glbMaxPosList=glbMaxPosList(~badMax);
glbMaxValList=glbMaxValList(~badMax);
[h,p] = ttest(glbMaxPosList);

glbMaxPosMean  = mean(glbMaxPosList);
glbMaxPosSEM   = std(glbMaxPosList); % this is the bootstrap estimate for the standard error of the mean given by the std of the distributions of the means in the bootstrap samples.
glbMaxPosCI95  = prctile(glbMaxPosList,[2.5 97.5]);

glbMaxValMean  = mean(glbMaxValList);
glbMaxValSEM   = std(glbMaxValList);
glbMaxValCI95  = prctile(glbMaxValList,[2.5 97.5]);

display(['value:    mean+-SEM (95% conf. intv.)= ',num2str(glbMaxValMean),' +- ',num2str(glbMaxValSEM),' (',num2str(glbMaxValCI95),')']);

display(['position: mean+-SEM (95% conf. intv.)= ',num2str(glbMaxPosMean),' +- ',num2str(glbMaxPosSEM),' (',num2str(glbMaxPosCI95),')']);
display(['p-value t-test = ',num2str(p)]);





% for these parameters
% normVar=1;
% tBtwFrms=150;
% aveType='none'; % first checks 'none', 'nanmean', 'mean' makes little difference
% maxLag =round(7200/tBtwFrms); % round(7200/tBtwFrms) means a maxLag of 2h

% I got from:
% random drawing of numSet datasets (Bootstrap), 250 times:
% all data (none rejected):
% mean+-95%-error= -0.0038131 +- 0.043399
% p-value t-test = 0.86349
%----------------------------------------
% assembling data (none rejected): 
% mean+-95%-error= -0.96007 +- 0.066015
% p-value t-test = 2.5374e-80
%----------------------------------------
% disassembly data (12 rejected): 
% mean+-95%-error= 0.32311 +- 0.14592
% p-value t-test = 2.1273e-05


% random drawing of 50% of the data (ohne zuruecklegen):
% all data:
% mean+-95%-error= 0.05301 +- 0.075215
% p-value t-test = 0.16907
%----------------------------------------
% assembling data: 
% mean+-95%-error= -1.0178 +- 0.26991
% p-value t-test = 3.815e-07
%----------------------------------------
% disassembly data: 
% mean+-95%-error= 0.6541 +- 0.32602
% p-value t-test = 0.0025422
%----------------------------------------
