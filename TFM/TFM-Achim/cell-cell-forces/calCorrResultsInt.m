function [corrResults]=calCorrResultsInt(corrSets,maxLag,opt1,normVar,tBtwFrms,aveType,opt2)
% INPUT 
% opt1:  'usefn': take only network forces into account. At the moment, flag
%                is only used to select the correct fm values, one could also use it to
%                select the correct resF values. At the moment, all resF
%                values are taken to get a more reliable value for the resF
%                variance.
%       'usefc': take only cluster forces into account.
%       'usefm': take network forces when possible, cluster forces
%                otherwise.
%    'tBtwFrms': all frame rates will be transformed to this user defined
%                time interval between frames. Of course, this value has to
%                be set very carefully. If it is not given it is assumed
%                that all data have been acquired at the SAME frame rate.
% aveType:       In case tBtwFrms is not empty the user has a few opt1ions
%                to choose how the time points are picked and bined into
%                the bins set by tBtwFrms.
%     'none'   : The default: If there are multiple measurments in the same 
%                bin, then ONLY the last time point in this bin is taken.
%     'mean'   : If there are multiple measurments in the same bin, then
%                the average over all measurements in this bin is taken. If
%                the number of interfaces changes in this time bin, then,
%                the cell is considered to have the minimal connectivity it
%                ever assums during this time interval. This is because, 
%                the simple mean, at an interface that at least once assumes
%                a NaN during the time intervall, will yield a NaN. So
%                information is lost. In such a situation, the finite
%                values at the stable edges will be biased by all values
%                although the cell might change its connectivity. The sum
%                over Fi_tot is underestimated in this situation.
%     'nanmean': If there are multiple measurments in the same bin, then
%                the average over all measurements in this bin is taken. If
%                the number of interfaces changes in this time bin, then,
%                the cell is considered to have the maximal connectivity it
%                ever assums during this time interval. The sum
%                over Fi_tot is overstimated in this situation.
%                COMMENT: I slightly prefer the 'nanmean'.
%**************************************************************************
% 1) Find all frames for this cell where the network has not changed      *
%**************************************************************************

% previously I have checked:
% 1) the frame has been analyzed.
% 2) the cell has not been lost.
% 3) the node still exists.
% 4) the degree = the number of edges is the same.
% 5) all edges are identical.

% 1)-4) can be accounted for findCells. 
% 5) could be checked in here if needed, but not sure if necessary at all.

% For 'usefn', 'deg'=2, 'myo'=1, 'errs'=0 this function will yield the same
% results as the old TFM_part_8_corrAnalysis for cell 4 at frame 20 of:
% /orchestra/groups/lccb-mechanics/analysis/Rosa/TFM/clusters/2010_08_12_TFM_8kPa_5p_05p_10AshMYOIIA_hp94/wellA4_32

doPlotTraces=0;

% calculate the correlation coefficients:
if nargin<2 || isempty(maxLag)
    maxLag=5;
end

if nargin<3 || isempty(opt1)
    opt1='usefm';
end

if nargin<4 || isempty(normVar)
    normVar=0;
end

if nargin<6 || isempty(aveType)
    aveType='none';
end

if nargin<7 || isempty(opt2)
    opt2='useItot';
end

%**************************************************************************
% 1) Pull the right data in new fields of the data structure. 
%**************************************************************************
%Once when this is done, we can simply rely on this field to calculate the
% correlations:
corrSets=pullCorrVals(corrSets,opt1,opt2);


%**************************************************************************
% 2) adapt the frame rates:    
%**************************************************************************
if nargin>=5 && ~isempty(tBtwFrms)
    corrSets=adaptFrameRates(corrSets,tBtwFrms,aveType);
end


%**************************************************************************
% 3) Extract the time course of forces and intensities                    *
%**************************************************************************
% run through all clusters and collect the data.
for idx=1:length(corrSets)
    edgF(idx).observations = corrSets(idx).fcorr;
    edgI(idx).observations = corrSets(idx).Icorr;
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

if doPlotTraces
    for k=1:length(edgF)
        figure()
        plot(1:length(edgF(k).observations),edgF(k).observations,'b');
        hold on;
        plot(1:length(edgI(k).observations),-edgI(k).observations,'r');
        hold off;
        set(gca,'fontsize',16)
    end
end

display(['Avg. XXXX correlation cFI: ',num2str(cFI(maxLag+1)),'+-(',num2str(cFI_std(maxLag+1)),')']);

figure()
plotmatrix([vertcat(edgF.observations) vertcat(edgI.observations)])
title('F / I')


figure()
title('The autco-rrelation for cFF')
plot(-maxLag:1:maxLag,cFF)
ylim([-1 1])
xlim([-maxLag maxLag])
title('Autocorrelation F F')
xlabel('dframes')
ylabel('corr')


figure()
title('The auto-correlation for cII')
plot(-maxLag:1:maxLag,cII)
ylim([-1 1])
xlim([-maxLag maxLag])
title('Autocorrelation I I')
xlabel('dframes')
ylabel('corr')


% generate the smoothing spline:
x=-maxLag:1:maxLag;
y=cFI;
p=0.2;
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

glbMaxVal=[];
glbMaxPos=[];
if ~isempty(maxPos)
    maxVals=fnval(sp,maxPos);
    % the global maximum:
    [glbMaxVal,id]=max(maxVals);
    glbMaxPos=maxPos(id);
end
glbMaxPos
glbMaxVal
    

figure()
title('The cross correlation for cFI')
errorbar(-maxLag:1:maxLag,cFI,cFI_std,'k')
hold on
fnplt(sp,'r');
plot(-maxLag:1:maxLag,cFI,'k')
%xlim([-maxLag maxLag])
xlim([-15 15])
ylim([ 0 0.5])
title('F / I')
xlabel('dframes')
ylabel('corr')
box on
set(gca,'LineWidth',2,'FontSize',20)
hold off

% make the bootstrap analysis:
[glbMaxPosMean,glbMaxValMean,glbMaxPosSEM,glbMaxValSEM,glbMaxPosCI95,glbMaxValCI95,glbMaxPosList,glbMaxValList]=perfJackKnife(corrSets,maxLag,normVar);


% plot the inset:
xmin=-4;
xmax= 4;
ymin= 0.2;
ymax= 0.55;
figure()
title('The cross correlation for cFI zoom up')
% have to add a second fake point in the limbo to get the nice errorbars...
% but I dont under stand why
%errorbarxy([glbMaxPosMean glbMaxPosMean+10],[glbMaxValMean glbMaxValMean+1.5],[glbMaxPosSTD   0],[glbMaxValSTD 0]  ,[],[],'sr','r')
errorbarxy([glbMaxPosMean glbMaxPosMean+5 ],[glbMaxValMean glbMaxValMean+1  ],[abs(glbMaxPosCI95(1)-glbMaxPosMean) 0],[abs(glbMaxValCI95(1)-glbMaxValMean) 0],[abs(glbMaxPosCI95(2)-glbMaxPosMean) 0],[abs(glbMaxValCI95(2)-glbMaxValMean) 0],'sr','r')
hold on
fnplt(sp,'r');
plot([0 0],[0 1],'k','LineWidth',2);
xlim([xmin xmax])
ylim([ymin ymax])
title('F / I')
xlabel('dframes')
ylabel('corr')
box on
set(gca,'LineWidth',2,'FontSize',20)
hold off

xmin=-3;
xmax=3;
figure()
hist(glbMaxPosList,linspace(xmin,xmax,101))
box on
set(gca,'LineWidth',2,'FontSize',20)
xlim([xmin xmax])

save('currentWS.mat')

corrResults.cFI     = cFI;
corrResults.cFI_std = cFI_std;
corrResults.cFF     = cFF;
corrResults.cFF_std = cFF_std;
corrResults.cII     = cII;
corrResults.cII_std = cII_std;
corrResults.par.maxLag   = maxLag;
corrResults.par.opt1     = opt1;
corrResults.par.opt2     = opt2;
corrResults.par.tBtwFrms = tBtwFrms;
corrResults.par.aveType  = aveType;
corrResults.par.normVar  = normVar;


% function[tc1,tc2]=slimTC(tc1,tc2)
% for i=1:length(tc1)
%     for t=1:length(tc1(i).observations)
%         % If either one isnan, set both to NaN:
%         if isnan(tc1(i).observations(t)) || isnan(tc2(i).observations(t))
%             tc1(i).observations(t)=NaN;
%             tc2(i).observations(t)=NaN;            
%         end
%     end
% end


function vecSum=nanSum(vec1,vec2)
if size(vec1,1)~= size(vec2,1) || size(vec1,2)~= size(vec2,2)
    error('Vectors must have the same length!')
end

for row=1:size(vec1,1)
    if ~isnan(vec1(row)) && ~isnan(vec2(row))
        vecSum(row,1)=vec1(row) + vec2(row);
    elseif isnan(vec1(row)) && isnan(vec2(row))
        vecSum(row,1)=NaN;
    elseif ~isnan(vec1(row))
        vecSum(row,1)=vec1(row);
    elseif ~isnan(vec2(row))
        vecSum(row,1)=vec2(row);
    else
        error('something went wrong');
    end
end

function vec=mat2vec(mat)
vec=mat(:);