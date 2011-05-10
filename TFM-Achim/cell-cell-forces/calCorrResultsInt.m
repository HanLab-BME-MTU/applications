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
    opt1='useItot';
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



figure()
title('The cross correlation for cFI')
errorbar(-maxLag:1:maxLag,cFI,cFI_std,'r')
hold on
plot(-maxLag:1:maxLag,cFI,'k')
% ylim([-1 1])
xlim([-maxLag maxLag])
title('F / I')
xlabel('dframes')
ylabel('corr')
hold off

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