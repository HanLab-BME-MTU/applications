function [corrResults]=calCorrResults(corrSets,maxLag,opt,normVar,tBtwFrms,aveType)
% INPUT 
% opt:  'usefn': take only network forces into account. At the moment, flag
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
% aveType:       In case tBtwFrms is not empty the user has a few options
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

if nargin<3 || isempty(opt)
    opt='usefm';
end

if nargin<4 || isempty(normVar)
    normVar=0;
end

if nargin<6 || isempty(aveType)
    aveType='none';
end

%**************************************************************************
% 1) Pull the right data in new fields of the data structure. 
%**************************************************************************
%Once when this is done, we can simply rely on this field to calculate the
% correlations:
corrSets=pullCorrVals(corrSets,opt);

    

% this shouldn't happen:
idx=1;
while idx<=length(corrSets)
    if length(corrSets(idx).edge)<2
        corrSets(idx)=[];
        display('This shouldn''t happen! Cell found with less than two interfaces!')
    else
        idx=idx+1;
    end
end


%**************************************************************************
% 2) adapt the frame rates:    
%**************************************************************************
if nargin>=5 && ~isempty(tBtwFrms)
    corrSets=adaptFrameRates(corrSets,tBtwFrms,aveType);
end


%**************************************************************************
% 3) Extract the time course of the forces                                *
%**************************************************************************
% run through all clusters and collect the data.
currEdge=1;
setId=1;
for idx=1:length(corrSets)
    % Does this make sense for cells with only two edges? Shouldn't these
    % values be then more or less the same?
%     if length(corrSets(idx).edge)>2
%         edgePerm=1:length(corrSets(idx).edge);
%     else
%          edgePerm=1;
%     end
    
    edgePerm=1:length(corrSets(idx).edge);
    edgePerm(end)=[]; % skip the last value since it is dependent on the previous results?! At least this is obviously true for numEdge=2;
    
    for currEdge=edgePerm
        % 1/2=the x/y-component:
        fi(setId,1).observations = corrSets(idx).edge(currEdge).fcorr(:,1); 
        fi(setId,2).observations = corrSets(idx).edge(currEdge).fcorr(:,2);

        % the x/y-components:
        fi_tot(setId,1).observations = [];
        fi_tot(setId,2).observations = [];

        % sum up all forces on the other interfaces:
        for edgeId=setxor(currEdge,1:length(corrSets(idx).edge))
            % the x/y-components:
            if isempty(fi_tot(setId,1).observations)
                fi_tot(setId,1).observations=corrSets(idx).edge(edgeId).fcorr(:,1);
                fi_tot(setId,2).observations=corrSets(idx).edge(edgeId).fcorr(:,2);
            else
                % Check for NaNs before adding the values:            
                fi_tot(setId,1).observations=nanSum(fi_tot(setId,1).observations,corrSets(idx).edge(edgeId).fcorr(:,1));
                fi_tot(setId,2).observations=nanSum(fi_tot(setId,2).observations,corrSets(idx).edge(edgeId).fcorr(:,2));
            end
        end
        % the x/y-components:
        f_res(setId,1).observations  =   corrSets(idx).resF(:,1);
        f_res(setId,2).observations  =   corrSets(idx).resF(:,2);

        setId=setId+1;
    end
end

cols=2;
% Calculate the autocorrelation:
cF1F1=NaN*zeros(cols,cols,2*maxLag+1);
for i=1:cols
    for j=1:cols %min(i+2,cols):cols
        out1=crossCorr(fi(:,i),fi(:,j),maxLag,normVar); % in the ideal case these entries are all -1!
        cF1F1(i,j,:)=out1(:,1);
        cF1F1_std(i,j,:)=out1(:,2);
    end
end

cFtFt=NaN*zeros(cols,cols,2*maxLag+1);
for i=1:cols
    for j=1:cols %min(i+2,cols):cols
        out1=crossCorr(fi_tot(:,i),fi_tot(:,j),maxLag,normVar); % in the ideal case these entries are all -1!
        cFtFt(i,j,:)=out1(:,1);
        cFtFt_std(i,j,:)=out1(:,2);
    end
end

cFrFr=NaN*zeros(cols,cols,2*maxLag+1);
for i=1:cols
    for j=1:cols %min(i+2,cols):cols
        out1=crossCorr(f_res(:,i),f_res(:,j),maxLag,normVar); % in the ideal case these entries are all -1!
        cFrFr(i,j,:)=out1(:,1);
        cFrFr_std(i,j,:)=out1(:,2);
    end
end

% Calculate the cross correlation:
cF1Ft=NaN*zeros(cols,cols,2*maxLag+1);
for i=1:cols
    for j=1:cols %min(i+2,cols):cols
        out1=crossCorr(fi(:,i),fi_tot(:,j),maxLag,normVar); % in the ideal case these entries are all -1!
        cF1Ft(i,j,:)=out1(:,1);
        cF1Ft_std(i,j,:)=out1(:,2);
    end
end

if doPlotTraces
    for k=1:length(fi(:,1))
        figure()
        plot(1:length(fi(k,1).observations),fi(k,1).observations,'b');
        hold on;
        plot(1:length(fi_tot(k,1).observations),-fi_tot(k,1).observations,'r');
        hold off;
        set(gca,'fontsize',16)
    end
end

cF1Fr=NaN*zeros(cols,cols,2*maxLag+1);
for i=1:cols
    for j=1:cols %min(i+2,cols):cols
        out1=crossCorr(fi(:,i),f_res(:,j),maxLag,normVar); % in the ideal case these entries are all -1!
        cF1Fr(i,j,:)=out1(:,1);
        cF1Fr_std(i,j,:)=out1(:,2);
    end
end

cFtFr=NaN*zeros(cols,cols,2*maxLag+1);
for i=1:cols
    for j=1:cols %min(i+2,cols):cols
        out1=crossCorr(fi_tot(:,i),f_res(:,j),maxLag,normVar); % in the ideal case these entries are all -1!
        cFtFr(i,j,:)=out1(:,1);
        cFtFr_std(i,j,:)=out1(:,2);
    end
end

cF1Ft;
cF1Fr;
cFtFr;

display(['Avg. XXXX correlation cF1Ft, mean all: ',num2str(nanmean(mat2vec(cF1Ft(:,:,maxLag+1)))),'; mean trace: ',num2str(trace(cF1Ft(:,:,maxLag+1))/2),' (+-',num2str(trace(cF1Ft_std(:,:,maxLag+1))/2),')']);
display(['Avg. XXXX correlation cF1Fr, mean all: ',num2str(nanmean(mat2vec(cF1Fr(:,:,maxLag+1)))),'; mean trace: ',num2str(trace(cF1Fr(:,:,maxLag+1))/2),' (+-',num2str(trace(cF1Fr_std(:,:,maxLag+1))/2),')']);
display(['Avg. XXXX correlation cFtFr, mean all: ',num2str(nanmean(mat2vec(cFtFr(:,:,maxLag+1)))),'; mean trace: ',num2str(trace(cFtFr(:,:,maxLag+1))/2),' (+-',num2str(trace(cFtFr_std(:,:,maxLag+1))/2),')']);

display(['Avg. Auto correlation cF1F1, mean all: ',num2str(nanmean(mat2vec(cF1F1(:,:,maxLag+1)))),'; mean trace: ',num2str(trace(cF1F1(:,:,maxLag+1))/2),' (+-',num2str(trace(cF1F1_std(:,:,maxLag+1))/2),')']);
display(['Avg. Auto correlation cFtFt, mean all: ',num2str(nanmean(mat2vec(cFtFt(:,:,maxLag+1)))),'; mean trace: ',num2str(trace(cFtFt(:,:,maxLag+1))/2),' (+-',num2str(trace(cFtFt_std(:,:,maxLag+1))/2),')']);
display(['Avg. Auto correlation cFrFr, mean all: ',num2str(nanmean(mat2vec(cFrFr(:,:,maxLag+1)))),'; mean trace: ',num2str(trace(cFrFr(:,:,maxLag+1))/2),' (+-',num2str(trace(cFrFr_std(:,:,maxLag+1))/2),')']);



figure()
plotmatrix([vertcat(fi(:,1).observations) vertcat(fi(:,2).observations) vertcat(fi_tot(:,1).observations) vertcat(fi_tot(:,2).observations)])
title('f_i / f_{i_{tot}}')

figure()
plotmatrix([vertcat(fi(:,1).observations) vertcat(fi(:,2).observations) vertcat(f_res(:,1).observations) vertcat(f_res(:,2).observations)])
title('f_i / f_{res}')

figure()
plotmatrix([vertcat(fi_tot(:,1).observations) vertcat(fi_tot(:,2).observations) vertcat(f_res(:,1).observations) vertcat(f_res(:,2).observations)])
title('f_{i_{tot}} / f_{res}')

figure()
title('The autco-rrelation for cF1F1')
for i=1:cols
    for j=1:cols
        subplot(cols,cols,(i-1)*cols+j)
        plot(-maxLag:1:maxLag,reshape(cF1F1(i,j,:),[],1))
        ylim([-1 1])
        xlim([-maxLag maxLag])
        title('Autocorrelation f_i f_i')
        xlabel('dframes')
        ylabel('corr')
    end
end

figure()
title('The auto-correlation for cFrFr')
for i=1:cols
    for j=1:cols
        subplot(cols,cols,(i-1)*cols+j)
        plot(-maxLag:1:maxLag,reshape(cFrFr(i,j,:),[],1))
        ylim([-1 1])
        xlim([-maxLag maxLag])
        title('Autocorrelation f_r f_r')
        xlabel('dframes')
        ylabel('corr')
    end
end

figure()
title('The cross correlation for cF1Ft')
for i=1:cols
    for j=1:cols
        subplot(cols,cols,(i-1)*cols+j)
        plot(-maxLag:1:maxLag,reshape(cF1Ft(i,j,:),[],1))
        ylim([-1 1])
        xlim([-maxLag maxLag])
        title('f_i / f_{i_{tot}}')
        xlabel('dframes')
        ylabel('corr')
    end
end


if doPlotTraces
    figure()
    plot(-maxLag:1:maxLag,reshape(cF1Ft(1,1,:),[],1))
    xlim([-maxLag maxLag])
    set(gca,'fontsize',16)
end



figure()
title('The cross correlation for cF1Fr')
for i=1:cols
    for j=1:cols
        subplot(cols,cols,(i-1)*cols+j)
        plot(-maxLag:1:maxLag,reshape(cF1Fr(i,j,:),[],1))
        ylim([-1 1])
        xlim([-maxLag maxLag])
        title('f_i / f_{res}')
        xlabel('dframes')
        ylabel('corr')
    end
end



figure()
title('The cross correlation for cFtFr')
for i=1:cols
    for j=1:cols
        subplot(cols,cols,(i-1)*cols+j)
        plot(-maxLag:1:maxLag,reshape(cFtFr(i,j,:),[],1))
        ylim([-1 1])
        xlim([-maxLag maxLag])
        title('f_{i_{tot}} / f_{res}')
        xlabel('dframes')
        ylabel('corr')
    end
end



corrResults.cF1Ft     = cF1Ft;
corrResults.cF1Ft_std = cF1Ft_std;
corrResults.cF1Fr     = cF1Fr;
corrResults.cF1Fr_std = cF1Fr_std;
corrResults.cFtFr     = cFtFr;
corrResults.cFtFr_std = cFtFr_std;
corrResults.cF1F1     = cF1F1;
corrResults.cF1F1_std = cF1F1_std;
corrResults.cFtFt     = cFtFt;
corrResults.cFtFt_std = cFtFt_std;
corrResults.cFrFr     = cFrFr;
corrResults.cFrFr_std = cFrFr_std;
corrResults.par.maxLag = maxLag;
corrResults.par.opt    = opt;
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


