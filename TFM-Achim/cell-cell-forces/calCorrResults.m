function [corrResults]=calCorrResults(corrSets,maxLag,opt,normVar)
% INPUT 
% opt:  'usefn': take only network forces into account. At the moment, flag
%                is only used to select the correct fm values, one could also use it to
%                select the correct resF values. At the moment, all resF
%                values are taken to get a more reliable value for the resF
%                variance.
%       'usefc': take only cluster forces into account.
%       'usefm': take network forces when possible, cluster forces
%                otherwise.

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

% now pull the right data in a new field of the data structure. Once when
% this is done, we can simmply rely on this field to calculate the
% correlations:
for idx=1:length(corrSets)
    flags=horzcat(corrSets(idx).edge.flag);
    checkVec=(sum(flags,2)==length(corrSets(idx).edge));
    for edgeId=1:length(corrSets(idx).edge)
        if strcmp(opt,'usefm')
            corrSets(idx).edge(edgeId).fcorr=corrSets(idx).edge(edgeId).fm;
        elseif strcmp(opt,'usefn')
            corrSets(idx).edge(edgeId).fcorr=NaN+zeros(length(checkVec),2);
            corrSets(idx).edge(edgeId).fcorr(checkVec,:)=corrSets(idx).edge(edgeId).fm(checkVec,:);
        elseif strcmp(opt,'usefc')
            corrSets(idx).edge(edgeId).fcorr=corrSets(idx).edge(edgeId).fc;
        else
            error('The given option is not supported');
        end
    end
end

% this shouldn't happen:
idx=1;
while idx<=length(corrSets)
    if length(corrSets(idx).edge)<2
        corrSets(idx)=[];
        display('This shouldn''t happen!')
    else
        idx=idx+1;
    end
end

%**************************************************************************
% 2) Extract the time course of the forces                                *
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
title('fi / fi{_tot}')

figure()
plotmatrix([vertcat(fi(:,1).observations) vertcat(fi(:,2).observations) vertcat(f_res(:,1).observations) vertcat(f_res(:,2).observations)])
title('fi / f_res')

figure()
plotmatrix([vertcat(fi_tot(:,1).observations) vertcat(fi_tot(:,2).observations) vertcat(f_res(:,1).observations) vertcat(f_res(:,2).observations)])
title('fi{_tot} / f{{_res}}')

figure()
title('The cross correlation for cF1Ft')
for i=1:cols
    for j=1:cols
        subplot(cols,cols,(i-1)*cols+j)
        plot(-maxLag:1:maxLag,reshape(cF1Ft(i,j,:),[],1))
        ylim([-1 1])
        xlim([-maxLag maxLag])
    end
end
title('fi / fi{_tot}')
xlabel('dframes')
ylabel('corr')


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
    end
end
title('fi / f{_res}')
xlabel('dframes')
ylabel('corr')


figure()
title('The cross correlation for cFtFr')
for i=1:cols
    for j=1:cols
        subplot(cols,cols,(i-1)*cols+j)
        plot(-maxLag:1:maxLag,reshape(cFtFr(i,j,:),[],1))
        ylim([-1 1])
        xlim([-maxLag maxLag])
    end
end
title('fi{_tot} / f{_res}')
xlabel('dframes')
ylabel('corr')


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


