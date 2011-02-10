function [corrResults]=calCorrResults(corrSets,maxLag,opt)
% INPUT 
% opt:  'usefn': take only network forces into account.
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

% calculate the correlation coefficients:
if nargin<2 || isempty(maxLag)
    maxLag=5;
end

if nargin<3 || isempty(opt)
    opt='usefm';
end

% now pull the right data in a new field of the data structure. Once when
% this is done, we can simmply rely on this field to calculate the
% correlations:
for idx=1:length(corrSets)
    for edgeId=1:length(corrSets(idx).edge)
        if strcmp(opt,'usefm')
            corrSets(idx).edge(edgeId).fcorr=corrSets(idx).edge(edgeId).fm;
        elseif strcmp(opt,'usefn')
            corrSets(idx).edge(edgeId).fcorr=corrSets(idx).edge(edgeId).fn;
        elseif strcmp(opt,'usefc')
            corrSets(idx).edge(edgeId).fcorr=corrSets(idx).edge(edgeId).fc;
        else
            error('The given option is not supported');
        end
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
    if length(corrSets(idx).edge)>2
        edgePerm=1:length(corrSets(idx).edge);
    else
        edgePerm=1;
    end
    
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
        out1=crossCorr(fi(:,i),fi(:,j),maxLag); % in the ideal case these entries are all -1!
        cF1F1(i,j,:)=out1(:,1);
        cF1F1_std(i,j,:)=out1(:,2);
    end
end

cFtFt=NaN*zeros(cols,cols,2*maxLag+1);
for i=1:cols
    for j=1:cols %min(i+2,cols):cols
        out1=crossCorr(fi_tot(:,i),fi_tot(:,j),maxLag); % in the ideal case these entries are all -1!
        cFtFt(i,j,:)=out1(:,1);
        cFtFt_std(i,j,:)=out1(:,2);
    end
end

cFrFr=NaN*zeros(cols,cols,2*maxLag+1);
for i=1:cols
    for j=1:cols %min(i+2,cols):cols
        out1=crossCorr(f_res(:,i),f_res(:,j),maxLag); % in the ideal case these entries are all -1!
        cFrFr(i,j,:)=out1(:,1);
        cFrFr_std(i,j,:)=out1(:,2);
    end
end


% Calculate the cross correlation:
cF1Ft=NaN*zeros(cols,cols,2*maxLag+1);
for i=1:cols
    for j=1:cols %min(i+2,cols):cols
        out1=crossCorr(fi(:,i),fi_tot(:,j),maxLag); % in the ideal case these entries are all -1!
        cF1Ft(i,j,:)=out1(:,1);
        cF1Ft_std(i,j,:)=out1(:,2);
    end
end

cF1Fr=NaN*zeros(cols,cols,2*maxLag+1);
for i=1:cols
    for j=1:cols %min(i+2,cols):cols
        out1=crossCorr(fi(:,i),f_res(:,j),maxLag); % in the ideal case these entries are all -1!
        cF1Fr(i,j,:)=out1(:,1);
        cF1Fr_std(i,j,:)=out1(:,2);
    end
end

cFtFr=NaN*zeros(cols,cols,2*maxLag+1);
for i=1:cols
    for j=1:cols %min(i+2,cols):cols
        out1=crossCorr(fi_tot(:,i),f_res(:,j),maxLag); % in the ideal case these entries are all -1!
        cFtFr(i,j,:)=out1(:,1);
        cFtFr_std(i,j,:)=out1(:,2);
    end
end

cF1Ft;
cF1Fr;
cFtFr;

display(['Avg. XXXX correlation cF1Ft:',num2str(nanmean(mat2vec(cF1Ft(:,:,maxLag+1))))]);
display(['Avg. XXXX correlation cF1Fr:',num2str(nanmean(mat2vec(cF1Fr(:,:,maxLag+1))))]);
display(['Avg. XXXX correlation cFtFr:',num2str(nanmean(mat2vec(cFtFr(:,:,maxLag+1))))]);

display(['Avg. Auto correlation cF1F1:',num2str(nanmean(mat2vec(cF1F1(:,:,maxLag+1))))]);
display(['Avg. Auto correlation cFtFt:',num2str(nanmean(mat2vec(cFtFt(:,:,maxLag+1))))]);
display(['Avg. Auto correlation cFrFr:',num2str(nanmean(mat2vec(cFrFr(:,:,maxLag+1))))]);



figure()
plotmatrix([vertcat(fi(:,1).observations) vertcat(fi(:,2).observations) vertcat(fi_tot(:,1).observations) vertcat(fi_tot(:,2).observations)])

figure()
plotmatrix([vertcat(fi(:,1).observations) vertcat(fi(:,2).observations) vertcat(f_res(:,1).observations) vertcat(f_res(:,2).observations)])

figure()
plotmatrix([vertcat(fi_tot(:,1).observations) vertcat(fi_tot(:,2).observations) vertcat(f_res(:,1).observations) vertcat(f_res(:,2).observations)])

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

corrResults.cF1Ft = cF1Ft;
corrResults.cF1Fr = cF1Fr;
corrResults.cFtFr = cFtFr;
corrResults.cF1F1 = cF1F1;
corrResults.cFtFt = cFtFt;
corrResults.cFrFr = cFrFr;
corrResults.par.maxLag = maxLag;
corrResults.par.opt    = opt;




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


