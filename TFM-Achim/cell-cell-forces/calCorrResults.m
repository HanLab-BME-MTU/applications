function [corrResults]=calCorrResults(corrSets)

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

%**************************************************************************
% 3) Extract the time course of the forces                                *
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
        % the x-component:
        fi(setId,1).observations = corrSets(idx).edge(currEdge).fc(:,1); % In future, take only fc if necessary. 
        % the y-component:
        fi(setId,2).observations = corrSets(idx).edge(currEdge).fc(:,2); % In future, take only fc if necessary.

        % the x/y-components:
        fi_tot(setId,1).observations = [];
        fi_tot(setId,2).observations = [];

        % sum up all forces on the other interfaces:
        for edgeId=setxor(currEdge,1:length(corrSets(idx).edge))
            % the x/y-components:
            if isempty(fi_tot(setId,1).observations)
                fi_tot(setId,1).observations=corrSets(idx).edge(edgeId).fc(:,1);
                fi_tot(setId,2).observations=corrSets(idx).edge(edgeId).fc(:,2);
            else
                % Check for NaNs before adding the values:            
                fi_tot(setId,1).observations=nanSum(fi_tot(setId,1).observations,corrSets(idx).edge(edgeId).fc(:,1));
                fi_tot(setId,2).observations=nanSum(fi_tot(setId,2).observations,corrSets(idx).edge(edgeId).fc(:,2));
            end
        end
        % the x/y-components:
        f_res(setId,1).observations  =   corrSets(idx).resF(:,1);
        f_res(setId,2).observations  =   corrSets(idx).resF(:,2);

        setId=setId+1;
    end
end

% determine maximal possible lag:
% numPts=sum(~isnan(X(:,1)));
% maxLag=floor(numPts/4);

% calculate the correlation coefficients:
maxLag=5;

% Correlation between forces at the interface:
intx=crossCorr(fi(:,1),fi_tot(:,1),maxLag);
inty=crossCorr(fi(:,2),fi_tot(:,2),maxLag);

% Correlation between forces at the interface:
resx=crossCorr(fi(:,1),f_res(:,1),maxLag);
resy=crossCorr(fi(:,2),f_res(:,2),maxLag);

cols=2;
cF1Ft=NaN*zeros(cols,cols,2*maxLag+1);
for i=1:cols
    for j=i:cols %min(i+2,cols):cols
        out1=crossCorr(fi(:,i),fi_tot(:,j),maxLag); % in the ideal case these entries are all -1!
        cF1Ft(i,j,:)=out1(:,1);
        cF1Ft_std(i,j,:)=out1(:,2);
    end
end

cF1Fr=NaN*zeros(cols,cols,2*maxLag+1);
for i=1:cols
    for j=i:cols %min(i+2,cols):cols
        out1=crossCorr(fi(:,i),f_res(:,j),maxLag); % in the ideal case these entries are all -1!
        cF1Fr(i,j,:)=out1(:,1);
        cF1Fr_std(i,j,:)=out1(:,2);
    end
end

cFtFr=NaN*zeros(cols,cols,2*maxLag+1);
for i=1:cols
    for j=i:cols %min(i+2,cols):cols
        out1=crossCorr(fi_tot(:,i),f_res(:,j),maxLag); % in the ideal case these entries are all -1!
        cFtFr(i,j,:)=out1(:,1);
        cFtFr_std(i,j,:)=out1(:,2);
    end
end

figure(1)

cF1Ft;
cF1Fr;
cFtFr;

display(['Avg. correlation cF1Ft:',num2str(nanmean(mat2vec(cF1Ft(:,:,maxLag+1))))]);
display(['Avg. correlation cF1Fr:',num2str(nanmean(mat2vec(cF1Fr(:,:,maxLag+1))))]);
display(['Avg. correlation cFtFr:',num2str(nanmean(mat2vec(cFtFr(:,:,maxLag+1))))]);


figure(1)
plotmatrix([vertcat(fi(:,1).observations) vertcat(fi(:,2).observations) vertcat(fi_tot(:,1).observations) vertcat(fi_tot(:,2).observations)])

figure(2)
plotmatrix([vertcat(fi(:,1).observations) vertcat(fi(:,2).observations) vertcat(f_res(:,1).observations) vertcat(f_res(:,2).observations)])

figure(3)
title('The cross correlation for cF1Ft')
for i=1:cols
    for j=i:cols
        subplot(cols,cols,(i-1)*cols+j)
        plot(-maxLag:1:maxLag,reshape(cF1Ft(i,j,:),[],1))
        ylim([-1 1])
        xlim([-maxLag maxLag])
    end
end

figure(4)
title('The cross correlation for cF1Fr')
for i=1:cols
    for j=i:cols
        subplot(cols,cols,(i-1)*cols+j)
        plot(-maxLag:1:maxLag,reshape(cF1Fr(i,j,:),[],1))
        ylim([-1 1])
        xlim([-maxLag maxLag])
    end
end

figure(5)
title('The cross correlation for cFtFr')
for i=1:cols
    for j=i:cols
        subplot(cols,cols,(i-1)*cols+j)
        plot(-maxLag:1:maxLag,reshape(cFtFr(i,j,:),[],1))
        ylim([-1 1])
        xlim([-maxLag maxLag])
    end
end

corrResults.cF1Ft = cF1Ft;
corrResults.cF1Fr = cF1Fr;
corrResults.cFtFr = cFtFr;
corrResults.maxLag= maxLag;




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


