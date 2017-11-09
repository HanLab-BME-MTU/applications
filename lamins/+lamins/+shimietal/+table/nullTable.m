function [T, data, pvalue, Ncells, Nfaces] = nullTable(stats)

labels = {'LAC', 'LB1', 'LB2'};

WT.idx = 20;
% channels for comparison
% Use LAC channel for LB nulls
WT.LB1 = 1;
WT.LB2 = 1;
% Use LB12 channel for LAC nulls
WT.LAC = 2;

idx = [2 7 12];

% goodLB2null = [stats(12).area(1).movie([1 2 4 5 6 8 9 10]).data];

fields = {'edgesPerFace','avgEdgeLength','perimeter','circularity','area'};

data = zeros(length(idx),length(fields));
pvalue = zeros(length(idx),length(fields));
rowNames = strcat(labels,' null / WT');
Ncells = zeros(length(idx),2);
Nfaces = zeros(length(idx),length(fields),1);

for i=1:length(idx)
    for j=1:length(fields)
        field = fields{j};
        Xi.idx = WT.idx;
        Xi.ch = WT.(labels{i});
        X = stats(Xi.idx).(field)(Xi.ch).all.data;
        Yi.idx = idx(i);
        if(i == 3)
            Y = [stats(Yi.idx).(field)(1).movie([1 2 4 5 6 8 9 10]).data];
        else
            Y = stats(Yi.idx).(field)(1).all.data;
        end
        data(i,j) = lamins.functions.qqscale(X,Y);
        pvalue(i,j) = ranksum(X,Y);
        Nfaces(i,j,1) = length(Y);
        Nfaces(i,j,2) = length(X);
    end
    if(isfield(stats(Xi.idx).(field)(Xi.ch),'movie'))
        Ncells(i,1) = length(stats(Yi.idx).(field)(1).movie);
        Ncells(i,2) = length(stats(Xi.idx).(field)(Xi.ch).movie);
    else
        Ncells(i,:) = NaN;
    end
end

Ncells(3,1) = 8;

T = array2table(data,'RowNames',rowNames,'Variablenames',fields);

end