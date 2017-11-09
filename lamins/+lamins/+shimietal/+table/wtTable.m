function [T, data, pvalue, Ncells, Nfaces] = wtTable(stats)

idx = [21 22 25 27 28 23]';
                
hist.labels = {  ...
   'LA'  'LB1';  ...
   'LA'  'LB2';  ...
   'LB2' 'LB1';  ...
   'LC'  'LB1';  ...
   'LC'  'LB2';  ...
   'LC'  'LA' };

hist.labels = reshape(hist.labels,[],2);

fields = {'edgesPerFace','avgEdgeLength','perimeter','circularity','area'};

rowLabels = strcat(hist.labels(:,1),'/',hist.labels(:,2));

data = zeros(length(idx),length(fields));
pvalue = zeros(size(data));
Ncells = zeros(length(idx),2);
Nfaces = zeros(length(idx),length(fields),1);

for i=1:length(idx)
    for j=1:length(fields)
        field = fields{j};
        Xi.idx = idx(i);
        X = stats(Xi.idx).(field)(1).all.data;
        Yi.idx = idx(i);
        Y = stats(Yi.idx).(field)(2).all.data;
        data(i,j) = lamins.functions.qqscale(X,Y);
        pvalue(i,j) = ranksum(X,Y);
        Nfaces(i,j,1) = length(X);
        Nfaces(i,j,2) = length(Y);
    end
    if(isfield(stats(Xi.idx).(field)(1),'movie'))
        Ncells(i,1) = length(stats(Xi.idx).(field)(1).movie);
        Ncells(i,2) = length(stats(Yi.idx).(field)(2).movie);
    else
        Ncells(i,:) = NaN;
    end
end

T = array2table(data,'RowNames',rowLabels,'Variablenames',fields);

end