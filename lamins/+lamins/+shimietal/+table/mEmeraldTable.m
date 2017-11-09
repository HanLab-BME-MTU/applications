function [T, data, pvalue, Ncells, Nfaces] = mEmeraldTable(stats)

immuno.A.idx = 21;
immuno.B1.idx = 21;
immuno.B2.idx = 28;
immuno.C.idx = 28;

immuno.A.ch = 1;
immuno.B1.ch = 2;
immuno.B2.ch = 2;
immuno.C.ch = 1;

mEmerald.A.idx = 29;
mEmerald.B1.idx = 30;
mEmerald.B2.idx = 31;
mEmerald.C.idx = 32;

mEmerald.A.ch = 1;
mEmerald.B1.ch = 1;
mEmerald.B2.ch = 1;
mEmerald.C.ch = 1;

fields = {'edgesPerFace','avgEdgeLength','perimeter','circularity','area'};

subtypes = {'A','B1','B2','C'};

data = zeros(length(subtypes),length(fields));
pvalue = zeros(size(data));
Ncells = zeros(length(subtypes),2);
Nfaces = zeros(length(subtypes),length(fields),1);

for i=1:length(subtypes)
    for j=1:length(fields)
        field = fields{j};
        Xi = mEmerald.(subtypes{i});
        X = stats(Xi.idx).(field)(Xi.ch).all.data;
        Yi = immuno.(subtypes{i});
        Y = stats(Yi.idx).(field)(Yi.ch).all.data;
        data(i,j) = lamins.functions.qqscale(X,Y);
        pvalue(i,j) = ranksum(X,Y);
        Nfaces(i,j,1) = length(X);
        Nfaces(i,j,2) = length(Y);
    end
    if(isfield(stats(Xi.idx).(field)(Xi.ch),'movie'))
        Ncells(i,1) = length(stats(Xi.idx).(field)(Xi.ch).movie);
        Ncells(i,2) = length(stats(Yi.idx).(field)(Yi.ch).movie);
    else
        Ncells(i,:) = NaN;
    end
end

T = array2table(data,'RowNames',subtypes,'Variablenames',fields);

end