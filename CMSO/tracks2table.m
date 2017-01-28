function T=tracks2table(tracks)
% P. Roudot 2016
varName={'trackId','t','X','Y','Z','A','dX','dY','dZ','dA'};
tableCell=cell(1,length(tracks));
parfor tIdx=1:length(tracks)
    if(isprop(tracks(tIdx),'z'))
        Z=tracks(tIdx).z;
        dZ=tracks(tIdx).dz;
    else
        Z=zeros(size((tracks(tIdx).z)));
        dZ=Z;
    end
    detArray=[tIdx*ones(size(Z')), tracks(tIdx).t', tracks(tIdx).x', tracks(tIdx).y',Z',tracks(tIdx).A',tracks(tIdx).dx', tracks(tIdx).dy',dZ',tracks(tIdx).dA'];
    tableCell{tIdx}=array2table(detArray,'variableNames',varName);
end
T=table();
for tIdx=1:length(tracks)
    T=[T;tableCell{tIdx}];
end