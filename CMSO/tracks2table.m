function T=tracks2table(tracks)
% P. Roudot 2016
varName={'trackId','t','X','Y','Z','A','dX','dY','dZ','dA'};
trackStartIdx=[1,cumsum([tracks.lifetime])+1];
preAllocArray=zeros(trackStartIdx(end)-1,length(varName));
for tIdx=1:length(tracks)
    if(isprop(tracks(tIdx),'z'))
        Z=tracks(tIdx).z;
        dZ=tracks(tIdx).dz;
    else
        Z=zeros(size((tracks(tIdx).z)));
        dZ=Z;
    end
    preAllocArray(trackStartIdx(tIdx):(trackStartIdx(tIdx+1)-1),:)=[tIdx*ones(size(Z')), tracks(tIdx).t', tracks(tIdx).x', tracks(tIdx).y',Z',tracks(tIdx).A',tracks(tIdx).dx', tracks(tIdx).dy',dZ',tracks(tIdx).dA'];
end
T=array2table(preAllocArray,'variableNames',varName);
