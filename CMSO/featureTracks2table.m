function T=featureTracks2table(tracks)
% P. Roudot 2016
varName={'trackId','t','f','featureId'};
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
    preAllocArray(trackStartIdx(tIdx):(trackStartIdx(tIdx+1)-1),:)=[tIdx*ones(size(Z')), tracks(tIdx).t', tracks(tIdx).f', tracks(tIdx).tracksFeatIndxCG'];
end
T=array2table(preAllocArray,'variableNames',varName);
