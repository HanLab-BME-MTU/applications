function [result] = segmentLbpImage(Ilbp,patchSize)
[sizeY,sizeX] = size(Ilbp);

for l = 0 : 9
    fun = @(block_struct) ...
        sum(sum(block_struct.data == l))/(patchSize^2);
    tmp = blockproc(Ilbp,[patchSize patchSize],fun);
    if ~exist('Ihist','var')
        Ihist = nan(size(tmp,1),size(tmp,2),10);
    end
    Ihist(:,:,l+1) = tmp;
%     Il = ['I' int2str(l)];
%     eval([Il '= tmp;']);
end

blocksHist = reshape(Ihist,size(Ihist,1)*size(Ihist,2),10);

% [pc,score,latent,tsquare] = princomp(blocksHist);

[kinds, kcenters] = kmeans(blocksHist,2);
kdistances = pdist2(blocksHist,kcenters);
[kWrongDist] = max(kdistances,[],2);
[kRightDist] = min(kdistances,[],2);
kdist = (kWrongDist - kRightDist)./kWrongDist;


result.IHist = Ihist;
result.pacthHist = blocksHist;
result.kmeans.labels = kinds;
result.kmeans.centers = kcenters;
result.kmeans.distances = kdistances;
result.kmeans.normDistances = kdist;
% result.pca.score = score;
% result.pca.pc = pc;
% result.pca.latent = latent;
end