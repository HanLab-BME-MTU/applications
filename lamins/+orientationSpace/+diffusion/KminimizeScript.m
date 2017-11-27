maxima = R.getRidgeOrientationLocalMaxima;
import orientationSpace.diffusion.*;
[ L.maxima, L.coords, L.response] = linearizeMaxima( R, maxima*2);
% [ d.maxima, d.coords, d.response, d.xgd ] = descendHalleyK( L.maxima, L.coords, L.response );
pool_size = 12;
block_size = ceil(size(L.maxima,2)/pool_size);
partition = [repmat(block_size,1,pool_size-1) mod(size(L.maxima,2),block_size)];
LC.maxima = mat2cell(L.maxima,1,partition);
LC.response = mat2cell(L.response,size(L.response,1),partition);
LC.coords.r = mat2cell(L.coords.r,1,partition);
LC.coords.c = mat2cell(L.coords.c,1,partition);
LC.coords.m = mat2cell(L.coords.m,1,partition);
LC.coords.K = mat2cell(L.coords.K,1,partition);
LC.coords.ind = mat2cell(L.coords.ind,1,partition);
LC.coords = cellfun(@(r,c,m,K,ind) struct('r',r,'c',c,'m',m,'K',K,'ind',ind),LC.coords.r,LC.coords.c,LC.coords.m,LC.coords.K,LC.coords.ind,'UniformOutput',false);
[ d.maxima, d.coords, d.response ] = parcellfun_progress(@orientationSpace.diffusion.descendHalleyK,LC.maxima,LC.coords,LC.response,'UniformOutput',false);
% d.maxima = [d.maxima{:}];
% d.response = [d.response{:}];
d.coords = [d.coords{:}];
% d.coords = struct('r',[d.coords.r],'c',[d.coords.c],'m',[d.coords.m],'K',[d.coords.K],'ind',[d.coords.ind]);
% clear LC
orgK = repmat({8},1,pool_size);

% [x_out,K_out,notDone] = refineBifurcation(d.maxima,d.coords.K,L.response,8);
[x_out,K_out,notDone] = parcellfun_progress(@orientationSpace.diffusion.refineBifurcation,d.maxima,{d.coords.K},LC.response,orgK,'UniformOutput',false);
x_out = [x_out{:}];
K_out = [K_out{:}];
notDone = [notDone{:}];
clear LC;

maximaF = NaN(1024,1024,7);
L.coords.ind = sub2ind([1024 1024 7],L.coords.r,L.coords.c,L.coords.m);
maximaF(L.coords.ind) = x_out;
[fr,fc] = find(sum(maximaF(:,:,1:end-1) == maximaF(:,:,2:end),3));
maximaF_repeats = cat(3,false(1024),maximaF(:,:,1:end-1) == maximaF(:,:,2:end));
nn = find(maximaF_repeats(L.coords.ind));
x_fixed = x_out;
K_fixed = K_out;
for n = nn
[x_fixed(n),K_fixed(n)] = orientationSpace.diffusion.newtonBPproto(R, n, L.coords, L.maxima(n), 8);
end
maximaFixed = NaN(1024,1024,7);
maximaFixed(L.coords.ind) = x_fixed;
maximaF_repeats = cat(3,maximaF(:,:,1:end-1) == maximaF(:,:,2:end),false(1024));
nn = find(maximaF_repeats(L.coords.ind));
for n = nn
[x_fixed(n),K_fixed(n)] = orientationSpace.diffusion.newtonBPproto(R, n, L.coords, L.maxima(n), 8);
end
save('maximaFixed.mat','maximaFixed')

maximaFixedN = NaN(1024,1024,7);
maximaFixedN(coords.ind) = x_fixed;
maximaFixedN = maximaFixedN/2;

K_fixedN = NaN(1024,1024,7);
K_fixedN(coords.ind) = K_fixed;

R3 = R.getResponseAtOrderFT(3,2);
maxima3 = R3.getRidgeOrientationLocalMaxima;

maximaFixedN(K_fixedN <= 3) = NaN;
maximaCombined = cat(3,maximaFixedN,maxima3);
maximaCombined = sort(maximaCombined,3);
max(max(sum(~isnan(maximaCombined),3)))
maximaCombined = maximaCombined(:,:,1:7);

maximaCombinedValues = R3.interpft1(maximaCombined);
maximaCombinedValues(isnan(maximaCombinedValues)) = -Inf;
[maximaCombinedValues,maximaCombined] = sortMatrices(maximaCombinedValues,maximaCombined,3,'descend');
maximaCombinedValues(isinf(maximaCombinedValues)) = -NaN;

% Minimal Bridging process

% nlms = R3.nonLocalMaximaSuppressionPrecise(maximaCombined);
[nlms,nlms_offset] = nonLocalMaximaSuppressionPrecise(R3.a,maximaCombined);


nlms_mip3 = nanmax(nlms,[],3);

% mask = imopen(imfill(nlms_mip3 > 10,'holes'),strel('disk',10));
mask = lamins.functions.maskFromSteerable(R3);

selected = nlms_mip3(mask & nlms_mip3 ~= 0);
[outliers,inliers] = detectOutliers(nlms_mip3(mask & nlms_mip3 ~= 0));
T = thresholdRosin(selected(inliers));

k = 1;
nms = nlms(:,:,1);
% figure; imshow(nms,[]);
nlms_binary = (nlms_mip3 > T) & mask;
nms_binary = (nms > T) & mask;
nms_skel = bwmorph(nms_binary,'skel',Inf);
nms_skel_thresh = nms_skel.* nlms_binary;
nms_skel_thresh_wo_bp = lamins.functions.bwRemoveBranchPoints(nms_skel_thresh);
nlms_fragments = nlms_binary & ~nms_skel_thresh_wo_bp;
nlms_fragments_cc = bwconncomp(nlms_fragments);
nms_skel_wo_bp_cc = bwconncomp(nms_skel_thresh_wo_bp);
out.bridges{k} = lamins.functions.minimalBridge(nlms_fragments_cc,nms_skel_wo_bp_cc,I);
nms_skel_bridged = nms_skel | out.bridges{k};
out.nms_skel_bridged_skel{k} = bwmorph(nms_skel_bridged,'skel',Inf);