function simulateDistance

L = imread('/orchestra/groups/lccb-interfil/Tropo/Data_for_paper/TM2_Actin/cell5/windowAnalysis/labels/labels_021.tif');
BW = imread('/orchestra/groups/lccb-interfil/Tropo/Data_for_paper/TM2_Actin/cell5/Actin/edge/cell_mask/mask_default21.tif');
D = double(bwdist(max(BW(:)) - BW)) * 67;

labels = 1:max(L(:));


idx = find(D < 5000 & L ~= 0);
n = numel(idx);
p = randperm(n);
idx = idx(p(1:ceil(.3 * n)));
[y x] = ind2sub(size(BW), idx);

%imshow(L, []); hold on;
%plot(x, y, 'r.');

expMean = mean(D(idx));

sectorIdx = arrayfun(@(l) idx(L(idx) == l), ...
    1:max(L(:)), 'UniformOutput', false);

sectorMean = cellfun(@(i) mean(D(i)), ...
    sectorIdx, 'UniformOutput', false);
sectorMean = horzcat(sectorMean{:});

plot(1:numel(labels), repmat(expMean, 1, numel(labels)), 'r'); hold on;
plot(1:numel(labels), sectorMean, 'b');

w = arrayfun(@(l) mean(D(L == l)), labels);
w = expMean - w;

plot(1:numel(labels), sectorMean + w, 'g');

legend({'Expected distance to edge', 'Sectorized average distance to edge', ...
    'Shifted sectorized distance to edge'});

% w = (2500/67) ./ w;
% 
% W = zeros(size(D));
% 
% for i = 1:numel(w)
%     l = labels(i);
%     W(L == l) = w(l);
% end
% 
% imshow(W, []);