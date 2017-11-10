function simulateDistance(BW, L)

D = double(bwdist(max(BW(:)) - BW)) * 67;

labels = 1:max(L(:));


idx = find(D < 5000 & L ~= 0);
n = numel(idx);
p = randperm(n);
idx = idx(p(1:ceil(.05 * n)));
[y x] = ind2sub(size(BW), idx);

imshow(BW, []); hold on;
plot(x, y, 'r.');

expMean = mean(D(idx));

sectorIdx = arrayfun(@(l) idx(L(idx) == l), ...
    1:max(L(:)), 'UniformOutput', false);

sectorMean = cellfun(@(i) mean(D(i)), ...
    sectorIdx, 'UniformOutput', false);
sectorMean = horzcat(sectorMean{:});

figure;
plot(1:numel(labels), repmat(expMean, 1, numel(labels)), 'r'); hold on;
plot(1:numel(labels), sectorMean, 'b');

w = arrayfun(@(l) mean(D(L == l)), labels);
w = expMean - w;

plot(1:numel(labels), sectorMean + w, 'g');
xlabel('# sector');
ylabel('distance to edge (nm)');
legend({'Expected', 'Sectorized', ...
    'Corrected sectorized'});

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