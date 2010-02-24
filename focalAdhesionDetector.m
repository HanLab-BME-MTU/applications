function R = focalAdhesionDetector(I, BW, sigmaPSF, lMin)

I = double(I);

% Compute edge detection
[Re, Te, NMSe] = steerableFiltering(I, 1, sigmaPSF);

% Compute ridge detection
[Rr, Tr, NMSr, FCr] = steerableFiltering(I, 4, sigmaPSF);

R = NMSe .* NMSr .* BW;

Ifprintf('Initial number of candidates = %d\n', nnz(R));

% Remove outlier scalar product between Tr and Te (th = cos(pi/6))
cte = cos(Te);
ste = sin(Te);
ctr = cos(Tr);
str = sin(Tr);
SP = abs(cte .* ctr + ste .* str);
R(SP >= cos(pi / 6)) = 0;
fprintf('Number of candidates after removing outlier scalar products = %d\n', nnz(R));

% Amplication of response along the ridge.
idx = find(R);
[Y X] = ind2sub(size(R), idx);
CT = ctr(idx);
CT2 = CT.^2;
ST = str(idx);
ST2 = ST.^2;

% This part is for M = 2.
B = zeros(numel(idx), size(FCr, 3));
a20 = 0.16286750396763996 * sigmaPSF;
a22 = -0.4886025119029199 * sigmaPSF;
B(:, 1) = (CT2 * a20 + ST2 * a22);
B(:, 2) = CT .* ST * 2 * (a20 - a22);
B(:, 3) = (CT2 * a22 + ST * a20);

R1 = zeros(size(idx));
R2 = zeros(size(idx));
imshow(I, []); hold on;
for l = 1:lMin
    plot(X + l * CT, Y + l * ST, 'r.');
    plot(X - l * CT, Y - l * ST, 'g.');
    for iFilter = 1:size(FCr, 3)
        R1 = R1 + B(:, iFilter) .* interp2(FCr(:, :, iFilter), X + l * CT, Y + l * ST);
        R2 = R2 + B(:, iFilter) .* interp2(FCr(:, :, iFilter), X - l * CT, Y - l * ST);
    end
end
R(idx) = R(idx) + max(R1 - R2, R2 - R1);

% threshold by 2 guaussian mixture model
idx = find(R > 0);
obj = gmdistribution.fit(R(idx), 2);
[~, imin] = min(obj.mu);
c = cluster(obj, R(idx));
R(idx(c == imin)) = 0;
fprintf('Number of candidates after classification = %d\n', nnz(R));

% Display
idx = find(R);
[X Y] = ind2sub(size(R), idx);
figure, imshow(I, []); hold on;
plot(Y, X, 'g.');
quiver(Y, X, ctr(idx), str(idx), 'r'); hold on;
