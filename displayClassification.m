function displayClassification(loDist, hiDist, pvalue, rootDirectory)

loDist = loDist * 1000 / 67;
hiDist = hiDist * 1000 / 67;

load([rootDirectory filesep 'globalCorrelations.mat']);

indDist1 = find(R{1}(:, 3) >= 75 & R{1}(:, 3) <= 150 & R{1}(:, 6) > 1); %#ok<USENS,COLND>
indDist2 = find(R{2}(:, 3) >= 75 & R{2}(:, 3) <= 150 & R{2}(:, 6) > 1);

X = vertcat(R{1}(indDist1, 4:5), R{2}(indDist2, 4:5)); %#ok<FNDSB>
obj = gmdistribution(mean(X), cov(X), 1);
D = mahal(obj,X);

% figure
% scatter(X(:,1),X(:,2),10,D(:,1),'.')
% colorbar;

indDist1 = find(R{1}(:, 3) >= loDist & R{1}(:, 3) <= hiDist & R{1}(:, 6) > 1);
indDist2 = find(R{2}(:, 3) >= loDist & R{2}(:, 3) <= hiDist & R{2}(:, 6) > 1);

Y = vertcat(R{1}(indDist1, :), R{2}(indDist2, :)); %#ok<FNDSB>
D = mahal(obj,Y(:, 4:5));

% figure,
% scatter(Y(:,4),Y(:,5),10,D(:,1),'.')
% colorbar;

dCrit = icdf('chi2', pvalue, 1);
ind1 = find(D <= dCrit);
ind2 = find(D > dCrit);

% figure
scatter(Y(ind1,4),Y(ind1,5),10,'b.'); hold on;
scatter(Y(ind2,4),Y(ind2,5),10,'r.');

figure
line(Y(ind1, 1), Y(ind1, 2), 'Color', 'b', 'LineStyle', '.');
line(Y(ind2, 1), Y(ind2, 2), 'Color', 'r', 'LineStyle', '.');
end