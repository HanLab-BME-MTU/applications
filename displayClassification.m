function displayClassification(loDist, hiDist, pvalue, rootDirectory)

load([rootDirectory filesep 'globalCorrelations.mat']);

% First version: with Mahalanobis distance

% indDist1 = find(R{1}(:, 3) >= 75 & R{1}(:, 3) <= 150 & R{1}(:, 6) > 1); %#ok<USENS,COLND>
% indDist2 = find(R{2}(:, 3) >= 75 & R{2}(:, 3) <= 150 & R{2}(:, 6) > 1);

% X = vertcat(R{1}(indDist1, 4:5), R{2}(indDist2, 4:5)); %#ok<FNDSB>
% obj = gmdistribution(mean(X), cov(X), 1);
% D = mahal(obj,X);
% 
% figure('Name', 'Homogenous population in [5-10]um');
% scatter(X(:,1),X(:,2),10,D(:,1),'.')
% colorbar;
% 
% indDist1 = find(R{1}(:, 3) >= loDist * 1000 / 67 & R{1}(:, 3) <= hiDist * 1000 / 67 & R{1}(:, 6) > 1);
% indDist2 = find(R{2}(:, 3) >= loDist * 1000 / 67 & R{2}(:, 3) <= hiDist * 1000 / 67 & R{2}(:, 6) > 1);
% 
% Y = vertcat(R{1}(indDist1, :), R{2}(indDist2, :)); %#ok<FNDSB>
% D = mahal(obj,Y(:, 4:5));
% 
% figure('Name', ['Population in [' num2str(loDist) '-' num2str(hiDist) ']um']);
% scatter(Y(:,4),Y(:,5),10,D(:,1),'.')
% colorbar;
% 
% dCrit = icdf('chi2', pvalue, 1);
% ind1 = find(D <= dCrit);
% ind2 = find(D > dCrit);
% 
% figure('Name', 'Classification of uncorrelated population');
% scatter(Y(ind1,4),Y(ind1,5),10,'b.'); hold on;
% scatter(Y(ind2,4),Y(ind2,5),10,'r.');
% 
% figure;
% line(Y(ind1, 1), Y(ind1, 2), 'Color', 'k', 'LineStyle', '.');
% line(Y(ind2, 1), Y(ind2, 2), 'Color', 'r', 'LineStyle', '.');
% legend('High speed Actin / low speed TM', 'Correlated Actin/TM speed');

% Second version: with linear regression model

indDist1 = find(R{1}(:, 3) >= 75 & R{1}(:, 3) <= 150 & R{1}(:, 6) > 1); %#ok<USENS,COLND>
indDist2 = find(R{2}(:, 3) >= 75 & R{2}(:, 3) <= 150 & R{2}(:, 6) > 1);

X = vertcat(R{1}(indDist1, 4:5), R{2}(indDist2, 4:5)); %#ok<FNDSB>
% Fit a line to the sample
p = polyfit(X(:, 1), X(:, 2), 1);
% Compute the distance of the sample to y = p(1) x + p(2)
D = (p(1) * X(:, 1) - X(:, 2) + p(2)) / sqrt(p(1)^2 + 1);
% figure('Name', 'Homogenous population in [5-10]um');
% scatter(X(:,1),X(:,2),10,D(:,1),'.')
% colorbar;
% Fit a Gaussian model to D
[mu, sigma] = normfit(D);
% Find the threshold according to the given pvalue
dCrit = icdf('normal', pvalue, mu, sigma);

indDist1 = find(R{1}(:, 3) >= loDist * 1000 / 67 & R{1}(:, 3) <= hiDist * 1000 / 67 & R{1}(:, 6) > 1);
indDist2 = find(R{2}(:, 3) >= loDist * 1000 / 67 & R{2}(:, 3) <= hiDist * 1000 / 67 & R{2}(:, 6) > 1);

Y = vertcat(R{1}(indDist1, :), R{2}(indDist2, :)); %#ok<FNDSB>
% Compute the distance of the sample to y = p(1) x + p(2)
D = (p(1) * Y(:, 4) - Y(:, 5) + p(2)) / sqrt(p(1)^2 + 1);
figure('Name', ['Population in [' num2str(loDist) '-' num2str(hiDist) ']um']);
scatter(Y(:,4),Y(:,5),10,abs(D(:,1)),'.')
% colorbar;

ind1 = find(abs(D) <= dCrit);
ind2 = find(abs(D) > dCrit);

%figure('Name', 'Classification of uncorrelated population');
%scatter(Y(ind1,4),Y(ind1,5),10,'k.'); hold on;
%scatter(Y(ind2,4),Y(ind2,5),10,'r.');

delta = dCrit / cos(atan(1 / p(1)));

line([0, 1200], [delta, 1200 * p(1) + delta], 'LineWidth', 2, 'Color', 'g'); 
line([(delta - 100) / p(1), 1400], [0, 1400 * p(1) - delta + 100], 'LineWidth', 2, 'Color', 'g'); 
 
% imshow(ones(1024, 548));hold on;
% line(Y(ind1, 2), fliplr(Y(ind1, 1)), 'Color', 'k', 'LineStyle', '.');
% line(Y(ind2, 2), fliplr(Y(ind2, 1)), 'Color', 'r', 'LineStyle', '.');
% legend('High speed Actin / low speed TM', 'Correlated Actin/TM speed');

end