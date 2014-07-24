function displayGlobalCorr(loDist, hiDist, varargin)

ps = 67 / 1000; % 67nm

if nargin == 3
    rootDirectory = varargin{1};
else
    % Ask for the root directory.
    rootDirectory = uigetdir('', 'Select a root directory:');

    if ~ischar(rootDirectory)
        return;
    end
end

if loDist > hiDist
    error('loDist must be lesser than hiDist');
end

load([rootDirectory filesep 'globalCorrelations.mat']);
indDist1 = find(R{1}(:, 3) >= loDist & R{1}(:, 3) <= hiDist); %#ok<USENS,COLND>
indDist2 = find(R{2}(:, 3) >= loDist & R{2}(:, 3) <= hiDist);

figure('Name', ['Speed Correlation within [' num2str(loDist * ps) '-' num2str(hiDist * ps) ']um']);

for i = 1:7
    subplot(2, 4, i);
    ind1 = find(R{1}(indDist1, 6) == i);
    ind2 = find(R{2}(indDist2, 6) == i);
    plot(R{1}(indDist1(ind1), 4), R{1}(indDist1(ind1), 5), 'r.'); hold on;
    plot(R{2}(indDist2(ind2), 4), R{2}(indDist2(ind2), 5), 'g.');    
    x = 1:max([R{1}(indDist1(ind1), 2); R{1}(indDist1(ind1), 3)]);
    line(x, x, 'Color', 'b');
    title(['lifetime = ' num2str(i) ' (' num2str(numel(ind1)) ')']);
    legend('Actin', 'TM', 'y = x');
end

ind1 = find(R{1}(indDist1, 6) >= 8);
ind2 = find(R{2}(indDist2, 6) >= 8);
subplot(2, 4, 8);
plot(R{1}(indDist1(ind1), 4), R{1}(indDist1(ind1), 5), 'r.'); hold on;
plot(R{2}(indDist2(ind2), 4), R{2}(indDist2(ind2), 5), 'g.');
x = 1:max([R{1}(indDist1(ind1), 2); R{1}(indDist1(ind1), 3)]);
line(x, x, 'Color', 'b');
title(['lifetime >= 8 (' num2str(numel(ind1)) ')']);
legend('Actin', 'TM', 'y = x');
end