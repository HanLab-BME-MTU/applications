function displayTMDelay(varargin)
if nargin > 0
    directory = varargin{1};
else
    % Ask for the root directory.
    directory = uigetdir('', 'Select a directory:');
end

if ~ischar(directory)
    return;
end

% Load the file

fileName = [directory filesep 'tmDelay.mat'];
load(fileName);

%figure('Name', 'Histogram of mean delay between TM and Actin birth');
%hist(tmBirths(:, 4), 50); %#ok<COLND>
%title('Histo of TM delay');
%xlabel('mean TM delay (frame)');

figure('Name', 'Histogram of TM birth distance from the edge (search radius = 300 nm)');

ind = find(tmBirths(:, 4) == -1); %#ok<COLND>
subplot(3, 3, 1); hist(tmBirths(ind, 3) * 76 / 1000, 20);
title('NO Actin track near TM birth');
text(1, 100, ['mode = ' num2str(mode(tmBirths(ind, 3)) * 76 / 1000) ' \mum'], 'Color', 'y');
xlabel('(\mum)');

for delay = 0:5
    ind = find(tmBirths(:, 4) == delay); %#ok<COLND>
    
    subplot(3, 3, delay + 4); hist(tmBirths(ind, 3) * 76 / 1000, 20);
    title(['Mean delay between TM and Actin birth = ' num2str(delay) ' frame.']);
    text(1, 100, ['mode = ' num2str(mode(tmBirths(ind, 3)) * 76 / 1000) ' \mum'], 'Color', 'y');
    xlabel('(\mum)');
end

ind = find(tmBirths(:, 4) > 5); %#ok<COLND>
subplot(3, 3, 3); hist(tmBirths(ind, 3) * 76 / 1000, 20);
title(['Mean delay between TM and Actin birth > ' num2str(delay) ' frame.']);
text(1, 100, ['mode = ' num2str(mode(tmBirths(ind, 3)) * 76 / 1000) ' \mum'], 'Color', 'y');
xlabel('(\mum)');

figure('Name', 'Location of TM birth in function of delay', 'Color', 'k');
set(gca, 'Color', 'k');
color = {'b', 'g', 'r', 'c', 'm', 'y'};
legendNames = cell(6, 1);
for delay = 0:5
    ind = find(tmBirths(:, 4) == delay); %#ok<COLND>
    legendNames{delay+1} = ['Delay == ' num2str(delay) ' (# = ' num2str(numel(ind)) ')'];
    line(tmBirths(ind, 2), tmBirths(ind, 1), 'Color', color{delay+1}, 'LineStyle', '.');
end
axis ij;
axis equal;
legend(legendNames, 'TextColor', 'white');
end