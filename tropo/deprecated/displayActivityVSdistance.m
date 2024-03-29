function displayActivityVSdistance(loDist, hiDist, varargin)

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

load([rootDirectory filesep 'windowAnalysis' filesep 'movieData.mat']);

ps = movieData.pixelSize_nm;
t = movieData.timeInterval_s;

load([movieData.activityVSdistance.directory filesep ...
    movieData.activityVSdistance.filename]);

% Convert into nm
activityVSdistance(1).distance = activityVSdistance(1).distance * ps; %#ok<NODEF>
activityVSdistance(2).distance = activityVSdistance(2).distance * ps;
% Convert into nm/s
activityVSdistance(1).activity = activityVSdistance(1).activity * ps / t;
activityVSdistance(2).activity = activityVSdistance(2).activity * ps / t;

ind1 = find(activityVSdistance(1).distance >= loDist & ...
    activityVSdistance(1).distance <= hiDist);
ind2 = find(activityVSdistance(2).distance >= loDist & ...
    activityVSdistance(2).distance <= hiDist);

a1 = activityVSdistance(1).activity(ind1);
a2 = activityVSdistance(2).activity(ind2);
d1 = activityVSdistance(1).distance(ind1);
d2 = activityVSdistance(2).distance(ind2);

minA = min([a1; a2]);
maxA = max([a1; a2]);
h = (maxA - minA) / 6;
g1 = zeros(size(d1));
g2 = zeros(size(d2));
labels = cell(12, 1);

for i = 1:6
    lo = minA + (i-1) * h;
    hi = minA + i * h;
    g1(a1 >= lo & a1 <= hi) = i;
    g2(a2 >= lo & a2 <= hi) = i;
    labels{2 * i - 1} = sprintf('[%.2f', lo);
    labels{2 * i} = sprintf('%.2f]', hi);
end

g1 = 2 * g1 -1;
g2 = 2 * g2;
boxplot([d1; d2], [g1; g2], ...
    'colors', 'rg', 'plotstyle', 'compact', ...
    'labelorientation', 'horizontal', ...
    'labels', labels);
