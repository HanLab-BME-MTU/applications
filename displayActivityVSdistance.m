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
activityVSdistance(1).distance = activityVSdistance(1).distance * ps;
activityVSdistance(2).distance = activityVSdistance(2).distance * ps;
% Convert into nm/s
activityVSdistance(1).activity = activityVSdistance(1).activity * ps / t;
activityVSdistance(2).activity = activityVSdistance(2).activity * ps / t;

ind1 = find(activityVSdistance(1).distance >= loDist & ...
    activityVSdistance(1).distance <= hiDist);
ind2 = find(activityVSdistance(2).distance >= loDist & ...
    activityVSdistance(2).distance <= hiDist);

figure('Name', ['Distance of Speckle to the Cell edge in Function of Cell Activity  [' ...
    num2str(loDist * ps) '-' num2str(hiDist * ps) ']um']);

plot(activityVSdistance(1).activity(ind1), activityVSdistance(1).distance(ind1), 'r.');
hold on;
plot(activityVSdistance(2).activity(ind2), activityVSdistance(2).distance(ind2), 'g.');
xlabel('protrusion/retraction speed (nm/s)');
ylabel('distance to the edge (nm)');
legend(activityVSdistance(1).name, activityVSdistance(2).name);