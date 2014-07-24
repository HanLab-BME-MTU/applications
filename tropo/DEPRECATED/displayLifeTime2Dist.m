function displayLifeTime2Dist(loDist, hiDist, varargin)

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

TMName = 'TM2';

if ~exist([rootDirectory filesep TMName], 'dir')
    TMName = 'TM4';

    if ~exist([rootDirectory filesep 'TM4'], 'dir')
        error('Unable to find TM directory');
    end
end

fileName = [rootDirectory filesep 'lifeTime2Dist.mat'];

if ~exist(fileName, 'file')
    error('Unable to find result file: redo the analysis.');
end

load(fileName);

figure;

actinInRange = find((actinLifeTime2Dist(:, 1) >= loDist & ...
    actinLifeTime2Dist(:, 1) <= hiDist) | ...
    actinLifeTime2Dist(:, 2) >= loDist & ...
    actinLifeTime2Dist(:, 2) <= hiDist); %#ok<COLND>

tmInRange = find((TMLifeTime2Dist(:, 1) >= loDist & ...
    TMLifeTime2Dist(:, 1) <= hiDist) | ...
    TMLifeTime2Dist(:, 2) >= loDist & ...
    TMLifeTime2Dist(:, 2) <= hiDist); %#ok<COLND>

nActin = hist(actinLifeTime2Dist(actinInRange, 3), 20);
nTM = hist(TMLifeTime2Dist(tmInRange, 3), 20);
bar([nActin / sum(nActin(:)) ; nTM / sum(nTM(:))]', 'group');

title(char(rootDirectory));
set(gcf, 'Name', ['Life Time in range ['...
    num2str(loDist * 76/1000) '-' num2str(hiDist * 76 /1000) '] um']);
legend('Actin Life Time', [TMName ' Life Time']);
xlabel('Life Time (frames)');
ylabel('Frequency');

end