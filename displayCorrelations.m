function displayCorrelations(loDist, hiDist, param, varargin)

if nargin == 4
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

if strcmp(param, 'speed')
    Param = 'Speeds';
elseif strcmp(param, 'kinetics')
    Param = 'Kinetics';    
else
    error('Invalid parameters');
end

TMName = 'TM2';

if ~exist([rootDirectory filesep TMName], 'dir')
    TMName = 'TM4';

    if ~exist([rootDirectory filesep 'TM4'], 'dir')
        error('Unable to find TM directory');
    end
end
    
load([rootDirectory filesep 'corr' Param '.mat']);
figure
ind = find(corrParam(:, 3) >= loDist & corrParam(:, 3) <= hiDist); %#ok<COLND>

plot(corrParam(ind, 1), corrParam(ind, 2), 'r.'); %#ok<COLND>
if strcmp(param, 'speed')
    xlabel(['Actin ' param, ' (nm/min)']);
    ylabel([TMName ' ' param, ' (nm/min)']);
else
    xlabel(['Actin ' param]);
    ylabel([TMName ' ' param]);
end
set(gcf, 'Name', ['Correlation Actin / ' TMName ' ' Param ' in range ['...
    num2str(loDist * 76/1000) '-' num2str(hiDist * 76 /1000) '] um']);
set(gcf, 'NumberTitle', 'off');
title(char(rootDirectory));
