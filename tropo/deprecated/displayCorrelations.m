function displayCorrelations(loDist, hiDist, varargin)

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

ps = movieData.pixelSize_nm / 1000; % 67nm

load([movieData.meanDistances.directory filesep ...
    movieData.meanDistances.filename]);

indDist1 = find(meanDistances{1}(:, 1) >= loDist & meanDistances{1}(:, 1) <= hiDist); %#ok<USENS,COLND>
indDist2 = find(meanDistances{2}(:, 1) >= loDist & meanDistances{2}(:, 1) <= hiDist);

clear meanDistances

figure('Name', ['Speed-Protrusion Correlation Within [' num2str(loDist * ps) '-' num2str(hiDist * ps) ']um']);

load([movieData.meanSpeeds.directory filesep ...
    movieData.meanSpeeds.filename]);

load([movieData.meanProtrusions.directory filesep ...
    movieData.meanProtrusions.filename]);

for i = 1:7
    subplot(2, 4, i);
    ind1 = find(meanSpeeds{1}(indDist1, 3) == i); %#ok<FNDSB,USENS>
    ind2 = find(meanSpeeds{2}(indDist2, 3) == i); %#ok<FNDSB>
    
    speed1 = vertcat(meanSpeeds{1}(ind1, 1), meanSpeeds{2}(ind2, 1));
    speed2 = vertcat(meanSpeeds{1}(ind1, 2), meanSpeeds{2}(ind2, 2));
    protrusion = vertcat(meanProtrusions{1}(ind1, 1), meanProtrusions{2}(ind2, 1)); %#ok<USENS>
    scatter(speed1,speed2,20,protrusion,'.');
    title(['lifetime = ' num2str(i) ' (' num2str(numel(ind1)) ')']);
    set(gca, 'Color', [.3 .3 .3]);
    xlabel('Actin speed (nm.min-1)');
    ylabel('TM speed (nm.min-1)');
    colorbar;
end
% load([rootDirectory filesep 'corr' Param '.mat']);
% figure
% ind = find(corrParam(:, 3) >= loDist & corrParam(:, 3) <= hiDist); %#ok<COLND>
% 
% plot(corrParam(ind, 1), corrParam(ind, 2), 'r.'); %#ok<COLND>
% if strcmp(param, 'speed')
%     xlabel(['Actin ' param, ' (nm/min)']);
%     ylabel([TMName ' ' param, ' (nm/min)']);
% else
%     xlabel(['Actin ' param]);
%     ylabel([TMName ' ' param]);
% end
% set(gcf, 'Name', ['Correlation Actin / ' TMName ' ' Param ' in range ['...
%     num2str(loDist * 76/1000) '-' num2str(hiDist * 76 /1000) '] um']);
% set(gcf, 'NumberTitle', 'off');
% title(char(rootDirectory));
