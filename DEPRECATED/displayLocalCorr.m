function R = displayLocalCorr(loDist, hiDist, varargin)
% actinEvents(i, 1) = lifetime of Actin track
% actinEvents(i, 2) = delay from TM birth
% actinEvents(i, 3) = mean speed of Actin track
% actinEvents(i, 4) = mean kinetics score along Actin track

% actinEvents(i, 5:6) = location of the TM birth
% actinEvents(i, 7) = distance of TM birth to the cell edge
% actinEvents(i, 8) = lifetime of TM track
% actinEvents(i, 9) = mean speed of TM track
% actinEvents(i, 10) = mean kinetics score along TM track

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

% Get all subdirectories containing Actin & TM.
paths = getDirectories(rootDirectory);

R = [];

for iMovie = 1:numel(paths)
    fileName = [paths{iMovie} filesep 'localCorrelations.mat'];
    
    if ~exist(fileName, 'file')
        continue;
    end
    
    load(fileName);

    indDist = find(actinEvents(:, 7) >= loDist & actinEvents(:, 7) <= hiDist); %#ok<COLND>
    
%     actinEvents(indDist, 3) = actinEvents(indDist, 3) - min(actinEvents(indDist, 3));
%     actinEvents(indDist, 3) = actinEvents(indDist, 3) / max(actinEvents(indDist, 3));
% 
%     actinEvents(indDist, 9) = actinEvents(indDist, 9) - min(actinEvents(indDist, 9));
%     actinEvents(indDist, 9) = actinEvents(indDist, 9) / max(actinEvents(indDist, 9));
    
    R = [R; mean(actinEvents(indDist, 3)) mean(actinEvents(indDist, 9))]; %#ok<AGROW>
end

R

plot(R(:, 1), R(:, 2), 'g.'); hold on;
% p = polyfit(R(:, 1), R(:, 2), 1);
% plot([0,1], p(1) * [0,1] + p(2), 'r');
% xlabel('Actin speed (nm.min-1)');
% ylabel('TM4 speed (nm.min-1)');
% title(rootDirectory);

% figure;
% 
% titles = {'Actin LifeTime',...
%     'Actin Delay Birth',...
%     'Actin mean speed',...
%     'Actin mean Kin Score',...
%     'distance to Edge',...
%     'TM lifeTime',...
%     'TM mean speed',...
%     'TM mean kinetics'};
% 
% for i = 1:4
%     subplot(2, 4, i);
%     hist(actinEvents(indDist, i));
%     title([titles{i} ' median = ' num2str(median(actinEvents(indDist, i)))]);    
% end
% for i = 5:8
%     subplot(2, 4, i);
%     hist(actinEvents(indDist, i + 2));
%     title([titles{i} ' median = ' num2str(median(actinEvents(indDist, i + 2)))]);
% end

% figure('Name', ['Actin Speed near TM within [' num2str(loDist * ps) '-' num2str(hiDist * ps) ']um']);
% 
% for i = 1:7
%     subplot(2, 4, i);
%     ind = find(actinEvents(indDist, 8) == i);
%     hist(actinEvents(indDist(ind), 3));
%     title(['TM lifetime = ' num2str(i) ' (' num2str(numel(ind)) ')']);
% end
% 
% subplot(2, 4, 8);
% ind = find(actinEvents(indDist, 8) >= 8);
% hist(actinEvents(indDist(ind), 3));
% title(['TM lifetime >= 8 (' num2str(numel(ind)) ')']);
% 
% figure('Name', ['Actin Kinetics near TM within [' num2str(loDist * ps) '-' num2str(hiDist * ps) ']um']);
% for i = 1:7
%     subplot(2, 4, i);
%     ind = find(actinEvents(indDist, 8) == i);
%     hist(actinEvents(indDist(ind), 1));
%     title(['TM lifetime = ' num2str(i) ' (' num2str(numel(ind)) ')']);
% end
% 
% subplot(2, 4, 8);
% ind = find(actinEvents(indDist, 8) >= 8);
% hist(actinEvents(indDist(ind), 1));
% title(['TM lifetime >= 8 (' num2str(numel(ind)) ')']);

end