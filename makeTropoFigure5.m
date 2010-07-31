function makeTropoFigure5(paths, outputDirectory)

colors = [
   0.983333333333333   1.000000000000000   0.800000000000000;
   0.360000000000000   0.630000000000000   0.900000000000000;
   0.060000000000000   0.330000000000000   0.600000000000000;
   0.700000000000000   0.245000000000000   0.245000000000000;
   0.550000000000000                   0                   0;
   0.250000000000000                   0                   0; ];

names = {'TM2','TM4','TM5NM1'};

for iTM = 1:3
    % Load Movie Data
    fileName = [paths{iTM} filesep 'movieData.mat'];
    if ~exist(fileName, 'file')
        error(['Unable to locate ' fileName]);
    end
    load(fileName);

    nFrames = movieData.labels.nFrames;
    pixelSize = movieData.pixelSize_nm;
    timeInterval = movieData.timeInterval_s;
    
    % Read the MPM
    load(fullfile(movieData.fsmDirectory{1}, 'tack', 'mpm.mat'));
    
    % Correct and use this function (return the 2 speckles classes)
    track_length_plotter(MPM, outputDirectory);    
    
%     nrows = size(MPM,1); %#ok<NODEF>
%     trackID = (1:nrows)';
%     trackMask = MPM(:,1:2:end) ~= 0;
%     
%     % Remove any track that begins at 1st frame
%     startAtFirstFrame = trackMask(:,1);
%     trackMask(:,1) = false;
%     for iFrame = 2:nFrames
%         idxDead = trackID(~trackMask(:,iFrame));
%         trackMask(:,iFrame) = trackMask(:,iFrame) & ~startAtFirstFrame;
%         startAtFirstFrame(idxDead) = false;
%     end
%     % Remove any track that ends at last frame
%     endAtLastFrame = trackMask(:,end);
%     trackMask(:,end) = false;
%     for iFrame = nFrames-1:-1:1
%         idxDead = trackID(~trackMask(:,iFrame));
%         trackMask(:,iFrame) = trackMask(:,iFrame) & ~endAtLastFrame;
%         endAtLastFrame(idxDead) = false;
%     end
%     
%     % Compute pairwise distance
%     D = sqrt((MPM(:,3:2:end-3)-MPM(:,5:2:end-1)).^2 + ...
%         (MPM(:,4:2:end-2)-MPM(:,6:2:end)).^2);
% 
%     % Set pairwise distance of 1 point long track to 0
%     D(trackMask(:, 2:end-1) & ~trackMask(:, 3:end)) = 0;
%     D = [zeros(nrows,1) D zeros(nrows,1)];
% 
%     D(~trackMask) = 0;
%     
%     accuLT = zeros(nrows,1);
%     accuV = zeros(nrows,1);
%     
%     LT = [];
%     V = []; 
%     
%     for iFrame = 1:nFrames
%         idxLive = trackID(trackMask(:,iFrame));
%         idxDead = trackID(~trackMask(:,iFrame));
%         
%         % Accumulate lifetime
%         accuLT(idxLive) = accuLT(idxLive) + 1;
%         
%         % Accumulate velocity
%         accuV(idxLive) = accuV(idxLive) + D(idxLive,iFrame);
%         
%         % Gather lifetime and velocity
%         idx = trackID(~trackMask(:,iFrame) & accuLT > 1);
%         
%         LT = vertcat(LT, (accuLT(idx) - 1) * timeInterval);
%         V = vertcat(V, (accuV(idx) ./ (accuLT(idx) - 1)) * pixelSize * 60 / timeInterval);
%         
%         % Reset
%         accuLT(idxDead) = 0;
%         accuV(idxDead) = 0;
%     end
%     
%     hFig = figure('Visible', 'off');
%     set(gca, 'FontName', 'Helvetica', 'FontSize', 18);
%     set(gcf, 'Position', [680 678 560 400], 'PaperPositionMode', 'auto');
% 
%     scatter(V, LT, 12, colors(2,:), 'Marker', '.');
%     legend(names{iTM}); legend('boxoff');
%     
%     xlabel('Average Edge Velocity (nm/min)');
%     ylabel('Lifetime (s)');
% 
%     fileName = [outputDirectory filesep 'Fig5_A' num2str(iTM) '.eps'];
%     print(hFig, '-depsc', fileName);
%     fixEpsFile(fileName);
%     close(hFig);
end

