function makeTropoFigure6(paths)

paths = paths([3 9 13 15]);

nMovies = numel(paths);

legendNames = {'TM2','TM4','TM5NM1', 'TM5NM1'};

averageDensityTotal = cell(1, 3);

for iMovie = 1:3
    % Load Movie data
    load(fullfile(paths{iMovie}, 'movieData.mat'));
    
    % Load density scores
    load(fullfile(movieData.density.directory, movieData.density.channelDirectory{1}, 'densityScores.mat'));
    
    distFromEdge = vertcat(densityScores(:).distFromEdge) * movieData.pixelSize_nm; %#ok<NODEF>
    averageDensity = vertcat(densityScores(:).averageDensity) * movieData.pixelSize_nm;
    protrusionState = vertcat(densityScores(:).protrusionState);
    %minMaxDensity = vertcat(densityScores(:).minMaxDensity) * movieData.pixelSize_nm;
    
    maxDistFromEdge = min(15000,max(distFromEdge));
    
    dist = 0:500:maxDistFromEdge;
    
    averageDensityTotal{iMovie} = zeros(numel(dist)-1,1);
    
    for i = 1:numel(dist)-1
        averageDensityTotal{iMovie}(i) = mean(averageDensity(distFromEdge > dist(i) & distFromEdge <= dist(i+1) & protrusionState == 3));
    end
end

maxLength = max(cellfun(@numel, averageDensityTotal));

for iMovie = 1:3
    averageDensityTotal{iMovie}(maxLength+1) = 0;
end

averageDensityTotal = horzcat(averageDensityTotal{:});

averageDensityTotal(averageDensityTotal == 0) = NaN;

plot(dist, averageDensityTotal);

legend(legendNames);
xlabel('Distance away from cell edge during retraction (nm)');
ylabel('Average neightbor distance (nm)');

% nMovies = numel(paths);
% 
% scrsz = get(0,'ScreenSize');
% h = figure('Position',scrsz);

% mapNames = {'averageDensityMap*.mat', 'minMaxDensityMap*.mat'};
% nMaps = numel(mapNames);
% movieNames = {'averageDensityMapMovie.mov', 'minMaxDensityMapMovie.mov'};
% titles = {'Average Speckle Distance', 'Min/Max Speckle Distance'};
% 
% for iMovie = 1:nMovies
%     % Load movieData
%     load(fullfile(paths{iMovie}, 'movieData.mat'));
%     
%     if ~checkMovieDensity(movieData)
%         disp('Density needs to be computed!');
%         continue;
%     end
%     
%     nFrames = movieData.labels.nFrames;
%     nChannels = numel(movieData.density.channelDirectory);
%     imSize = movieData.imSize';
%     
%     imageRange = [1 movieData.imSize(1); 1 movieData.imSize(2)];
%     textDeltaCoord = min(diff(imageRange,[],2)) / 20;
%         
%     % Get image file list
%     imagePaths = cellfun(@(x) fullfile(x, 'crop'), movieData.fsmDirectory, 'UniformOutput' ,false);
%     imageFiles = cellfun(@(x) dir([x filesep '*.tif']), imagePaths, 'UniformOutput', false);
% 
%     % Get speckle file list
%     specklePaths = cellfun(@(x) fullfile(x, 'tack', 'locMax'), movieData.fsmDirectory, 'UniformOutput', false);
%     speckleFiles = cellfun(@(x) dir([x filesep '*.mat']), specklePaths, 'UniformOutput', false);
%     
%     % Get map file location
%     mapPaths = cellfun(@(x) fullfile(movieData.density.directory, x), movieData.density.channelDirectory, 'UniformOutput', false);
%     
%     for iMap = 1:nMaps
%         
%         mapFiles = cellfun(@(x) dir([x filesep mapNames{iMap}]), mapPaths, 'UniformOutput', false);
%         
%         % Get the min and max value of the every channel
%         minValue = +inf;
%         maxValue = -inf;
%         
%         for iFrame = 1:nFrames
%             for iChannel = 1:nChannels
%                 map = load(fullfile(mapPaths{iChannel}, mapFiles{iChannel}(iFrame).name));
%                 s = fieldnames(map);
%                 map = map.(s{1});
%                 
%                 minValue = min(minValue, min(nonzeros(map(:))));
%                 maxValue = max(maxValue, max(nonzeros(map(:))));
%             end
%         end
%         
%         clf;
%         
%         MakeQTMovie('start',fullfile(movieData.density.directory, movieNames{iMap}));
%         
%         colorMap = colormap('jet');
%         set(h,'Colormap',colormap('gray'));
%             
%         for iFrame = 1:nFrames
%             for iChannel = 1:nChannels
%                 subplot(1,4, 2*iChannel - 1);
%                                
%                 ima = imread(fullfile(imagePaths{iChannel}, imageFiles{iChannel}(iFrame).name));
%                 
%                 % Display image
%                 imagesc(ima),axis image,axis off;
%                 
%                 hold on;
%                 
%                 % Display speckles
%                 load(fullfile(specklePaths{iChannel}, speckleFiles{iChannel}(iFrame).name));
%                 ind = find(locMax ~= 0);
%                 nPoints = numel(ind);
%                 [Y X] = ind2sub(size(locMax), ind);
%                 plot(X,Y,'r.');
%                 
%                 % Display Time count
%                 text(imageRange(1,1)+textDeltaCoord,imageRange(2,1)+...
%                     textDeltaCoord,num2str(iFrame),'Color','white');
%                 text(imageRange(2,2)-2*textDeltaCoord,imageRange(2,1)+...
%                     textDeltaCoord,movieData.density.channelDirectory{iChannel},...
%                     'Color','white');
%                 
%                 % Display speckle count
%                 text(imageRange(1,1)+textDeltaCoord,imageRange(1,2)-...
%                     textDeltaCoord,['Speckle count = ' num2str(nPoints)],...
%                     'Color','white');
%                 
%                 subplot(1,4,2*iChannel);
%                 
%                 map = load(fullfile(mapPaths{iChannel}, mapFiles{iChannel}(iFrame).name));
%                 s = fieldnames(map);
%                 map = map.(s{1});
%                 
%                 map = (map - minValue) / (maxValue - minValue);
%                 indexedIma = gray2ind(map,size(colorMap,1)) + 1;
%                 coloredIma = zeros([imSize(1) imSize(2),3]);
%                 coloredIma(:,:,1) = reshape(colorMap(indexedIma(:),1), imSize);
%                 coloredIma(:,:,2) = reshape(colorMap(indexedIma(:),2), imSize);
%                 coloredIma(:,:,3) = reshape(colorMap(indexedIma(:),3), imSize);
%                 
%                 imshow(coloredIma);
%                 hold on;
%                 text(imageRange(1,1)+textDeltaCoord,imageRange(2,1)+...
%                     textDeltaCoord,num2str(iFrame),'Color','white');
%                 text(imageRange(2,2)-2*textDeltaCoord,imageRange(2,1)+...
%                     textDeltaCoord,movieData.density.channelDirectory{iChannel},...
%                     'Color','white');
%                 
%                 title(titles{iMap});
%             end
%             
%             %add frame to movie if movie is saved
%             MakeQTMovie('addfigure');
%         end
%         
%         %finish movie
%         MakeQTMovie('finish');
%     
%     end
% end
% 
% close(h);
