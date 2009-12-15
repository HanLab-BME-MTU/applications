function displayActivityVSGeometry(rootDirectory)

load([rootDirectory filesep 'windowAnalysis' filesep 'movieData.mat']);

% Load every windows
load([movieData.windows.directory filesep ...
    movieData.windows.fileName]);

% Get all the windows of the 1st band of every frame.
windows = allWinPoly(1, :, :);
clear allWinPoly;

e = cellfun(@(x) ~(isempty(x) || any(isnan(x(:)))), {windows(:).outerBorder}) & ...
    cellfun(@(x) ~(isempty(x) || any(isnan(x(:)))), {windows(:).innerBorder});

% Compute for every windows the convexity
d2Out = cellfun(@(x) norm(x(:, 1) - x(:, end)), {windows(e).outerBorder});
d2In = cellfun(@(x) norm(x(:, 1) - x(:, end)), {windows(e).innerBorder});

convexity = sqrt(d2Out ./ d2In);

% Load the protrusion samples.
load([movieData.protrusion.directory filesep ...
    movieData.protrusion.samples.fileName]);
            
% Add another frame
avgProt = horzcat(protrusionSamples.averageNormalComponent,...
    protrusionSamples.averageNormalComponent(:, end));

ind = convexity ~= 0;
convexity = convexity(ind);
avgProt = avgProt(ind);

%hist(convexity(avgProt > 0), 50);
%figure, hist(convexity(avgProt < 0), 50);
%figure, plot(avgProt(:), convexity, '.');
if kstest2(convexity(avgProt > 0), convexity(avgProt < 0)) == 0
    disp('OK');
else
    disp('Not OK');
end

end