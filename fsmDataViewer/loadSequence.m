function newSettings = loadSequence(settings)

newSettings = settings;

h = waitbar(0, 'Load sequence...');

numFrames = settings.numFrames;

% Load the first image
I = loadImage(settings, 1);

waitbar(1 / numFrames);

% Allocate the sequence
[m n c] = size(I);
newSettings.sequence = zeros([m n c numFrames], class(I));
newSettings.sequence(:, :, :, 1) = I;
clear I;

% Load the rest of the sequence

for iFrame = 2:numFrames
    newSettings.sequence(:, :, :, iFrame) = loadImage(settings, iFrame);
    
    waitbar(iFrame / numFrames);
end

close(h);

end