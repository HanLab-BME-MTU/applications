function makeTropoFigure5(paths, outputDirectory)

nBands = 4;
distIntervals = [0 .5 1 2 5 10 15] * 1000;

for iTM = 1:numel(paths)
    % Load Movie Data
    fileName = [paths{iTM} filesep 'movieData.mat'];
    if ~exist(fileName, 'file')
        error(['Unable to locate ' fileName]);
    end
    load(fileName);

    nFrames = movieData.labels.nFrames;
    pixelSize = movieData.pixelSize_nm;
    timeInterval = movieData.timeInterval_s;
    
    % Read distance transforms
    distToEdge = zeros(movieData.imSize(1),movieData.imSize(2),nFrames);
    
    for iFrame = 1:nFrames
        fileName = fullfile(bwdistPath, bwdistFiles(iFrame).name);
        tmp = load(fileName);
        distToEdge(:,:,iFrame) = tmp.distToEdge;
    end
    
    % Read the MPMs
    load(fullfile(movieData.fsmDirectory{1}, 'tack', 'mpm.mat'));
    mpmTM = MPM;
    load(fullfile(movieData.fsmDirectory{1}, 'tack', 'mpm.mat'));
    mpmActin = MPM;

    for iBand = 1:nBands
        % Get the lifetime and average velocity of TM tracks that belong to
        % the ith band.
        [lifeTimeTM avgSpeedTM] = getTrackParamsFromMPM(mpmTM,distToEdge,...
            distIntervals(iBand:iBand+1) / pixelSize, 2, true);

        % Get the lifetime and average velocity of Actin tracks that belong
        % to the ith band.
        [lifeTimeActin avgSpeedActin] = getTrackParamsFromMPM(mpmActin,...
            distToEdge,distIntervals(iBand:iBand+1) / pixelSize, 2, true);
    end
end

colors = [
   0.983333333333333   1.000000000000000   0.800000000000000;
   0.360000000000000   0.630000000000000   0.900000000000000;
   0.060000000000000   0.330000000000000   0.600000000000000;
   0.700000000000000   0.245000000000000   0.245000000000000;
   0.550000000000000                   0                   0;
   0.250000000000000                   0                   0; ];

