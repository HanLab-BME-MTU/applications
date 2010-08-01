function makeTropoFigure5(paths, outputDirectory)

nBands = 4;
dLims = [0 .5 1 2 +inf] * 1000;

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
    
    % Read the list of distance transforms
    bwdistPath = movieData.bwdist.directory;
    bwdistFiles = dir([bwdistPath filesep '*.mat']);

    % Read distance transforms
    distToEdge = zeros(movieData.imSize(1),movieData.imSize(2),nFrames);
    
    for iFrame = 1:nFrames
        fileName = fullfile(bwdistPath, bwdistFiles(iFrame).name);
        tmp = load(fileName);
        distToEdge(:,:,iFrame) = tmp.distToEdge * pixelSize;
    end
    
    % Read the TM MPM
    load(fullfile(movieData.fsmDirectory{1}, 'tack', 'mpm.mat'));
    
    % Compute lifetime and average speed of TM tracks within each bands
    % defined by dLims
    [lifeTimeTM, avgSpeedTM] = getLifeTimeAndAvgSpeed(MPM,distToEdge,dLims);
    
    % Read the Actin MPM    
    load(fullfile(movieData.fsmDirectory{2}, 'tack', 'mpm.mat'));
    
    % Compute lifetime and average speed of TM tracks within each bands
    % defined by dLims
    [lifeTimeActin, avgSpeedActin] = getLifeTimeAndAvgSpeed(MPM,distToEdge,dLims);
end

colors = [
   0.983333333333333   1.000000000000000   0.800000000000000;
   0.360000000000000   0.630000000000000   0.900000000000000;
   0.060000000000000   0.330000000000000   0.600000000000000;
   0.700000000000000   0.245000000000000   0.245000000000000;
   0.550000000000000                   0                   0;
   0.250000000000000                   0                   0; ];

