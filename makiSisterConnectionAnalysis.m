function analysisStruct = makiSisterConnectionAnalysis(jobType,analysisStruct,verbose)
%MAKISISTERCONNECTIONANALYSIS analyzes in multiple ways the connection between sisters
%
%SYNOPSIS analysisStruct = makiSisterConnectionAnalysis(jobType,analysisStruct,verbose)
%
%INPUT  jobType: 1: test job
%                2: hercules run (default)
%       analysisStruct: Structure with field movies indicating the movies
%                       to be analyzed. Optional. If not input, GUI to load
%                       movies is launched.
%       verbose: 1 to make plots, 0 otherwise. Optional. Default: 0.
%
%OUTPUT analysisStruct: Same as input but with addition field
%           .sisterConnection: 
%
%Khuloud Jaqaman, July 2007

%% input
if nargin < 1 || isempty(jobType)
    jobType = 2;
end

if nargin < 3 || isempty(verbose)
    verbose = 0;
end

%interactively obtain analysisStruct if not input
if nargin < 2 || isempty(analysisStruct) || ~isfield(analysisStruct,'movies')
    analysisStruct = makiCollectMovies(jobType);
end

%extract fileName and directory name from analysisStruct
fileName = analysisStruct.fileName;
dir2SaveRes = analysisStruct.filePath;

%load dataStruct's belonging to the movies in analysisStruct
moviesList = analysisStruct.movies;
numMovies = size(moviesList,1);
for iMovie = numMovies : -1 : 1
    dataStruct(iMovie,1) = makiLoadDataFile(fullfile(moviesList{iMovie,2},moviesList{iMovie,1}));
end

%find number of sisters in each movie
numSisters = zeros(numMovies,1);
for iMovie = 1 : numMovies
    numSisters(iMovie) = length(dataStruct(iMovie).sisterList);
    
    % when the length of sisterList == 1, it could either mean that there
    % is no sister at all, or that there is only 1 sister. Check here which
    % one is the case
    if numSisters(iMovie) == 1
        if isempty(dataStruct(iMovie).sisterList.distances)
            % no sister
            numSisters(iMovie) = 0;
            if verbose 
                disp(sprintf('no sisters in %s',fullfile(moviesList{iMovie,2},moviesList{iMovie,1})));
            end
        end
    end
end
numSistersTot = sum(numSisters);

%% collect angles and distances

%initialize global sister index and structures
iGlobal = 0;
sisterDist(1:numSistersTot,1) = struct('observations',[]);
sisterVel(1:numSistersTot,1) = struct('observations',[]);
angleNormal(1:numSistersTot,1) = struct('observations',[]);
angularVel(1:numSistersTot,1) = struct('observations',[]);

%go over all movies
for iMovie = 1 : numMovies

    %copy fields out of dataStruct(iMovie)
    sisterList = dataStruct(iMovie).sisterList;
    tracks = dataStruct(iMovie).tracks;
    frameAlignment = dataStruct(iMovie).frameAlignment;
    planeFit = dataStruct(iMovie).planeFit;
    timeLapse = dataStruct(iMovie).dataProperties.timeLapse;

    %determine frames where there is a plane
    phase = vertcat(planeFit.phase)';
    framesWithPlane = find(phase ~= 'e');

    %go over all sisters in movie
    if numSisters(iMovie) > 0
        for iSister = 1 : numSisters(iMovie)

            %increase global index of sisters by 1
            iGlobal = iGlobal + 1;

            %find track indices
            tracksIndx = sisterList(1).trackPairs(iSister,1:2);

            %determine frame where each track starts
            trackStart = [tracks(tracksIndx(1)).seqOfEvents(1,1) ...
                tracks(tracksIndx(2)).seqOfEvents(1,1)];

            %find number of frames and frames where pair "exists"
            goodFrames = ~isnan(sisterList(iSister).distances(:,1));
            numFrames = length(goodFrames);
            goodFrames = find(goodFrames);

            %find feature indices making up sisters
            sisterIndx1 = NaN(numFrames,1);
            sisterIndx2 = NaN(numFrames,1);
            sisterIndx1(goodFrames) = tracks(tracksIndx(1)).tracksFeatIndxCG(goodFrames-trackStart(1)+1);
            sisterIndx2(goodFrames) = tracks(tracksIndx(2)).tracksFeatIndxCG(goodFrames-trackStart(2)+1);

            %get aligned sister coordinates
            coords1 = NaN(numFrames,6);
            coords2 = NaN(numFrames,6);
            for iFrame = goodFrames'
                coords1(iFrame,:) = frameAlignment(iFrame).alignedCoord(sisterIndx1(iFrame),:);
                coords2(iFrame,:) = frameAlignment(iFrame).alignedCoord(sisterIndx2(iFrame),:);
            end

            %calculate vector between sisters
            sisterVec = [coords2(:,1:3)-coords1(:,1:3) sqrt(coords1(:,4:6).^2+coords2(:,4:6).^2)];

            %calculate angle with normal
            angleWithNorm = NaN(numFrames,1);
            for iFrame = framesWithPlane
                angleWithNorm(iFrame) = acos(abs(sisterVec(iFrame,1)./norm(sisterVec(iFrame,1:3))));
            end

            %calculate angle between consecutive frames
            angularDisp = NaN(numFrames-1,1);
            for iFrame = 1 : numFrames - 1
                angularDisp(iFrame) = acos(abs(sisterVec(iFrame,1:3) * sisterVec(iFrame+1,1:3)' ...
                    / norm(sisterVec(iFrame,1:3)) / norm(sisterVec(iFrame+1,1:3))));
            end

            sisterDist(iGlobal).observations = sisterList(iSister).distances; %um
            sisterVel(iGlobal).observations = diff(sisterList(iSister).distances(:,1)) ...
                * 1000 / timeLapse; %nm/s
            angleNormal(iGlobal).observations = angleWithNorm * 180 / pi; %deg
            angularVel(iGlobal).observations = angularDisp * 180 / pi / timeLapse; %deg/s

        end %(for iSister = 1 : numSisters(iMovie))
    end

end %(for iMovie = 1 : numMovies)
    
%% means and stds

%distance
allValues = vertcat(sisterDist.observations);
allValues = allValues(:,1);
sisterDistMeanStd = [nanmean(allValues) nanstd(allValues)];

%positive velocity
allValues = vertcat(sisterVel.observations);
allValues = allValues(allValues(:,1)>0,1);
sisterVelPosMeanStd = [nanmean(allValues) nanstd(allValues)];

%negative velocity
allValues = vertcat(sisterVel.observations);
allValues = abs(allValues(allValues(:,1)<0,1));
sisterVelNegMeanStd = [nanmean(allValues) nanstd(allValues)];

%angle with normal
allValues = vertcat(angleNormal.observations);
allValues = allValues(:,1);
angleNormalMeanStd = [nanmean(allValues) nanstd(allValues)];

%angular velocity
allValues = vertcat(angularVel.observations);
allValues = allValues(:,1);
angularVelMeanStd = [nanmean(allValues) nanstd(allValues)];

%% autocorrelation

%define maximum lag
maxLag = 10;

%distance
sisterDistAutocorr = autoCorr(sisterDist,maxLag);

%velocity
sisterVelAutocorr = autoCorr(sisterVel,maxLag);

%angle with normal
angleNormalAutocorr = autoCorr(angleNormal,maxLag);

%angular velocity
angularVelAutocorr = autoCorr(angularVel,maxLag);

%% ARMA

%define model orders to test
modelOrder = [0 3; 0 3; 0 3];

%call ARMA analysis function
% sisterDistArma = armaxFitKalman(sisterDist,[],modelOrder,'tl');
sisterDistArma = [];

%call ARMA analysis function
% sisterVelArma = armaxFitKalman(sisterVel,[],modelOrder,'tl');
sisterVelArma = [];

%call ARMA analysis function
% angleNormalArma = armaxFitKalman(angleNormal,[],modelOrder,'tl');
angleNormalArma = [];

%call ARMA analysis function
% angularVelArma = armaxFitKalman(angularVel,[],modelOrder,'tl');
angularVelArma = [];

%% angle vs. distance

%collect all angles and distances
allAngles = vertcat(angleNormal.observations);
allDistances = vertcat(sisterDist.observations);
goodIndx = find(~isnan(allAngles(:,1)) & ~isnan(allDistances(:,1)));
allAngles = allAngles(goodIndx,1);
allDistances = allDistances(goodIndx,1);

%for now, only plot at the end

%% output to analysisStruct

%store results in one structure
meanStd = struct('distance',sisterDistMeanStd,'posRateChangeDist',...
    sisterVelPosMeanStd,'negRateChangeDist',sisterVelNegMeanStd,...
    'angleWithNormal',angleNormalMeanStd,'angularVel',angularVelMeanStd);
autocorr = struct('distance',sisterDistAutocorr,'rateChangeDist',...
    sisterVelAutocorr,'angleWithNormal',angleNormalAutocorr,'angularVel',...
    angularVelAutocorr);
arma = struct('distance',sisterDistArma,'rateChangeDist',...
    sisterVelArma,'angleWithNormal',angleNormalArma,'angularVel',...
    angularVelArma);

sisterConnection = struct('meanStd',meanStd,'autocorr',autocorr,...
    'arma',arma);

%check whether current analysisStruct already has the sisterConnection field
fieldExists = isfield(analysisStruct,'sisterConnection');

%store results in analysisStruct
analysisStruct.sisterConnection = sisterConnection;

%if sisterConnection field already existed, add 1 to the version number in
%the file name where analysisStruct will be stored
if fieldExists
    [versionNum,fileBody] = makiGetVersion(fileName);
    fileName = [fileBody '_' num2str(versionNum+1) '.mat'];
    analysisStruct.fileName = fileName;
end

%save analysisStruct
save(fullfile(dir2SaveRes,fileName),'analysisStruct');

%for now, since we cannot run the ARMA stuff at MBL but must run it at
%Scripps, save variables in a .mat file to ship them to Scripps
save(fullfile(dir2SaveRes,['arma_' fileName]),'sisterDist','sisterVel',...
    'angleNormal','angularVel','modelOrder');

%% plots

if verbose

    %get time between frames
    timeLapse = round(dataStruct(1).dataProperties.timeLapse);
    
    %get number of frames in each movie
    numFrames = dataStruct(1).dataProperties.movieSize(end);

    %% distance stuff %%
    
    %open figure and write title
    figure('Name',[fileName ' - distances'],'NumberTitle','off');
    
    %plot a sample of trajectories
    
    %create subplot 1
    subplot(2,1,1);
    hold on;
    
    %put all distances together in one matrix
    distanceMat = [sisterDist.observations];
    distanceMat = distanceMat(:,1:2:end);
    
    %plot distance over time for all sisters
    plot((0:numFrames-1)*timeLapse,distanceMat);

    %set axes limit
    axis([0 (numFrames-1)*timeLapse 0 max(distanceMat(:))+1]);
    
    %write axes labels
    xlabel('Time (s)');
    ylabel('Sister separation (um)');
    
    %write averaging information
    text(timeLapse,max(distanceMat(:))+0.6,sprintf([' Sister separation' ...
        ' (um): %4.2f +- %4.2f \n Rate of change of sister separation ' ...
        '(nm/s): \n +ve: %4.2f +- %4.2f, -ve: %4.2f +- %4.2f'],...
        sisterDistMeanStd(1),sisterDistMeanStd(2),sisterVelPosMeanStd(1),...
        sisterVelPosMeanStd(2),sisterVelNegMeanStd(1),sisterVelNegMeanStd(2)));
    
    %hold off figure
    hold off
    
    %plot autocorrelation functions

    %create subplot 2
    subplot(2,1,2); 
    hold on;

    %plot the distance and velocity autocorrelations
    plot((0:maxLag)*timeLapse,sisterDistAutocorr,'k','marker','.');
    plot((0:maxLag)*timeLapse,sisterVelAutocorr,'r','marker','.');

    %set axes limit
    axis([0 maxLag*timeLapse min(0,1.1*min([sisterDistAutocorr;sisterVelAutocorr])) 1.1]);
    
    %write axes labels
    xlabel('Lag (s)');
    ylabel('Autocorrelation');

    %write legend
    text(1*timeLapse,0.9,sprintf([' Black: Sister separation \n Red: ' ...
        'Rate of change of sister separation']));

    %hold off figure
    hold off

    %% angle stuff %%
    
    %open figure and write title
    figure('Name',[fileName ' - angles'],'NumberTitle','off');
    
    %plot a sample of time series of angle with normal
    
    %create subplot 1
    subplot(2,2,1);
    hold on;
    
    %put all angles with normal together in one matrix
    angleMat = [angleNormal.observations];
    
    %plot angles with normal over time for all sisters
    plot((0:numFrames-1)*timeLapse,angleMat);

    %set axes limit
    axis([0 (numFrames-1)*timeLapse 0 90]);
    
    %write axes labels
    xlabel('Time (s)');
    ylabel('Angle with normal to metaphase plate (degrees)');
    
    %write averaging information
    text(timeLapse,80,sprintf('angle (degrees): %4.2f +- %4.2f',...
        angleNormalMeanStd(1),angleNormalMeanStd(2)));
    
    %hold off figure
    hold off
    
    %plot autocorrelation function of angle with normal

    %create subplot 3
    subplot(2,2,3); 
    hold on;

    %plot the autocorrelation of angle with normal
    plot((0:maxLag)*timeLapse,angleNormalAutocorr,'k','marker','.');

    %set axes limit
    axis([0 maxLag*timeLapse min(0,1.1*min(angleNormalAutocorr)) 1.1]);
    
    %write axes labels
    xlabel('Lag (s)');
    ylabel('Autocorrelation of angle with normal');

    %hold off figure
    hold off

    %plot a sample of time series of angular velocity
    
    %create subplot 2
    subplot(2,2,2);
    hold on;
    
    %put all angular velocities together in one matrix
    angleMat = [angularVel.observations];
    
    %plot angular velocities over time for all sisters
    plot((0:numFrames-2)*timeLapse,angleMat);

    %set axes limit
    axis([0 (numFrames-2)*timeLapse 0 1.1*max(angleMat(:))]);
    
    %write axes labels
    xlabel('Time (s)');
    ylabel('Angular velocity (degrees/s)');
    
    %write averaging information
    text(timeLapse,1.05*max(angleMat(:)),sprintf('angular velocity (degrees/s): %4.2f +- %4.2f',...
        angularVelMeanStd(1),angularVelMeanStd(2)));
    
    %hold off figure
    hold off
    
    %plot autocorrelation function of angle with normal

    %create subplot 4
    subplot(2,2,4); 
    hold on;

    %plot the autocorrelation of angular velocity
    plot((0:maxLag)*timeLapse,angularVelAutocorr,'k','marker','.');

    %set axes limit
    axis([0 maxLag*timeLapse min(0,1.1*min(angularVelAutocorr)) 1.1]);
    
    %write axes labels
    xlabel('Lag (s)');
    ylabel('Autocorrelation of angular velocity');

    %hold off figure
    hold off

    %% angle vs. distance %%

    %open figure and write title
    figure('Name',[fileName ' - angle vs. distance'],'NumberTitle','off');
    
    %plot angle vs. distance as a scatter plot
    plot(allDistances,allAngles,'k.')
        
    %write axes labels
    xlabel('Sister separation (um)');
    ylabel('Angle with normal (degrees)');
    
end

%% ~~~ the end ~~~ %%

