function analysisStruct = makiSisterMotionCoupling(jobType,analysisStruct,verbose)
%MAKISISTERMOTIONCOUPLING looks for coupling in the motion of sister kinetochores
%
%SYNOPSIS analysisStruct = makiSisterConnectionAnalysis(jobType,analysisStruct,verbose)
%
%INPUT  jobType: 1: test job
%                2: hercules run (default)
%       analysisStruct: Structure with fields
%               .fileName: Name of file where analysisStruct is stored.
%               .filePath: Directory where analysisStruct is stored.
%               .movies  : nx2 cell array of the names of the n selected movies.
%                          First column: file name, second column: file path.
%                       Optional. If not input, GUI to load movies is launched.
%       verbose: 1 to make plots, 0 otherwise. Optional. Default: 0.
%
%OUTPUT analysisStruct: Same as input but with addition field
%           .motionCoupling: 
%
%Khuloud Jaqaman, August 2007

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

%% get aligned sister coordinates

%initialize global sister index and structures
iGlobal = 0;
projectionSis1(1:numSistersTot,1) = struct('observations',[]);
projectionSis2(1:numSistersTot,1) = struct('observations',[]);
angleSis1(1:numSistersTot,1) = struct('observations',[]);
angleSis2(1:numSistersTot,1) = struct('observations',[]);

%go over all movies
for iMovie = 1 : numMovies

    %copy fields out of dataStruct(iMovie)
    sisterList = dataStruct(iMovie).sisterList;
    tracks = dataStruct(iMovie).tracks;
    frameAlignment = dataStruct(iMovie).frameAlignment;

    %get number of sisters in movie
    numSisters = length(sisterList);
    
    %go over all sisters in movie
    if numSisters(iMovie) > 0

        for iSister = 1 : numSisters

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

            %calculate vector between sisters and normalize it
            sisterVec = [coords2(:,1:3)-coords1(:,1:3) sqrt(coords1(:,4:6).^2+coords2(:,4:6).^2)];
            sisterVec = sisterVec ./ repmat(sqrt(sum(sisterVec(:,1:3).^2,2)),1,6);

            %calculate the displacement between frames
            coords1Diff = [coords1(2:end,1:3)-coords1(1:end-1,1:3) ...
                sqrt(coords1(2:end,4:6).^2+coords1(1:end-1,4:6).^2)];
            coords2Diff = [coords2(2:end,1:3)-coords2(1:end-1,1:3) ...
                sqrt(coords2(2:end,4:6).^2+coords2(1:end-1,4:6).^2)];

            %calculate projection of displacements on vector connecting sisters
            projection1 = coords1Diff(:,1:3) * sisterVec(1:end-1,1:3)';
            projection1 = diag(projection1);
            projection2 = coords2Diff(:,1:3) * sisterVec(1:end-1,1:3)';
            projection2 = diag(projection2);

            %calculate angle between displacements and vector connecting sisters
            angle1 = acos(projection1 ./ sqrt(sum(coords1Diff(:,1:3).^2,2))) * 180 / pi;
            angle2 = acos(projection2 ./ sqrt(sum(coords2Diff(:,1:3).^2,2))) * 180 / pi;

            %store projections and angles in their structures
            projectionSis1(iGlobal).observations = projection1; %um
            projectionSis2(iGlobal).observations = projection2; %um
            angleSis1(iGlobal).observations = angle1; %deg
            angleSis2(iGlobal).observations = angle2; %deg

        end %(for iSister = 1 : numSisters)
        
    end
    
end %(for iMovie = 1 : numMovies)
    
%% cross-correlation of motion

%define maximum lag
maxLag = 10;

%projections
projectionCrosscorr = crossCorr(projectionSis1,projectionSis2,maxLag);

%angles
angleCrosscorr = crossCorr(angleSis1,angleSis2,maxLag);

%% ARMAX

%this part will be done later - ARMAX still needs some work

%% output to analysisStruct

%store results in one structure
crosscorr = struct('projections',projectionCrosscorr,'angles',...
    angleCrosscorr);

sisterMotionCoupling = struct('crosscorr',crosscorr);

%check whether current analysisStruct already has the sisterConnection field
fieldExists = isfield(analysisStruct,'sisterMotionCoupling');

%store results in analysisStruct
analysisStruct.sisterMotionCoupling = sisterMotionCoupling;

%if sisterMotionCoupling field already existed, add 1 to the version number in
%the file name where analysisStruct will be stored
if fieldExists
    [versionNum,fileBody] = makiGetVersion(fileName);
    fileName = [fileBody '_' num2str(versionNum+1) '.mat'];
    analysisStruct.fileName = fileName;
end

%save analysisStruct
save(fullfile(dir2SaveRes,fileName),'analysisStruct');

%% plots

if verbose

    %get time between frames
    timeLapse = round(dataStruct(1).dataProperties.timeLapse);
    
    %open figure and write title
    figure('Name',[fileName ' - motion coupling'],'NumberTitle','off');
    
    %create subplot 1
    subplot(2,1,1);
    hold on;
    
    %plot projection crosscorrelation
    plot((-maxLag:maxLag)*timeLapse,projectionCrosscorr,'k','marker','.');

    %set axes limit
    axis([-maxLag*timeLapse maxLag*timeLapse min(0,1.1*min(projectionCrosscorr)) 1.1]);
    
    %write axes labels
    xlabel('Lag (s)');
    ylabel('Projection crosscorrelation');

    %hold off figure
    hold off

    %create subplot 2
    subplot(2,1,2);
    hold on;
    
    %plot angle crosscorrelation
    plot((-maxLag:maxLag)*timeLapse,angleCrosscorr,'k','marker','.');

    %set axes limit
    axis([-maxLag*timeLapse maxLag*timeLapse min(0,1.1*min(projectionCrosscorr)) 1.1]);
    
    %write axes labels
    xlabel('Lag (s)');
    ylabel('Angle crosscorrelation');

    %hold off figure
    hold off
    
end

%% ~~~ the end ~~~ %%

