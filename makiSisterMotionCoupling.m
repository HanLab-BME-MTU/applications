function analysisStruct = makiSisterMotionCoupling(jobType,analysisStruct,verbose)
%MAKISISTERMOTIONCOUPLING looks for coupling in the motion of sister kinetochores
%
%SYNOPSIS analysisStruct = makiSisterMotionCoupling(jobType,analysisStruct,verbose)
%
%INPUT  jobType: string which can take the values:
%               'TEST','HERCULES','DANUSER','MERLADI',SWEDLOW' or MCAINSH
%       analysisStruct: Structure with fields
%               .fileName: Name of file where analysisStruct is stored.
%               .filePath: Directory where analysisStruct is stored.
%               .movies  : nx2 cell array of the names of the n selected movies.
%                          First column: file name, second column: file path.
%                       Optional. If not input, GUI to load movies is launched.
%       verbose: 1 to make plots, 0 otherwise. Optional. Default: 0.
%
%OUTPUT analysisStruct: Same as input but with additional field
%           .motionCoupling: 
%
%Khuloud Jaqaman, August 2007

%% input
if nargin < 1 || isempty(jobType)
    jobType = 'TEST';
end

if nargin < 3 || isempty(verbose)
    verbose = 0;
end

%interactively obtain analysisStruct if not input
if nargin < 2 || isempty(analysisStruct) || ~isfield(analysisStruct,'movies')
    analysisStruct = makiCollectMovies(jobType);
end

%convert file paths in analysisStruct from identifier to real path
analysisStruct = makiMakeAnalysisPlatformIndependent(analysisStruct);

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

%% collect sister projections and angles

%initialize global sister index and structures
iGlobalInlier = 0;
iGlobalUnaligned = 0;
iGlobalLagging = 0;

projectionSis1Inlier(1:numSistersTot,1) = struct('observations',[]);
projectionSis2Inlier(1:numSistersTot,1) = struct('observations',[]);
angleSis1Inlier(1:numSistersTot,1) = struct('observations',[]);
angleSis2Inlier(1:numSistersTot,1) = struct('observations',[]);

projectionSis1Unaligned(1:numSistersTot,1) = struct('observations',[]);
projectionSis2Unaligned(1:numSistersTot,1) = struct('observations',[]);
angleSis1Unaligned(1:numSistersTot,1) = struct('observations',[]);
angleSis2Unaligned(1:numSistersTot,1) = struct('observations',[]);

projectionSis1Lagging(1:numSistersTot,1) = struct('observations',[]);
projectionSis2Lagging(1:numSistersTot,1) = struct('observations',[]);
angleSis1Lagging(1:numSistersTot,1) = struct('observations',[]);
angleSis2Lagging(1:numSistersTot,1) = struct('observations',[]);

%go over all movies
for iMovie = 1 : numMovies

    %copy fields out of dataStruct(iMovie)
    sisterList = dataStruct(iMovie).sisterList;
    tracks = dataStruct(iMovie).tracks;
    frameAlignment = dataStruct(iMovie).frameAlignment;
    updatedClass = dataStruct(iMovie).updatedClass;

    if numSisters(iMovie) > 0

        %sister type = 0 if all kinetochores are inliers
        %sister type = 1 if some kinetochores are unaligned
        %sister type = 2 if some kinetochores are lagging
        sisterType = zeros(numSisters(iMovie),1);
        sisterType(updatedClass(1).sistersUnaligned) = 1;
        sisterType(updatedClass(1).sistersLagging) = 2;
        
        %go over all sisters in movie
        for iSister = 1 : numSisters(iMovie)

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

            %store projections and angles based on the sister type
            switch sisterType(iSister)

                case 0

                    %increase global index of inlier sisters by 1
                    iGlobalInlier = iGlobalInlier + 1;

                    %store information
                    projectionSis1Inlier(iGlobalInlier).observations = projection1; %um
                    projectionSis2Inlier(iGlobalInlier).observations = projection2; %um
                    angleSis1Inlier(iGlobalInlier).observations = angle1; %deg
                    angleSis2Inlier(iGlobalInlier).observations = angle2; %deg

                case 1

                    %increase global index of unaligned sisters by 1
                    iGlobalUnaligned = iGlobalUnaligned + 1;

                    %store information
                    projectionSis1Unaligned(iGlobalUnaligned).observations = projection1; %um
                    projectionSis2Unaligned(iGlobalUnaligned).observations = projection2; %um
                    angleSis1Unaligned(iGlobalUnaligned).observations = angle1; %deg
                    angleSis2Unaligned(iGlobalUnaligned).observations = angle2; %deg

                case 2

                    %increase global index of unaligned sisters by 1
                    iGlobalLagging = iGlobalLagging + 1;

                    %store information
                    projectionSis1Lagging(iGlobalLagging).observations = projection1; %um
                    projectionSis2Lagging(iGlobalLagging).observations = projection2; %um
                    angleSis1Lagging(iGlobalLagging).observations = angle1; %deg
                    angleSis2Lagging(iGlobalLagging).observations = angle2; %deg
                    
            end %(switch sisterType(iSister))

        end %(for iSister = 1 : numSisters(iMovie) )

    end %(if numSisters(iMovie) > 0)

end %(for iMovie = 1 : numMovies)

%define cell array of sister labels
label{1,1} = 'Inlier';
label{2,1} = 'Unaligned';
label{3,1} = 'Lagging';
for i = 1 : 3
    eval(['label{i,2} = iGlobal' label{i,1} ' ~= 0;']);
end
goodLabel = find(vertcat(label{:,2}))';

%remove unused entries from structures
for iLabel = 1 : 3

    eval(['projectionSis1' label{iLabel,1} ' = projectionSis1' label{iLabel,1} ...
        '(1:iGlobal' label{iLabel,1} ');'])
    eval(['projectionSis2' label{iLabel,1} ' = projectionSis2' label{iLabel,1} ...
        '(1:iGlobal' label{iLabel,1} ');'])
    eval(['angleSis1' label{iLabel,1} ' = angleSis1' label{iLabel,1} ...
        '(1:iGlobal' label{iLabel,1} ');'])
    eval(['angleSis2' label{iLabel,1} ' = angleSis2' label{iLabel,1} ...
        '(1:iGlobal' label{iLabel,1} ');'])

end
    
%% cross-correlation of motion

%define maximum lag
maxLag = 10;

%initialization
for iLabel = 1 : 3
    eval(['projectionCrosscorr' label{iLabel,1} ' = [];'])
    eval(['angleCrosscorr' label{iLabel,1} ' = [];'])
end

%calculation
for iLabel = goodLabel

    %projections
    eval(['projectionCrosscorr' label{iLabel,1} ' = crossCorr(projectionSis1' ...
        label{iLabel,1} ',projectionSis2' label{iLabel,1} ',maxLag);'])

    %angles
    eval(['angleCrosscorr' label{iLabel,1} ' = crossCorr(angleSis1' ...
        label{iLabel,1} ',angleSis2' label{iLabel,1} ',maxLag);'])
    
end

%% ARMAX

%this part will be done later - ARMAX still needs some work

%% output to analysisStruct

for iLabel = 1 : 3

    %store results in one structure
    eval(['crosscorr = struct(''projections'',projectionCrosscorr' label{iLabel,1} ','...
        '''angles'',angleCrosscorr' label{iLabel,1} ');']);

    eval([label{iLabel,1} ' = struct(''crosscorr'',crosscorr);'])

end

sisterMotionCoupling = struct('Inlier',Inlier,'Unaligned',Unaligned,...
    'Lagging',Lagging);

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
analysisStruct = makiMakeAnalysisPlatformIndependent(analysisStruct);
save(fullfile(dir2SaveRes,fileName),'analysisStruct');

%% plots

if verbose

    %get time between frames
    timeLapse = round(dataStruct(1).dataProperties.timeLapse);
    
    for iLabel = goodLabel

        %open figure and write title
        figure('Name',[fileName ' - Motion coupling - ' label{iLabel,1}],'NumberTitle','off');

        %create subplot 1
        subplot(2,1,1);
        hold on;

        %plot projection crosscorrelation
        eval(['projectionCrosscorr = projectionCrosscorr' label{iLabel,1} ';']);
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
        eval(['angleCrosscorr = angleCrosscorr' label{iLabel,1} ';']);
        plot((-maxLag:maxLag)*timeLapse,angleCrosscorr,'k','marker','.');

        %set axes limit
        axis([-maxLag*timeLapse maxLag*timeLapse min(0,1.1*min(angleCrosscorr)) 1.1]);

        %write axes labels
        xlabel('Lag (s)');
        ylabel('Angle crosscorrelation');

        %hold off figure
        hold off
        
    end
    
end

%% ~~~ the end ~~~ %%

