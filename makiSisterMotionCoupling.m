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
%REMARKS Code is not applicable to anaphase movies/frames
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
analysisStruct = makiMakeAnalysisPlatformIndependent(analysisStruct,jobType);

%extract fileName and directory name from analysisStruct
fileName = analysisStruct.fileName;
dir2SaveRes = analysisStruct.filePath;

%load dataStruct's belonging to the movies in analysisStruct
moviesList = analysisStruct.movies;
numMovies = size(moviesList,1);
for iMovie = numMovies : -1 : 1
    dataStruct(iMovie,1) = makiLoadDataFile(jobType,fullfile(moviesList{iMovie,2},moviesList{iMovie,1}));
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

%define cell array of sister labels
label{1,1} = 'Inlier';
label{2,1} = 'Unaligned';
label{3,1} = 'Lagging';

%initialize global sister index and structures
for iLabel = 1 : 3

    eval(['iGlobal' label{iLabel,1} ' = 0;'])

    eval(['projectionSis1' label{iLabel,1} '(1:numSistersTot,1) = struct(''observations'',[]);'])
    eval(['projectionSis2' label{iLabel,1} '(1:numSistersTot,1) = struct(''observations'',[]);'])
    eval(['angleSis1' label{iLabel,1} '(1:numSistersTot,1) = struct(''observations'',[]);'])
    eval(['angleSis2' label{iLabel,1} '(1:numSistersTot,1) = struct(''observations'',[]);'])
    
end

%go over all movies
for iMovie = 1 : numMovies

    %copy fields out of dataStruct(iMovie)
    sisterList = dataStruct(iMovie).sisterList;
    tracks = dataStruct(iMovie).tracks;
    frameAlignment = dataStruct(iMovie).frameAlignment;
    updatedClass = dataStruct(iMovie).updatedClass;
    numFramesMovie = length(updatedClass);

    %find frame where anaphase starts (if it starts at all)
    framePhase = vertcat(updatedClass.phase);
    firstFrameAna = find(framePhase=='a',1,'first');
    if isempty(firstFrameAna)
        firstFrameAna = numFramesMovie + 1;
    end

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
            goodFrames = goodFrames(goodFrames < firstFrameAna);

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
            iLabel = sisterType(iSister) + 1;

            %increase global index of sister type by 1
            eval(['iGlobal' label{iLabel,1} ' = iGlobal' label{iLabel,1} ' + 1;'])

            %store information
            eval(['projectionSis1' label{iLabel,1} '(iGlobal' label{iLabel,1} ...
                ').observations = projection1;']) %um
            eval(['projectionSis2' label{iLabel,1} '(iGlobal' label{iLabel,1} ...
                ').observations = projection2;']) %um
            eval(['angleSis1' label{iLabel,1} '(iGlobal' label{iLabel,1} ...
                ').observations = angle1;']) %um
            eval(['angleSis2' label{iLabel,1} '(iGlobal' label{iLabel,1} ...
                ').observations = angle2;']) %um
            
        end %(for iSister = 1 : numSisters(iMovie) )

    end %(if numSisters(iMovie) > 0)

end %(for iMovie = 1 : numMovies)

%store number of sisters per category
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

%initialization
for iLabel = 1 : 3
    eval(['projectionCrosscorr' label{iLabel,1} ' = [];'])
    eval(['angleCrosscorr' label{iLabel,1} ' = [];'])
    eval(['projection1Autocorr' label{iLabel,1} ' = [];'])
    eval(['angle1Autocorr' label{iLabel,1} ' = [];'])
    eval(['projection2Autocorr' label{iLabel,1} ' = [];'])
    eval(['angle2Autocorr' label{iLabel,1} ' = [];'])
end

%define maximum lag
maxLag = 10;

%calculation
for iLabel = goodLabel

    %between sisters
    
    %projections
    eval(['projectionCrosscorr' label{iLabel,1} ' = crossCorr(projectionSis1' ...
        label{iLabel,1} ',projectionSis2' label{iLabel,1} ',maxLag);'])

    %angles
    eval(['angleCrosscorr' label{iLabel,1} ' = crossCorr(angleSis1' ...
        label{iLabel,1} ',angleSis2' label{iLabel,1} ',maxLag);'])
    
    %sister with itself
    
    %projections
    eval(['projection1Autocorr' label{iLabel,1} ' = autoCorr(projectionSis1' ...
        label{iLabel,1} ',maxLag);'])

    %angles
    eval(['angle1Autocorr' label{iLabel,1} ' = autoCorr(angleSis1' ...
        label{iLabel,1} ',maxLag);'])
    
    %projections
    eval(['projection2Autocorr' label{iLabel,1} ' = autoCorr(projectionSis2' ...
        label{iLabel,1} ',maxLag);'])

    %angles
    eval(['angle2Autocorr' label{iLabel,1} ' = autoCorr(angleSis2' ...
        label{iLabel,1} ',maxLag);'])
    
end

%% ARMAX

%this part will be done later - ARMAX still needs some work

%% output to analysisStruct

for iLabel = 1 : 3

    %store results in one structure
    eval(['crosscorr = struct(''projections'',projectionCrosscorr' label{iLabel,1} ','...
        '''angles'',angleCrosscorr' label{iLabel,1} ');']);
    eval(['autocorr = struct(''projectionsSis1'',projection1Autocorr' label{iLabel,1} ...
        ',''anglesSis1'',angle1Autocorr' label{iLabel,1} ...
        ',''projectionsSis2'',projection2Autocorr' label{iLabel,1} ...
        ',''anglesSis2'',angle2Autocorr' label{iLabel,1} ');']);
    eval(['numSistersCat = iGlobal' label{iLabel,1} ';']);

    eval([label{iLabel,1} ' = struct('...
        '''numSisters'',numSistersCat,'...
        '''crosscorr'',crosscorr,'...
        '''autocorr'',autocorr);'])

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
analysisStruct = makiMakeAnalysisPlatformIndependent(analysisStruct,jobType);
save(fullfile(dir2SaveRes,fileName),'analysisStruct');

%% plots

if verbose

    %get time between frames
    timeLapse = round(dataStruct(1).dataProperties.timeLapse);
    
    for iLabel = goodLabel

        %open figure and write title
        figure('Name',[fileName(1:end-4) ' - Motion coupling - ' label{iLabel,1}],'NumberTitle','off');

        %create subplot 1
        subplot(2,1,1);
        hold on;

        %plot projection crosscorrelation
        eval(['projectionCrosscorr = projectionCrosscorr' label{iLabel,1} ';']);
        eval(['autocorr1 = projection1Autocorr' label{iLabel,1} ';']);
        eval(['autocorr2 = projection2Autocorr' label{iLabel,1} ';']);
        plot((-maxLag:maxLag)*timeLapse,projectionCrosscorr(:,1),'k','marker','.');
        plot((0:maxLag)*timeLapse,autocorr1(:,1),'g:','marker','.');
        plot((0:maxLag)*timeLapse,autocorr2(:,1),'r:','marker','.');

        %set axes limit
        minVal = min([projectionCrosscorr(:,1); autocorr1(:,1); autocorr2(:,1)]);
        axis([-maxLag*timeLapse maxLag*timeLapse min(0,1.1*minVal) 1.1]);

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
        eval(['autocorr1 = angle1Autocorr' label{iLabel,1} ';']);
        eval(['autocorr2 = angle2Autocorr' label{iLabel,1} ';']);
        plot((-maxLag:maxLag)*timeLapse,angleCrosscorr(:,1),'k','marker','.');
        plot((0:maxLag)*timeLapse,autocorr1(:,1),'g:','marker','.');
        plot((0:maxLag)*timeLapse,autocorr2(:,1),'r:','marker','.');

        %set axes limit
        minVal = min([angleCrosscorr(:,1); autocorr1(:,1); autocorr2(:,1)]);
        axis([-maxLag*timeLapse maxLag*timeLapse min(0,1.1*minVal) 1.1]);

        %write axes labels
        xlabel('Lag (s)');
        ylabel('Angle crosscorrelation');

        %hold off figure
        hold off
        
    end
    
end

%% ~~~ the end ~~~ %%

