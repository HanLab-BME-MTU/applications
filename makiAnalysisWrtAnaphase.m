function analysisStruct = makiAnalysisWrtAnaphase(jobType,analysisStruct,verbose)
%MAKIANALYSISWRTANAPHASE analyzes sister behavior just before anaphase onset
%
%SYNOPSIS analysisStruct = makiAnalysisWrtAnaphase(jobType,analysisStruct,verbose)
%
%INPUT  jobType: string which can take the values:
%               'TEST', 'HERCULES', 'DANUSER', 'MERALDI', 'SWEDLOW' or
%               'MCAINSH'
%       analysisStruct: Structure with fields
%               .fileName: Name of file where analysisStruct is stored.
%               .filePath: Directory where analysisStruct is stored.
%               .movies  : nx2 cell array of the names of the n selected movies.
%                          First column: file name, second column: file path.
%                       Optional. If not input, GUI to load movies is launched.
%       verbose: 1 to make plots, 0 otherwise. Optional. Default: 0.
%
%OUTPUT analysisStruct: Same as input but with additional field
%           .wrtAnaphaseOnset:
%
%Khuloud Jaqaman, November 2007

%% input
if nargin < 1 || isempty(jobType)
    jobType = 'DANUSER';
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

%load dataStructs belonging to the movies in analysisStruct
moviesList = analysisStruct.movies;
numMovies = size(moviesList,1);
for iMovie = numMovies : -1 : 1
    dataStruct(iMovie,1) = makiLoadDataFile(jobType,fullfile(moviesList{iMovie,2},moviesList{iMovie,1}));
end

%find anaphase onset for each movie
firstFrameAna = NaN(numMovies,1);
for iMovie = 1 : numMovies
    framePhase = vertcat(dataStruct(iMovie).updatedClass.phase);
    firstFrameAnaTmp = find(framePhase=='a',1,'first');
    if ~isempty(firstFrameAnaTmp)
        firstFrameAna(iMovie) = firstFrameAnaTmp;
    end
end

%keep only movies where there is anaphase and it starts after frame 10
goodMovies = find(firstFrameAna > 10);
dataStruct = dataStruct(goodMovies);
firstFrameAna = firstFrameAna(goodMovies);
numAnaFramesSave = 30;
firstFrameSave = firstFrameAna + numAnaFramesSave;
numMovies = length(dataStruct);

%get number of frames in movies (assume that all movies have the same
%number of frames)
%add 30 to look at the 30 frames after anaphase onset also
numFramesAll = dataStruct(1).dataProperties.movieSize(end) + numAnaFramesSave;

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

%% collect sister information, aligned with anaphase onset

%define cell array of sister labels
label{1,1} = 'Inlier';
label{2,1} = 'Unaligned';
label{3,1} = 'Lagging';

%initialize global sister index and structures
for iLabel = 1 : 3

    eval(['iGlobal' label{iLabel,1} ' = 0;'])

    eval(['sisterDist' label{iLabel,1} ' = NaN(numSistersTot,numFramesAll);'])
    eval(['sisterVel' label{iLabel,1} ' = NaN(numSistersTot,numFramesAll);'])
    eval(['angleNormal' label{iLabel,1} ' = NaN(numSistersTot,numFramesAll);'])
    eval(['angularVel' label{iLabel,1} ' = NaN(numSistersTot,numFramesAll);'])
    eval(['projectionSis1' label{iLabel,1} ' = NaN(numSistersTot,numFramesAll);'])
    eval(['projectionSis2' label{iLabel,1} ' = NaN(numSistersTot,numFramesAll);'])
    eval(['angleSis1' label{iLabel,1} ' = NaN(numSistersTot,numFramesAll);'])
    eval(['angleSis2' label{iLabel,1} ' = NaN(numSistersTot,numFramesAll);'])

end

%go over all movies
for iMovie = 1 : numMovies

    %copy fields out of dataStruct(iMovie)
    sisterList = dataStruct(iMovie).sisterList;
    tracks = dataStruct(iMovie).tracks;
    planeFit = dataStruct(iMovie).planeFit;
    frameAlignment = dataStruct(iMovie).frameAlignment;
    updatedClass = dataStruct(iMovie).updatedClass;
    timeLapse = round(dataStruct(iMovie).dataProperties.timeLapse);

    %determine frames where there is a plane
    framesWithPlane = [];
    for t = 1 : length(planeFit)
        if ~isempty(planeFit(t).planeVectors)
            framesWithPlane = [framesWithPlane t];
        end
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
            
            %consider this pair only if it exists before and after anaphase
            %onset
            if any(goodFrames >= firstFrameAna(iMovie)) && ...
                    any(goodFrames < firstFrameAna(iMovie))

                %find feature indices making up sisters
                sisterIndx1 = NaN(numFrames,1);
                sisterIndx2 = NaN(numFrames,1);
                sisterIndx1(goodFrames) = tracks(tracksIndx(1))...
                    .tracksFeatIndxCG(goodFrames-trackStart(1)+1);
                sisterIndx2(goodFrames) = tracks(tracksIndx(2))...
                    .tracksFeatIndxCG(goodFrames-trackStart(2)+1);

                %get aligned sister coordinates
                coords1 = NaN(numFrames,6);
                coords2 = NaN(numFrames,6);
                for iFrame = goodFrames'
                    coords1(iFrame,:) = frameAlignment(iFrame)...
                        .alignedCoord(sisterIndx1(iFrame),:);
                    coords2(iFrame,:) = frameAlignment(iFrame)...
                        .alignedCoord(sisterIndx2(iFrame),:);
                end

                %calculate vector between sisters
                sisterVec = [coords2(:,1:3)-coords1(:,1:3) ...
                    sqrt(coords1(:,4:6).^2+coords2(:,4:6).^2)];

                %get distance between sisters
                sisterDist = sqrt(sum(sisterVec(:,1:3).^2,2)); %um
                
                %normalize vector between sisters
                sisterVec = sisterVec ./ repmat(sisterDist,1,6);

                %calculate rate of change of sister distance
                sisterVel = diff(sisterDist) * 1000 / timeLapse; %nm/s

                %calculate angle with normal
                angleWithNorm = NaN(numFrames,1);
                for iFrame = framesWithPlane
                    angleWithNorm(iFrame) = acos(abs(sisterVec(iFrame,1))); %radians
                end
                angleWithNorm = angleWithNorm * 180 / pi; %degrees

                %calculate angle between consecutive frames
                angularDisp = NaN(numFrames-1,1);
                for iFrame = 1 : numFrames - 1
                    angularDisp(iFrame) = acos(abs(sisterVec(iFrame,1:3) ...
                        * sisterVec(iFrame+1,1:3)')); %radians
                end
                angularDisp = angularDisp * 180 / pi / timeLapse; %degrees

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

                %store information based on the sister type
                iLabel = sisterType(iSister) + 1;

                %increase global index of sister type by 1
                eval(['iGlobal' label{iLabel,1} ' = iGlobal' label{iLabel,1} ' + 1;'])

                numFrames2SaveFrame = firstFrameSave(iMovie) - goodFrames(end);
                if numFrames2SaveFrame >= 0
                    lastFrameStoreSis = goodFrames(end);
                    firstFrameStoreSis = 1;
                    lastFrameStoreGlob = numFramesAll - numFrames2SaveFrame;
                    firstFrameStoreGlob = lastFrameStoreGlob - lastFrameStoreSis + 1;
                else
                    lastFrameStoreSis = goodFrames(end) + numFrames2SaveFrame;
                    firstFrameStoreSis = 1;
                    lastFrameStoreGlob = numFramesAll;
                    firstFrameStoreGlob = lastFrameStoreGlob - lastFrameStoreSis + 1;
                end

                %store information
                eval(['sisterDist' label{iLabel,1} '(iGlobal' label{iLabel,1} ...
                    ',firstFrameStoreGlob:lastFrameStoreGlob) '...
                    '= sisterDist(firstFrameStoreSis:lastFrameStoreSis);']) %um
                eval(['sisterVel' label{iLabel,1} '(iGlobal' label{iLabel,1} ...
                    ',firstFrameStoreGlob:lastFrameStoreGlob-1) '...
                    '= sisterVel(firstFrameStoreSis:lastFrameStoreSis-1);']) %nm/s
                eval(['angleNormal' label{iLabel,1} '(iGlobal' label{iLabel,1} ...
                    ',firstFrameStoreGlob:lastFrameStoreGlob) '...
                    '= angleWithNorm(firstFrameStoreSis:lastFrameStoreSis);']) %deg
                eval(['angularVel' label{iLabel,1} '(iGlobal' label{iLabel,1} ...
                    ',firstFrameStoreGlob:lastFrameStoreGlob-1) '...
                    '= angularDisp(firstFrameStoreSis:lastFrameStoreSis-1);']) %deg/s
                eval(['projectionSis1' label{iLabel,1} '(iGlobal' label{iLabel,1} ...
                    ',firstFrameStoreGlob:lastFrameStoreGlob-1) '...
                    '= projection1(firstFrameStoreSis:lastFrameStoreSis-1);']) %um
                eval(['projectionSis2' label{iLabel,1} '(iGlobal' label{iLabel,1} ...
                    ',firstFrameStoreGlob:lastFrameStoreGlob-1) '...
                    '= projection2(firstFrameStoreSis:lastFrameStoreSis-1);']) %um
                eval(['angleSis1' label{iLabel,1} '(iGlobal' label{iLabel,1} ...
                    ',firstFrameStoreGlob:lastFrameStoreGlob-1) '...
                    '= angle1(firstFrameStoreSis:lastFrameStoreSis-1);']) %deg
                eval(['angleSis2' label{iLabel,1} '(iGlobal' label{iLabel,1} ...
                    ',firstFrameStoreGlob:lastFrameStoreGlob-1) '...
                    '= angle2(firstFrameStoreSis:lastFrameStoreSis-1);']) %deg
                
            end %(if any(goodFrames >= firstFrameAna(iMovie)) && ...)

        end %(for iSister = 1 : numSisters(iMovie) )

    end %(if numSisters(iMovie) > 0)

end %(for iMovie = 1 : numMovies)

%store number of sisters per category
for i = 1 : 3
    eval(['label{i,2} = iGlobal' label{i,1} ' ~= 0;']);
end
goodLabel = find(vertcat(label{:,2}))';

%remove unused rows from matrices
for iLabel = 1 : 3
    eval(['sisterDist' label{iLabel,1} ' = sisterDist' label{iLabel,1} ...
        '(1:iGlobal' label{iLabel,1} ',:);'])
    eval(['sisterVel' label{iLabel,1} ' = sisterVel' label{iLabel,1} ...
        '(1:iGlobal' label{iLabel,1} ',:);'])
    eval(['angleNormal' label{iLabel,1} ' = angleNormal' label{iLabel,1} ...
        '(1:iGlobal' label{iLabel,1} ',:);'])
    eval(['angularVel' label{iLabel,1} ' = angularVel' label{iLabel,1} ...
        '(1:iGlobal' label{iLabel,1} ',:);'])
    eval(['projectionSis1' label{iLabel,1} ' = projectionSis1' label{iLabel,1} ...
        '(1:iGlobal' label{iLabel,1} ',:);'])
    eval(['projectionSis2' label{iLabel,1} ' = projectionSis2' label{iLabel,1} ...
        '(1:iGlobal' label{iLabel,1} ',:);'])
    eval(['angleSis1' label{iLabel,1} ' = angleSis1' label{iLabel,1} ...
        '(1:iGlobal' label{iLabel,1} ',:);'])
    eval(['angleSis2' label{iLabel,1} ' = angleSis2' label{iLabel,1} ...
        '(1:iGlobal' label{iLabel,1} ',:);'])
end

%% mean and std of distances and angles

%initialization
for iLabel = 1 : 3
    eval(['distanceMeanStd' label{iLabel,1} ' = NaN(2,numFramesAll);'])
    eval(['posVelMeanStd' label{iLabel,1} ' = NaN(2,numFramesAll);'])
    eval(['negVelMeanStd' label{iLabel,1} ' = NaN(2,numFramesAll);'])
    eval(['angleNormalMeanStd' label{iLabel,1} ' = NaN(2,numFramesAll);'])
    eval(['angularVelMeanStd' label{iLabel,1} ' = NaN(2,numFramesAll);'])
end

%calculation
for iLabel = goodLabel

    %go over all time intervals from anaphase onset
    for iFrame = 1 : numFramesAll
        
        %distance
        eval(['values = sisterDist' label{iLabel,1} '(:,iFrame);'])
        values = values(~isnan(values));
        if length(values) >= 15
            [valuesMean,valuesStd] = robustMean(values);
            if isinf(valuesStd)
                valuesMean = nanmean(values);
                valuesStd = nanstd(values);
            end
        else
            valuesMean = nanmean(values);
            valuesStd = nanstd(values);
        end
        eval(['distanceMeanStd' label{iLabel,1} ...
            '(:,iFrame) = [valuesMean valuesStd]'';'])

        %positive velocity
        eval(['values = sisterVel' label{iLabel,1} '(:,iFrame);'])
        values = values(values >= 0);
        if length(values) >= 15
            [valuesMean,valuesStd] = robustMean(values);
            if isinf(valuesStd)
                valuesMean = nanmean(values);
                valuesStd = nanstd(values);
            end
        else
            valuesMean = nanmean(values);
            valuesStd = nanstd(values);
        end
        eval(['posVelMeanStd' label{iLabel,1} ...
            '(:,iFrame) = [valuesMean valuesStd]'';'])

        %negative velocity
        eval(['values = sisterVel' label{iLabel,1} '(:,iFrame);'])
        values = values(values < 0);
        if length(values) >= 15
            [valuesMean,valuesStd] = robustMean(values);
            if isinf(valuesStd)
                valuesMean = nanmean(values);
                valuesStd = nanstd(values);
            end
        else
            valuesMean = nanmean(values);
            valuesStd = nanstd(values);
        end
        eval(['negVelMeanStd' label{iLabel,1} ...
            '(:,iFrame) = [valuesMean valuesStd]'';'])

        %angle with normal
        eval(['values = angleNormal' label{iLabel,1} '(:,iFrame);'])
        values = values(~isnan(values));
        if length(values) >= 15
            [valuesMean,valuesStd] = robustMean(values);
            if isinf(valuesStd)
                valuesMean = nanmean(values);
                valuesStd = nanstd(values);
            end
        else
            valuesMean = nanmean(values);
            valuesStd = nanstd(values);
        end
        eval(['angleNormalMeanStd' label{iLabel,1} ...
            '(:,iFrame) = [valuesMean valuesStd]'';'])

        %angular velocity
        eval(['values = angularVel' label{iLabel,1} '(:,iFrame);'])
        values = values(~isnan(values));
        if length(values) >= 15
            [valuesMean,valuesStd] = robustMean(values);
            if isinf(valuesStd)
                valuesMean = nanmean(values);
                valuesStd = nanstd(values);
            end
        else
            valuesMean = nanmean(values);
            valuesStd = nanstd(values);
        end
        eval(['angularVelMeanStd' label{iLabel,1} ...
            '(:,iFrame) = [valuesMean valuesStd]'';'])

    end

end

%% cross-correlation of projections and angles

%initialization
for iLabel = 1 : 3
    eval(['projectionCrosscorr' label{iLabel,1} ' = NaN(1,numFramesAll);'])
    eval(['angleCrosscorr' label{iLabel,1} ' = NaN(1,numFramesAll);'])
end

%calculation
for iLabel = goodLabel

    %go over all time intervals from anaphase onset
    for iFrame = 1 : numFramesAll

        %get cross-correlation of projections for this interval
        eval(['values1 = projectionSis1' label{iLabel,1} '(:,iFrame);'])
        eval(['values2 = projectionSis2' label{iLabel,1} '(:,iFrame);'])
        crosscorrTmp = nancov(values1,values2);
        crosscorrTmp = crosscorrTmp(1,2) / sqrt(crosscorrTmp(1,1) ...
            * crosscorrTmp(2,2));
        eval(['projectionCrosscorr' label{iLabel,1} '(iFrame) = crosscorrTmp;'])

        %get cross-correlation of angles for this interval
        eval(['values1 = angleSis1' label{iLabel,1} '(:,iFrame);'])
        eval(['values2 = angleSis2' label{iLabel,1} '(:,iFrame);'])
        crosscorrTmp = nancov(values1,values2);
        crosscorrTmp = crosscorrTmp(1,2) / sqrt(crosscorrTmp(1,1) ...
            * crosscorrTmp(2,2));
        eval(['angleCrosscorr' label{iLabel,1} '(iFrame) = crosscorrTmp;'])

    end

end

%% output to analysisStruct

for iLabel = 1 : 3

    %store results in one structure
    eval(['distribution = struct(''distance'',sisterDist' label{iLabel,1} ','...
        '''rateChangeDist'',sisterVel' label{iLabel,1} ','...
        '''angleWithNormal'',angleNormal' label{iLabel,1} ','...
        '''angularVel'',angularVel' label{iLabel,1} ');']);
    eval(['meanStd = struct(''distance'',distanceMeanStd' label{iLabel,1} ','...
        '''posRateChangeDist'',posVelMeanStd' label{iLabel,1} ','...
        '''negRateChangeDist'',negVelMeanStd' label{iLabel,1} ','...
        '''angleWithNormal'',angleNormalMeanStd' label{iLabel,1} ','...
        '''angularVel'',angularVelMeanStd' label{iLabel,1} ');']);
    eval(['crosscorr = struct(''projections'',projectionCrosscorr' label{iLabel,1} ','...
        '''angles'',angleCrosscorr' label{iLabel,1} ');']);
    eval(['numSistersCat = iGlobal' label{iLabel,1} ';']);

    eval([label{iLabel,1} ' = struct('...
        '''numSisters'',numSistersCat,'...
        '''distribution'',distribution,'...
        '''meanStd'',meanStd,'...
        '''crosscorr'',crosscorr);'])

end

wrtAnaphaseOnset = struct('Inlier',Inlier,'Unaligned',Unaligned,...
    'Lagging',Lagging);

%check whether current analysisStruct already has the sisterConnection field
fieldExists = isfield(analysisStruct,'wrtAnaphaseOnset');

%store results in analysisStruct
analysisStruct.wrtAnaphaseOnset = wrtAnaphaseOnset;

%if wrtAnaphaseOnset field already existed, add 1 to the version number in
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

numFrames = dataStruct(1).dataProperties.movieSize(end);

if verbose

    for iLabel = goodLabel

        %% distance stuff %%

        %open figure and write title
        figure('Name',[fileName ' - Sister distance vs. anaphase onset - ' ...
            label{iLabel,1}],'NumberTitle','off');

        %create subplot 1
        subplot(2,1,1);
        hold on;

        %plot mean and std of distance over time
        eval(['distanceMat = distanceMeanStd' label{iLabel,1} ';'])
        plot((-numFrames+1:numAnaFramesSave)*timeLapse,distanceMat(1,:),...
            'k','marker','.');
        plot((-numFrames+1:numAnaFramesSave)*timeLapse,distanceMat(2,:),...
            'k:','marker','.');
        
        %plot line highlighting anaphase onset
        plot([0 0],[0 max(distanceMat(:))+1],'r--');

        %set axes limit
        axis([-(numFrames-1)*timeLapse numAnaFramesSave*timeLapse ...
            0 max(distanceMat(:))+1]);

        %write axes labels
        xlabel('Time wrt anaphase onset (s)');
        ylabel('Sister separation (\mum)');

        %write legend
        text((-numFrames+2)*timeLapse,0.8*max(distanceMat(:)),...
            sprintf(' Solid: Mean \n Dashed: Std'));

        %hold off figure
        hold off

        %create subplot 2
        subplot(2,1,2);
        hold on;

        %plot mean and std pos. velocity over time
        eval(['velocityMat = posVelMeanStd' label{iLabel,1} ';'])
        plot((-numFrames+1:numAnaFramesSave)*timeLapse,velocityMat(1,:),...
            'k','marker','+');
        plot((-numFrames+1:numAnaFramesSave)*timeLapse,velocityMat(2,:),...
            'k:','marker','+');
        maxVel1 = max(velocityMat(:));

        %plot mean and std of neg. velocity over time
        eval(['velocityMat = negVelMeanStd' label{iLabel,1} ';'])
        plot((-numFrames+1:numAnaFramesSave)*timeLapse,-velocityMat(1,:),...
            'k','marker','.');
        plot((-numFrames+1:numAnaFramesSave)*timeLapse,velocityMat(2,:),...
            'k:','marker','.');
        maxVel2 = max(velocityMat(:));

        %plot line highlighting anaphase onset
        maxVel = max(maxVel1,maxVel2);
        plot([0 0],[0 maxVel+1],'r--');

        %set axes limit
        axis([-(numFrames-1)*timeLapse numAnaFramesSave*timeLapse ...
            0 maxVel+1]);

        %write axes labels
        xlabel('Time wrt anaphase onset (s)');
        ylabel('Rate change sister separation (nm/s)');

        %write legend
        text((-numFrames+2)*timeLapse,0.8*maxVel,...
            sprintf(' ''+'': Positive rate \n ''.'': Negative rate'));

        %hold off figure
        hold off

        %% angle stuff %%

        %open figure and write title
        figure('Name',[fileName ' - Sister angle vs. anaphase onset - ' ...
            label{iLabel,1}],'NumberTitle','off');

        %create subplot 1
        subplot(2,1,1);
        hold on;

        %plot mean angle with normal over time
        %indicate std by dashed lines
        eval(['angleMat = angleNormalMeanStd' label{iLabel,1} ';'])
        plot((-numFrames+1:numAnaFramesSave)*timeLapse,angleMat(1,:),'k','marker','.');
        plot((-numFrames+1:numAnaFramesSave)*timeLapse,angleMat(2,:),'k:','marker','.');

        %plot line highlighting anaphase onset
        plot([0 0],[0 max(angleMat(:))+1],'r--');

        %set axes limit
        axis([-(numFrames-1)*timeLapse numAnaFramesSave*timeLapse ...
            0 max(angleMat(:))+1]);

        %write axes labels
        xlabel('Time wrt anaphase onset (s)');
        ylabel('Angle with normal to metaphase plate (degrees)');

        %write legend
        text((-numFrames+2)*timeLapse,0.8*max(angleMat(:)),...
            sprintf(' Solid: Mean \n Dashed: Std'));

        %hold off figure
        hold off

        %create subplot 2
        subplot(2,1,2);
        hold on;

        %plot angular velocity over time
        %indicate std by dashed lines
        eval(['angleMat = angularVelMeanStd' label{iLabel,1} ';'])
        plot((-numFrames+1:numAnaFramesSave)*timeLapse,angleMat(1,:),...
            'k','marker','.');
        plot((-numFrames+1:numAnaFramesSave)*timeLapse,angleMat(2,:),...
            'k:','marker','.');

        %plot line highlighting anaphase onset
        plot([0 0],[0 max(angleMat(:))+1],'r--');

        %set axes limit
        axis([-(numFrames-1)*timeLapse numAnaFramesSave*timeLapse ...
            0 max(angleMat(:))+1]);

        %write axes labels
        xlabel('Time wrt anaphase onset (s)');
        ylabel('Sister angular velocity (degrees/s)');

        %write legend
        text((-numFrames+2)*timeLapse,0.8*max(angleMat(:)),...
            sprintf(' Solid: Mean \n Dashed: Std'));

        %hold off figure
        hold off

        %% cross-correlation stuff %%

        %open figure and write title
        figure('Name',[fileName ' - Motion coupling vs. anaphase onset - ' ...
            label{iLabel,1}],'NumberTitle','off');

        %create subplot 1
        subplot(2,1,1);
        hold on;

        %plot projection cross-correlation over time
        eval(['crossCorrTmp = projectionCrosscorr' label{iLabel,1} ';'])
        plot((-numFrames+1:numAnaFramesSave)*timeLapse,crossCorrTmp,'k',...
            'marker','.');

        %plot line highlighting anaphase onset
        plot([0 0],[min(0,1.1*min(crossCorrTmp)) 1.1],'r--');

        %set axes limit
        axis([-(numFrames-1)*timeLapse numAnaFramesSave*timeLapse ...
            min(0,1.1*min(crossCorrTmp)) 1.1]);

        %write axes labels
        xlabel('Time wrt anaphase onset (s)');
        ylabel('Projection cross-correlation');

        %hold off figure
        hold off

        %create subplot 2
        subplot(2,1,2);
        hold on;

        %plot projection cross-correlation over time
        eval(['crossCorrTmp = angleCrosscorr' label{iLabel,1} ';'])
        plot((-numFrames+1:numAnaFramesSave)*timeLapse,crossCorrTmp,'k',...
            'marker','.');

        %plot line highlighting anaphase onset
        plot([0 0],[min(0,1.1*min(crossCorrTmp)) 1.1],'r--');

        %set axes limit
        axis([-(numFrames-1)*timeLapse numAnaFramesSave*timeLapse ...
            min(0,1.1*min(crossCorrTmp)) 1.1]);

        %write axes labels
        xlabel('Time wrt anaphase onset (s)');
        ylabel('Angle cross-correlation');

        %hold off figure
        hold off

    end

end

%% ~~~ the end ~~~ %%

