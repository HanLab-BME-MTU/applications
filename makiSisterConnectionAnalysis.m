function analysisStruct = makiSisterConnectionAnalysis(jobType,analysisStruct,verbose)
%MAKISISTERCONNECTIONANALYSIS analyzes in multiple ways the connection between sisters
%
%SYNOPSIS analysisStruct = makiSisterConnectionAnalysis(jobType,analysisStruct,verbose)
%
%INPUT  jobType: string which can take the values:
%               'TEST','HERCULES','DANUSER','MERLADI',SWEDLOW' or MCAINSH
%       analysisStruct: Structure with field movies indicating the movies
%                       to be analyzed. Optional. If not input, GUI to load
%                       movies is launched.
%       verbose: 1 to make plots, 0 otherwise. Optional. Default: 0.
%
%OUTPUT analysisStruct: Same as input but with additional field
%           .sisterConnection: 
%
%Khuloud Jaqaman, July 2007

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
    numSisters(iMovie,1) = length(dataStruct(iMovie).sisterList);
    
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

%get time between frames
timeLapse = round(dataStruct(1).dataProperties.timeLapse);

%% collect sister angles and distances

%initialize global sister index and structures
iGlobalInlier = 0;
iGlobalUnaligned = 0;
iGlobalLagging = 0;

sisterDistInlier(1:numSistersTot,1) = struct('observations',[]);
sisterVelInlier(1:numSistersTot,1) = struct('observations',[]);
angleNormalInlier(1:numSistersTot,1) = struct('observations',[]);
angularVelInlier(1:numSistersTot,1) = struct('observations',[]);

sisterDistUnaligned(1:numSistersTot,1) = struct('observations',[]);
sisterVelUnaligned(1:numSistersTot,1) = struct('observations',[]);
angleNormalUnaligned(1:numSistersTot,1) = struct('observations',[]);
angularVelUnaligned(1:numSistersTot,1) = struct('observations',[]);

sisterDistLagging(1:numSistersTot,1) = struct('observations',[]);
sisterVelLagging(1:numSistersTot,1) = struct('observations',[]);
angleNormalLagging(1:numSistersTot,1) = struct('observations',[]);
angularVelLagging(1:numSistersTot,1) = struct('observations',[]);

%go over all movies
for iMovie = 1 : numMovies

    %copy fields out of dataStruct(iMovie)
    sisterList = dataStruct(iMovie).sisterList;
    tracks = dataStruct(iMovie).tracks;
    frameAlignment = dataStruct(iMovie).frameAlignment;
    planeFit = dataStruct(iMovie).planeFit;
    timeLapse = dataStruct(iMovie).dataProperties.timeLapse;
    updatedClass = dataStruct(iMovie).updatedClass;

    %determine frames where there is a plane
    phase = vertcat(planeFit.phase)';
    framesWithPlane = find(phase ~= 'e');

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

            %store sister information based on the sister type
            switch sisterType(iSister)
                
                case 0

                    %increase global index of inlier sisters by 1
                    iGlobalInlier = iGlobalInlier + 1;

                    %store information
                    sisterDistInlier(iGlobalInlier).observations = sisterList(iSister).distances; %um
                    sisterVelInlier(iGlobalInlier).observations = diff(sisterList(iSister).distances(:,1)) ...
                        * 1000 / timeLapse; %nm/s
                    angleNormalInlier(iGlobalInlier).observations = angleWithNorm * 180 / pi; %deg
                    angularVelInlier(iGlobalInlier).observations = angularDisp * 180 / pi / timeLapse; %deg/s
                    
                case 1
                    
                    %increase global index of unaligned sisters by 1
                    iGlobalUnaligned = iGlobalUnaligned + 1;

                    %store information
                    sisterDistUnaligned(iGlobalUnaligned).observations = sisterList(iSister).distances; %um
                    sisterVelUnaligned(iGlobalUnaligned).observations = diff(sisterList(iSister).distances(:,1)) ...
                        * 1000 / timeLapse; %nm/s
                    angleNormalUnaligned(iGlobalUnaligned).observations = angleWithNorm * 180 / pi; %deg
                    angularVelUnaligned(iGlobalUnaligned).observations = angularDisp * 180 / pi / timeLapse; %deg/s
                    
                case 2
                    
                    %increase global index of lagging sisters by 1
                    iGlobalLagging = iGlobalLagging + 1;

                    %store information
                    sisterDistLagging(iGlobalLagging).observations = sisterList(iSister).distances; %um
                    sisterVelLagging(iGlobalLagging).observations = diff(sisterList(iSister).distances(:,1)) ...
                        * 1000 / timeLapse; %nm/s
                    angleNormalLagging(iGlobalLagging).observations = angleWithNorm * 180 / pi; %deg
                    angularVelLagging(iGlobalLagging).observations = angularDisp * 180 / pi / timeLapse; %deg/s
                    
            end %(switch sisterType(iSister))

        end %(for iSister = 1 : numSisters(iMovie))
        
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

    eval(['sisterDist' label{iLabel,1} ' = sisterDist' label{iLabel,1} ...
        '(1:iGlobal' label{iLabel,1} ');'])
    eval(['sisterVel' label{iLabel,1} ' = sisterVel' label{iLabel,1} ...
        '(1:iGlobal' label{iLabel,1} ');'])
    eval(['angleNormal' label{iLabel,1} ' = angleNormal' label{iLabel,1} ...
        '(1:iGlobal' label{iLabel,1} ');'])
    eval(['angularVel' label{iLabel,1} ' = angularVel' label{iLabel,1} ...
        '(1:iGlobal' label{iLabel,1} ');'])

end

%% distributions and some distribution parameters

%initialization
for iLabel = 1 : 3
    eval(['sisterDistDistr' label{iLabel,1} ' = [];'])
    eval(['sisterDistParam' label{iLabel,1} ' = [];'])
    eval(['sisterVelDistr' label{iLabel,1} ' = [];'])
    eval(['sisterVelPosParam' label{iLabel,1} ' = [];'])
    eval(['sisterVelNegParam' label{iLabel,1} ' = [];'])
    eval(['angleNormalDistr' label{iLabel,1} ' = [];'])
    eval(['angleNormalParam' label{iLabel,1} ' = [];'])
    eval(['angularVelDistr' label{iLabel,1} ' = [];'])
    eval(['angularVelParam' label{iLabel,1} ' = [];'])
end

%calculation
for iLabel = goodLabel

    %distance
    eval(['allValues = vertcat(sisterDist' label{iLabel,1} '.observations);']);
    allValues = allValues(:,1);
    allValues = allValues(~isnan(allValues));
    eval(['sisterDistDistr' label{iLabel,1} ' = allValues;']);
    eval(['sisterDistParam' label{iLabel,1} ...
        ' = [mean(allValues) std(allValues) min(allValues) ' ...
        'prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) '...
        'max(allValues)];']);

    %overall velocity
    eval(['allValues = vertcat(sisterVel' label{iLabel,1} '.observations);']);
    allValues = allValues(:,1);
    allValues = allValues(~isnan(allValues));
    eval(['sisterVelDistr' label{iLabel,1} ' = allValues;']);
    
    %positive velocity
    eval(['allValues = vertcat(sisterVel' label{iLabel,1} '.observations);']);
    allValues = allValues(allValues(:,1)>0,1);
    allValues = allValues(~isnan(allValues));
    eval(['sisterVelPosParam' label{iLabel,1} ...
        ' = [mean(allValues) std(allValues) min(allValues) ' ...
        'prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) '...
        'max(allValues)];']);

    %negative velocity
    eval(['allValues = vertcat(sisterVel' label{iLabel,1} '.observations);']);
    allValues = abs(allValues(allValues(:,1)<0,1));
    allValues = allValues(~isnan(allValues));
    eval(['sisterVelNegParam' label{iLabel,1} ...
        ' = [mean(allValues) std(allValues) min(allValues) ' ...
        'prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) '...
        'max(allValues)];']);

    %angle with normal
    eval(['allValues = vertcat(angleNormal' label{iLabel,1} '.observations);']);
    allValues = allValues(:,1);
    allValues = allValues(~isnan(allValues));
    eval(['angleNormalDistr' label{iLabel,1} ' = allValues;']);
    eval(['angleNormalParam' label{iLabel,1} ...
        ' = [mean(allValues) std(allValues) min(allValues) ' ...
        'prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) '...
        'max(allValues)];']);

    %angular velocity
    eval(['allValues = vertcat(angularVel' label{iLabel,1} '.observations);'])
    allValues = allValues(:,1);
    allValues = allValues(~isnan(allValues));
    eval(['angularVelDistr' label{iLabel,1} ' = allValues;']);
    eval(['angularVelParam' label{iLabel,1} ...
        ' = [mean(allValues) std(allValues) min(allValues) ' ...
        'prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) '...
        'max(allValues)];']);
    
end

%% temporal trends

%initialization
for iLabel = 1 : 3
    eval(['sisterDistTrend' label{iLabel,1} ' = [];'])
    eval(['sisterVelTrend' label{iLabel,1} ' = [];'])
    eval(['angleNormalTrend' label{iLabel,1} ' = [];'])
    eval(['angularVelTrend' label{iLabel,1} ' = [];'])
end

%calculation
for iLabel = goodLabel

    %get number of sisters in this category
    eval(['numSistersCat = iGlobal' label{iLabel,1} ';']);
    
    %distance
    trendTmp = NaN(numSistersCat,3);
    for iSister = 1 : numSistersCat
        eval(['allValues = sisterDist' label{iLabel,1} '(iSister).observations(:,1);'])
        indxAvail = find(~isnan(allValues));
        [lineParam,S] = polyfit(indxAvail*timeLapse,allValues(indxAvail),1);
        varCovMat = (inv(S.R)*inv(S.R)')*S.normr^2/S.df;
        slopeStd = sqrt(varCovMat(1));
        testStat = lineParam(1)/slopeStd;
        if testStat > 0
            pValue = 1 - tcdf(testStat,S.df);
        else
            pValue = tcdf(testStat,S.df);            
        end
        trendTmp(iSister,:) = [lineParam(1) slopeStd pValue];
    end
    eval(['sisterDistTrend' label{iLabel} ' = trendTmp;']);
    
    %absolute velocity
    trendTmp = NaN(numSistersCat,3);
    for iSister = 1 : numSistersCat
        eval(['allValues = sisterVel' label{iLabel,1} '(iSister).observations(:,1);'])
        allValues = abs(allValues);
        indxAvail = find(~isnan(allValues));
        [lineParam,S] = polyfit(indxAvail*timeLapse,allValues(indxAvail),1);
        varCovMat = (inv(S.R)*inv(S.R)')*S.normr^2/S.df;
        slopeStd = sqrt(varCovMat(1));
        testStat = lineParam(1)/slopeStd;
        if testStat > 0
            pValue = 1 - tcdf(testStat,S.df);
        else
            pValue = tcdf(testStat,S.df);            
        end
        trendTmp(iSister,:) = [lineParam(1) slopeStd pValue];
    end
    eval(['sisterVelTrend' label{iLabel} ' = trendTmp;']);
    
    %angle with normal
    trendTmp = NaN(numSistersCat,3);
    for iSister = 1 : numSistersCat
        eval(['allValues = angleNormal' label{iLabel,1} '(iSister).observations(:,1);'])
        indxAvail = find(~isnan(allValues));
        [lineParam,S] = polyfit(indxAvail*timeLapse,allValues(indxAvail),1);
        varCovMat = (inv(S.R)*inv(S.R)')*S.normr^2/S.df;
        slopeStd = sqrt(varCovMat(1));
        testStat = lineParam(1)/slopeStd;
        if testStat > 0
            pValue = 1 - tcdf(testStat,S.df);
        else
            pValue = tcdf(testStat,S.df);            
        end
        trendTmp(iSister,:) = [lineParam(1) slopeStd pValue];
    end
    eval(['angleNormalTrend' label{iLabel} ' = trendTmp;']);
    
    %angular velocity
    trendTmp = NaN(numSistersCat,3);
    for iSister = 1 : numSistersCat
        eval(['allValues = angularVel' label{iLabel,1} '(iSister).observations(:,1);'])
        indxAvail = find(~isnan(allValues));
        [lineParam,S] = polyfit(indxAvail*timeLapse,allValues(indxAvail),1);
        varCovMat = (inv(S.R)*inv(S.R)')*S.normr^2/S.df;
        slopeStd = sqrt(varCovMat(1));
        testStat = lineParam(1)/slopeStd;
        if testStat > 0
            pValue = 1 - tcdf(testStat,S.df);
        else
            pValue = tcdf(testStat,S.df);            
        end
        trendTmp(iSister,:) = [lineParam(1) slopeStd pValue];
    end
    eval(['angularVelTrend' label{iLabel} ' = trendTmp;']);
    
end

%% autocorrelation

%define maximum lag
maxLag = 10;

%initialization
for iLabel = 1 : 3
    eval(['sisterDistAutocorr' label{iLabel,1} ' = [];'])
    eval(['sisterVelAutocorr' label{iLabel,1} ' = [];'])
    eval(['angleNormalAutocorr' label{iLabel,1} ' = [];'])
    eval(['angularVelAutocorr' label{iLabel,1} ' = [];'])
end

%calculation
for iLabel = goodLabel

    %distance
    eval(['sisterDistAutocorr' label{iLabel,1} ...
        ' = autoCorr(sisterDist' label{iLabel,1} ',maxLag);'])

    %velocity
    eval(['sisterVelAutocorr' label{iLabel,1} ...
        ' = autoCorr(sisterVel' label{iLabel,1} ',maxLag);'])

    %angle with normal
    eval(['angleNormalAutocorr' label{iLabel,1} ...
        ' = autoCorr(angleNormal' label{iLabel,1} ',maxLag);'])

    %angular velocity
    eval(['angularVelAutocorr' label{iLabel,1} ...
        ' = autoCorr(angularVel' label{iLabel,1} ',maxLag);'])
    
end

%% ARMA

%define model orders to test
modelOrder = [0 0; 0 0; -1 -1];

%initialization
for iLabel = 1 : 3
    eval(['sisterDistArma' label{iLabel,1} ' = [];'])
    eval(['sisterVelArma' label{iLabel,1} ' = [];'])
    eval(['angleNormalArma' label{iLabel,1} ' = [];'])
    eval(['angularVelArma' label{iLabel,1} ' = [];'])
end

% %calculation
% for iLabel = goodLabel
%
%     %call ARMA analysis function for sister distance
%     eval(['sisterDistArma' label{iLabel,1} ' = armaxFitKalman(sisterDist' ...
%         label{iLabel,1} ',[],modelOrder,''tl'');'])
%
%     %call ARMA analysis function for sister velocity
%     eval(['sisterVelArma' label{iLabel,1} ' = armaxFitKalman(sisterVel' ...
%         label{iLabel,1} ',[],modelOrder,''tl'');'])
%
%     %call ARMA analysis function for angle with normal to plane
%     eval(['angleNormalArma' label{iLabel,1} ' = armaxFitKalman(angleNormal' ...
%         label{iLabel,1} ',[],modelOrder,''tl'');'])
%
%     %call ARMA analysis function for angular velocity
%     eval(['angularVelArma' label{iLabel,1} ' = armaxFitKalman(angularVel' ...
%         label{iLabel,1} ',[],modelOrder,''tl'');'])
%
% end

%% angle vs. distance

%initialization
for iLabel = 1 : 3
    eval(['allAngles' label{iLabel,1} ' = [];'])
    eval(['allDistances' label{iLabel,1} ' = [];'])
end

%calculation
for iLabel = goodLabel

    %collect all angles and distances
    eval(['allAngles = vertcat(angleNormal' label{iLabel,1} '.observations);'])
    eval(['allDistances = vertcat(sisterDist' label{iLabel,1} '.observations);'])
    goodIndx = find(~isnan(allAngles(:,1)) & ~isnan(allDistances(:,1)));
    eval(['allAngles' label{iLabel,1} ' = allAngles(goodIndx,1);'])
    eval(['allDistances' label{iLabel,1} ' = allDistances(goodIndx,1);'])
    
end

%for now, only plot at the end

%% output to analysisStruct

for iLabel = 1 : 3

    %store results in one structure
    eval(['distribution = struct(''distance'',sisterDistDistr' label{iLabel,1} ','...
        '''rateChangeDist'',sisterVelDistr' label{iLabel,1} ','...
        '''angleWithNormal'',angleNormalDistr' label{iLabel,1} ','...
        '''angularVel'',angularVelDistr' label{iLabel,1} ');']);
    eval(['meanStdMin25P50P75PMax = struct(''distance'',sisterDistParam' label{iLabel,1} ','...
        '''posRateChangeDist'',sisterVelPosParam' label{iLabel,1} ','...
        '''negRateChangeDist'',sisterVelNegParam' label{iLabel,1} ','...
        '''angleWithNormal'',angleNormalParam' label{iLabel,1} ','...
        '''angularVel'',angularVelParam' label{iLabel,1} ');']);
    eval(['temporalTrend = struct(''distance'',sisterDistTrend' label{iLabel,1} ','...
        '''absoluteRateChangeDist'',sisterVelTrend' label{iLabel,1} ','...
        '''angleWithNormal'',angleNormalTrend' label{iLabel,1} ','...
        '''angularVel'',angularVelTrend' label{iLabel,1} ');']);
    eval(['autocorr = struct(''distance'',sisterDistAutocorr' label{iLabel,1} ','...
        '''rateChangeDist'',sisterVelAutocorr' label{iLabel,1} ','...
        '''angleWithNormal'',angleNormalAutocorr' label{iLabel,1} ','...
        '''angularVel'',angularVelAutocorr' label{iLabel,1} ');']);
    eval(['arma = struct(''distance'',sisterDistArma' label{iLabel,1} ','...
        '''rateChangeDist'',sisterVelArma' label{iLabel,1} ','...
        '''angleWithNormal'',angleNormalArma' label{iLabel,1} ','...
        '''angularVel'',angularVelArma' label{iLabel,1} ');']);
    eval(['numSistersCat = iGlobal' label{iLabel,1} ';']);

    eval([label{iLabel,1} ' = struct(''numSisters'',numSistersCat,'...
        '''distribution'',distribution,'...
        '''meanStdMin25P50P75PMax'',meanStdMin25P50P75PMax,'...
        '''temporalTrend'',temporalTrend,'...
        '''autocorr'',autocorr,''arma'',arma);'])

end

sisterConnection = struct('Inlier',Inlier,'Unaligned',Unaligned,...
    'Lagging',Lagging);

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
analysisStruct = makiMakeAnalysisPlatformIndependent(analysisStruct);
save(fullfile(dir2SaveRes,fileName),'analysisStruct');

%% plots

if verbose

    %get number of frames in each movie
    numFrames = dataStruct(1).dataProperties.movieSize(end);

    %% distance stuff %%
    
    for iLabel = goodLabel

        %open figure and write title
        figure('Name',[fileName(1:end-4) ' - Distances - ' label{iLabel,1}],'NumberTitle','off');

        %plot a sample of trajectories

        %create subplot 1
        subplot(2,1,1);
        hold on;

        %put all distances together in one matrix
        eval(['distanceMat = [sisterDist' label{iLabel,1} '.observations];'])
        distanceMat = distanceMat(:,1:2:end);

        %plot distance over time for all sisters
        plot((0:numFrames-1)*timeLapse,distanceMat);

        %set axes limit
        axis([0 (numFrames-1)*timeLapse 0 max(distanceMat(:))+1]);

        %write axes labels
        xlabel('Time (s)');
        ylabel('Sister separation (\mum)');

        %         %write averaging information
        %         eval(['text(timeLapse,max(distanceMat(:))+0.6,sprintf(['' Sister separation''' ...
        %             ''' (um): %4.2f +- %4.2f \n Rate of change of sister separation ''' ...
        %             '''(nm/s): \n +ve: %4.2f +- %4.2f, -ve: %4.2f +- %4.2f''],'...
        %             'sisterDistMeanStd' label{iLabel,1} '(1),sisterDistMeanStd' ...
        %             label{iLabel,1} '(2),sisterVelPosMeanStd' label{iLabel,1} '(1),'...
        %             'sisterVelPosMeanStd' label{iLabel,1} '(2),sisterVelNegMeanStd' ...
        %             label{iLabel,1} '(1),sisterVelNegMeanStd' label{iLabel,1} '(2)));']);

        %hold off figure
        hold off

        %plot autocorrelation functions

        %create subplot 2
        subplot(2,1,2);
        hold on;

        %plot the distance and velocity autocorrelations
        eval(['sisterDistAutocorr = sisterDistAutocorr' label{iLabel,1} ';']);
        plot((0:maxLag)*timeLapse,sisterDistAutocorr,'k','marker','.');
        eval(['sisterVelAutocorr = sisterVelAutocorr' label{iLabel,1} ';']);
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
        
    end %(for iLabel = goodLabel)

    %% angle stuff %%
    
    for iLabel = goodLabel

        %open figure and write title
        figure('Name',[fileName(1:end-4) ' - Angles - ' label{iLabel,1}],'NumberTitle','off');

        %plot a sample of time series of angle with normal

        %create subplot 1
        subplot(2,2,1);
        hold on;

        %put all angles with normal together in one matrix
        eval(['angleMat = [angleNormal' label{iLabel,1} '.observations];']);

        %plot angles with normal over time for all sisters
        plot((0:numFrames-1)*timeLapse,angleMat);

        %set axes limit
        axis([0 (numFrames-1)*timeLapse 0 90]);

        %write axes labels
        xlabel('Time (s)');
        ylabel('Angle with normal to metaphase plate (degrees)');

        %         %write averaging information
        %         eval(['text(timeLapse,80,sprintf(''angle (degrees): %4.2f +- %4.2f'','...
        %             'angleNormalMeanStd' label{iLabel,1} '(1),angleNormalMeanStd' ...
        %             label{iLabel,1} '(2)));']);

        %hold off figure
        hold off

        %plot autocorrelation function of angle with normal

        %create subplot 3
        subplot(2,2,3);
        hold on;

        %plot the autocorrelation of angle with normal
        eval(['angleNormalAutocorr = angleNormalAutocorr' label{iLabel,1} ';']);
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
        eval(['angleMat = [angularVel' label{iLabel,1} '.observations];'])

        %plot angular velocities over time for all sisters
        plot((0:numFrames-2)*timeLapse,angleMat);

        %set axes limit
        axis([0 (numFrames-2)*timeLapse 0 1.1*max(angleMat(:))]);

        %write axes labels
        xlabel('Time (s)');
        ylabel('Angular velocity (degrees/s)');

        %         %write averaging information
        %         eval(['text(timeLapse,1.05*max(angleMat(:)),sprintf(''angular velocity (degrees/s): %4.2f +- %4.2f'','...
        %             'angularVelMeanStd' label{iLabel,1} '(1),angularVelMeanStd' label{iLabel,1} '(2)));'])

        %hold off figure
        hold off

        %plot autocorrelation function of angle with normal

        %create subplot 4
        subplot(2,2,4);
        hold on;

        %plot the autocorrelation of angular velocity
        eval(['angularVelAutocorr = angularVelAutocorr' label{iLabel,1} ';']);
        plot((0:maxLag)*timeLapse,angularVelAutocorr,'k','marker','.');

        %set axes limit
        axis([0 maxLag*timeLapse min(0,1.1*min(angularVelAutocorr)) 1.1]);

        %write axes labels
        xlabel('Lag (s)');
        ylabel('Autocorrelation of angular velocity');

        %hold off figure
        hold off

    end %(for iLabel = goodLabel)
        
    %% angle vs. distance %%

    for iLabel = goodLabel

        %open figure and write title
        figure('Name',[fileName(1:end-4) ' - Angle vs. Distance - ' label{iLabel,1}],'NumberTitle','off');

        %plot angle vs. distance as a scatter plot
        eval(['plot(allDistances' label{iLabel,1} ',allAngles' label{iLabel,1} ',''k.'')']);

        %write axes labels
        xlabel('Sister separation (\mum)');
        ylabel('Angle with normal (degrees)');
        
    end %(or iLabel = goodLabel)
    
    %% temporal trend stuff %%
    
    for iLabel = goodLabel

        %open figure and write title
        figure('Name',[fileName(1:end-4) ' - Temporal trends - ' label{iLabel,1}],'NumberTitle','off');

        %create subplot 1
        subplot(4,2,1);
        hold on;

        %plot histogram of distance temporal trend
        eval(['trend2plot = sisterDistTrend' label{iLabel,1} ';']);
        [x] = histogram(trend2plot(:,1));
        hist(trend2plot(:,1),length(x));
        
        %write axes labels and title
        xlabel('Sister separation slopes (\mum/s)');
        ylabel('# of occurance');
        
        %hold off figure
        hold off

        %create subplot 2
        subplot(4,2,2);
        hold on;

        %plot histogram of significant trends only
        trend2plot = trend2plot(trend2plot(:,3)<0.05,1);
        [x] = histogram(trend2plot);
        hist(trend2plot,length(x));
        
        %write axes labels and title
        xlabel('Significant slopes only');
        ylabel('# of occurance');
        
        %hold off figure
        hold off

        %create subplot 3
        subplot(4,2,3);
        hold on;

        %plot histogram of velocity temporal trend
        eval(['trend2plot = sisterVelTrend' label{iLabel,1} ';']);
        [x] = histogram(trend2plot(:,1));
        hist(trend2plot(:,1),length(x));
        
        %write axes labels and title
        xlabel('Absolute rate change sister separation slopes (nm/s/s)');
        ylabel('# of occurances');
        
        %hold off figure
        hold off

        %create subplot 4
        subplot(4,2,4);
        hold on;

        %plot histogram of significant trends only
        trend2plot = trend2plot(trend2plot(:,3)<0.05,1);
        [x] = histogram(trend2plot);
        hist(trend2plot,length(x));
        
        %write axes labels and title
        xlabel('Significant slopes only');
        ylabel('# of occurances');
        
        %hold off figure
        hold off

        %create subplot 5
        subplot(4,2,5);
        hold on;

        %plot histogram of angle with normal temporal trend
        eval(['trend2plot = angleNormalTrend' label{iLabel,1} ';']);
        [x] = histogram(trend2plot(:,1));
        hist(trend2plot(:,1),length(x));
        
        %write axes labels and title
        xlabel('Angle with normal slopes (degrees/s)');
        ylabel('# of occurances');
        
        %hold off figure
        hold off

        %create subplot 6
        subplot(4,2,6);
        hold on;

        %plot histogram of significant trends only
        trend2plot = trend2plot(trend2plot(:,3)<0.05,1);
        [x] = histogram(trend2plot);
        hist(trend2plot,length(x));
        
        %write axes labels and title
        xlabel('Significant slopes only');
        ylabel('# of occurances');
        
        %hold off figure
        hold off

        %create subplot 7
        subplot(4,2,7);
        hold on;

        %plot histogram of angular velocity temporal trend
        eval(['trend2plot = angularVelTrend' label{iLabel,1} ';']);
        [x] = histogram(trend2plot(:,1));
        hist(trend2plot(:,1),length(x));
        
        %write axes labels and title
        xlabel('Angular velocity slopes  (degrees/s/s)');
        ylabel('# of occurance');
        
        %hold off figure
        hold off

        %create subplot 8
        subplot(4,2,8);
        hold on;

        %plot histogram of angular velocity temporal trend
        trend2plot = trend2plot(trend2plot(:,3)<0.05,1);
        [x] = histogram(trend2plot);
        hist(trend2plot,length(x));
        
        %write axes labels and title
        xlabel('Significant slopes only');
        ylabel('# of occurance');
        
        %hold off figure
        hold off

    end %(or iLabel = goodLabel)    
    
end

%% ~~~ the end ~~~ %%

