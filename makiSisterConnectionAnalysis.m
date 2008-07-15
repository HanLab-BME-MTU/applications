function analysisStruct = makiSisterConnectionAnalysis(jobType,...
    analysisStruct,verbose,removeNetDisp,randomize,samplingPeriod)
%MAKISISTERCONNECTIONANALYSIS analyzes in multiple ways the connection between sisters
%
%SYNOPSIS analysisStruct = makiSisterConnectionAnalysis(jobType,...
%    analysisStruct,verbose,removeNetDisp,randomize,samplingPeriod)
%
%INPUT  jobType       : string which can take the values:
%                       'TEST', 'HERCULES', 'DANUSER', 'MERALDI',
%                       'SWEDLOW' or 'MCAINSH'
%       analysisStruct: Structure with field movies indicating the movies
%                       to be analyzed. Optional. If not input, GUI to load
%                       movies is launched.
%       verbose       : 1 to make plots, 0 otherwise. Optional. Default: 0.
%       removeNetDisp : 1 - Make center of each sister pair at zero
%                       throughout the whole movie, 0 otherwise.
%                       Optional. Default: 0.
%       randomize     : 1 - Randomize sister pairing, 0 otherwise.
%                       Optional. Default: 0.
%       samplingPeriod: 1 to keep sampling as is, 2 to downsample by taking
%                       every 2nd time point, 3 to downsample by taking
%                       every 3rd time point, etc.
%                       Optional. Default: 1.
%
%OUTPUT analysisStruct: Same as input but with additional field
%           .sisterConnection: 
%
%REMARKS Code is not applicable to anaphase movies/frames
%
%Khuloud Jaqaman, July 2007

%% input
if nargin < 1 || isempty(jobType)
    jobType = 'DANUSER';
end

if nargin < 3 || isempty(verbose)
    verbose = 0;
end

if nargin < 4 || isempty(removeNetDisp)
    removeNetDisp = 0;
end

if nargin < 5 || isempty(randomize)
    randomize = 0;
end

if nargin < 6 || isempty(samplingPeriod)
    samplingPeriod = 1;
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

%multiply number of sisters by samplingPeriod to get the effective number
%of sisters
numSistersTotEff = numSistersTot * samplingPeriod;

%get time between frames
timeLapse = round(2*dataStruct(1).dataProperties.timeLapse)/2 * samplingPeriod;

%% collect sister angles and distances

%define cell array of sister labels
label{1,1} = 'Inlier';
label{2,1} = 'Unaligned';
label{3,1} = 'Lagging';

%initialize global sister index and structures
for iLabel = 1 : 3

    eval(['iGlobal' label{iLabel,1} ' = 0;'])

    eval(['sisterDist' label{iLabel,1} '(1:numSistersTotEff,1) = struct(''observations'',[]);'])
    eval(['sisterVel' label{iLabel,1} '(1:numSistersTotEff,1) = struct(''observations'',[]);'])
    eval(['angleNormal' label{iLabel,1} '(1:numSistersTotEff,1) = struct(''observations'',[]);'])
    eval(['angularVel' label{iLabel,1} '(1:numSistersTotEff,1) = struct(''observations'',[]);'])
    
    eval(['movieStartIndx' label{iLabel,1} ' = zeros(numMovies,1);'])
    eval(['movieEndIndx' label{iLabel,1} ' = zeros(numMovies,1);'])
    
end

%go over all movies
for iMovie = 1 : numMovies
    
    %store the index of the first place where the sisters belonging to this
    %movie are stored
    movieStartIndxInlier(iMovie) = iGlobalInlier + 1;
    movieStartIndxUnaligned(iMovie) = iGlobalUnaligned + 1;
    movieStartIndxLagging(iMovie) = iGlobalLagging + 1;
    
    %if there are sister kinetochores in this movie
    if numSisters(iMovie) > 0

        %construct sisters with aligned coordinates
        %remove net displacement and randomize if requested
        sisterListTmp = makiConstructAlignedSisters(dataStruct(iMovie),...
            removeNetDisp,randomize);
        numSisters(iMovie) = length(sisterListTmp);

        %sister type = 0 if all kinetochores are inliers
        %sister type = 1 if some kinetochores are unaligned
        %sister type = 2 if some kinetochores are lagging
        sisterType = zeros(numSisters(iMovie),1);
        sisterType(dataStruct(iMovie).updatedClass(1).sistersUnaligned) = 1;
        sisterType(dataStruct(iMovie).updatedClass(1).sistersLagging) = 2;

        %subsample as requested
        for iSample = 1 : samplingPeriod

            %get coordinates with proper subsampling
            for iSister = 1 : numSisters(iMovie)
                sisterList(iSister).distanceAligned = sisterListTmp(iSister).distanceAligned(iSample:samplingPeriod:end,:);
                sisterList(iSister).coords1Aligned = sisterListTmp(iSister).coords1Aligned(iSample:samplingPeriod:end,:);
                sisterList(iSister).coords2Aligned = sisterListTmp(iSister).coords2Aligned(iSample:samplingPeriod:end,:);
            end
            
            %copy fields out of dataStruct(iMovie)
            planeFit = dataStruct(iMovie).planeFit(iSample:samplingPeriod:end);
            updatedClass = dataStruct(iMovie).updatedClass(iSample:samplingPeriod:end);
            numFramesMovie = length(planeFit);

            %determine frames where there is a plane
            framesWithPlane = [];
            for t = 1 : numFramesMovie
                if ~isempty(planeFit(t).planeVectors)
                    framesWithPlane = [framesWithPlane t];
                end
            end
            
            %find frame where anaphase starts (if it starts at all)
            framePhase = vertcat(updatedClass.phase);
            firstFrameAna = find(framePhase=='a',1,'first');
            if isempty(firstFrameAna)
                firstFrameAna = numFramesMovie + 1;
            end

            %go over all sisters in movie
            for iSister = 1 : numSisters(iMovie)

                %find number of frames and frames where pair "exists" before
                %anaphase
                goodFrames = ~isnan(sisterList(iSister).distanceAligned(:,1));
                numFrames = length(goodFrames);
                goodFrames = find(goodFrames);
                goodFrames = goodFrames(goodFrames < firstFrameAna);

                %get sister coordinates in good frames
                coords1 = NaN(numFrames,3);
                coords2 = NaN(numFrames,3);
                coords1Std = NaN(numFrames,3);
                coords2Std = NaN(numFrames,3);
                for iFrame = goodFrames'
                    coords1(iFrame,:) = sisterList(iSister).coords1Aligned(iFrame,1:3);
                    coords2(iFrame,:) = sisterList(iSister).coords2Aligned(iFrame,1:3);
                    coords1Std(iFrame,:) = sisterList(iSister).coords1Aligned(iFrame,4:6);
                    coords2Std(iFrame,:) = sisterList(iSister).coords2Aligned(iFrame,4:6);
                end

                %calculate vector between sisters
                sisterVec = coords2 - coords1; %um
                sisterVecVar = coords1Std.^2 + coords2Std.^2; %um^2

                %calculate distance between sisters
                sisterDistance = sqrt(sum(sisterVec.^2,2)); %um
                sisterDistStd = sqrt( sum( sisterVecVar .* sisterVec.^2 ,2) ) ...
                    ./ sisterDistance; %um
                
                %calculate sister velocity
                sisterVelocity = diff(sisterDistance) * 1000 / timeLapse; %nm/s
                sisterVelStd = sqrt( sum( [sisterDistStd(2:end) sisterDistStd(1:end-1)].^2 ,2) ) * 1000 / timeLapse; %nm/s

                %calculate angle with normal
                angleWithNorm = NaN(numFrames,1);
                for iFrame = framesWithPlane
                    angleWithNorm(iFrame) = acos(abs(sisterVec(iFrame,1)./norm(sisterVec(iFrame,1:3)))); %radians
                end

                %calculate angle between consecutive frames
                angularDisp = NaN(numFrames-1,1);
                for iFrame = 1 : numFrames - 1
                    angularDisp(iFrame) = acos(abs(sisterVec(iFrame,:) * sisterVec(iFrame+1,:)' ...
                        / norm(sisterVec(iFrame,:)) / norm(sisterVec(iFrame+1,:)))); %radians
                end
                
                %store sister information based on the sister type
                iLabel = sisterType(iSister) + 1;

                %increase global index of sister type by 1
                eval(['iGlobal' label{iLabel,1} ' = iGlobal' label{iLabel,1} ' + 1;'])

                %store information
                eval(['sisterDist' label{iLabel,1} '(iGlobal' label{iLabel,1} ...
                    ').observations = [sisterDistance sisterDistStd];']) %um
                eval(['sisterVel' label{iLabel,1} '(iGlobal' label{iLabel,1} ...
                    ').observations = [sisterVelocity sisterVelStd];']) %nm/s
                eval(['angleNormal' label{iLabel,1} '(iGlobal' label{iLabel,1} ...
                    ').observations = angleWithNorm * 180 / pi;']) %deg
                eval(['angularVel' label{iLabel,1} '(iGlobal' label{iLabel,1} ...
                    ').observations = angularDisp * 180 / pi / timeLapse;']) %deg/s

            end %(for iSister = 1 : numSisters(iMovie))

        end %(for iSample = 1 : samplingPeriod)

    end %(if numSisters(iMovie) > 0)

    %store the index of the last place where the sisters belonging to this
    %movie are stored
    movieEndIndxInlier(iMovie) = iGlobalInlier;
    movieEndIndxUnaligned(iMovie) = iGlobalUnaligned;
    movieEndIndxLagging(iMovie) = iGlobalLagging;
    
end %(for iMovie = 1 : numMovies)
    
%store number of sisters per category
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
    
    eval(['sisterDistIndParam' label{iLabel,1} ' = NaN(numMovies,7);'])
    eval(['sisterVelPosIndParam' label{iLabel,1} ' = NaN(numMovies,7);'])
    eval(['sisterVelNegIndParam' label{iLabel,1} ' = NaN(numMovies,7);'])
    eval(['angleNormalIndParam' label{iLabel,1} ' = NaN(numMovies,7);'])
    eval(['angularVelIndParam' label{iLabel,1} ' = NaN(numMovies,7);'])    
end

%calculation
for iLabel = goodLabel

    % overall %
    
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
    
    % individual cells %
    
    %distance
    for iMovie = 1 : numMovies
        eval(['allValues = vertcat(sisterDist' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie)).observations);']);
        if ~isempty(allValues)
            allValues = allValues(:,1);
            allValues = allValues(~isnan(allValues));
            if ~isempty(allValues)
                eval(['sisterDistIndParam' label{iLabel,1} '(iMovie,:)' ...
                    ' = [mean(allValues) std(allValues) min(allValues) ' ...
                    'prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) '...
                    'max(allValues)];']);
            end
        end
    end

    %positive velocity
    for iMovie = 1 : numMovies
        eval(['allValues = vertcat(sisterVel' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie)).observations);']);
        if ~isempty(allValues)
            allValues = allValues(allValues(:,1)>0,1);
            allValues = allValues(~isnan(allValues));
            if ~isempty(allValues)
                eval(['sisterVelPosIndParam' label{iLabel,1} '(iMovie,:)' ...
                    ' = [mean(allValues) std(allValues) min(allValues) ' ...
                    'prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) '...
                    'max(allValues)];']);
            end
        end
    end

    %negative velocity
    for iMovie = 1 : numMovies
        eval(['allValues = vertcat(sisterVel' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie)).observations);']);
        if ~isempty(allValues)
            allValues = abs(allValues(allValues(:,1)<0,1));
            allValues = allValues(~isnan(allValues));
            if ~isempty(allValues)
                eval(['sisterVelNegIndParam' label{iLabel,1} '(iMovie,:)' ...
                    ' = [mean(allValues) std(allValues) min(allValues) ' ...
                    'prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) '...
                    'max(allValues)];']);
            end
        end
    end

    %angle with normal
    for iMovie = 1 : numMovies
        eval(['allValues = vertcat(angleNormal' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie)).observations);']);
        if ~isempty(allValues)
            allValues = allValues(:,1);
            allValues = allValues(~isnan(allValues));
            if ~isempty(allValues)
                eval(['angleNormalIndParam' label{iLabel,1} '(iMovie,:)' ...
                    ' = [mean(allValues) std(allValues) min(allValues) ' ...
                    'prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) '...
                    'max(allValues)];']);
            end
        end
    end

    %angular velocity
    for iMovie = 1 : numMovies
        eval(['allValues = vertcat(angularVel' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie)).observations);']);
        if ~isempty(allValues)
            allValues = allValues(:,1);
            allValues = allValues(~isnan(allValues));
            if ~isempty(allValues)
                eval(['angularVelIndParam' label{iLabel,1} '(iMovie,:)' ...
                    ' = [mean(allValues) std(allValues) min(allValues) ' ...
                    'prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) '...
                    'max(allValues)];']);
            end
        end
    end
    
end

%% temporal trends

%initialization
for iLabel = 1 : 3
    eval(['sisterDistTrend' label{iLabel,1} ' = [];'])
    eval(['sisterVelTrend' label{iLabel,1} ' = [];'])
    eval(['angleNormalTrend' label{iLabel,1} ' = [];'])
    eval(['angularVelTrend' label{iLabel,1} ' = [];'])
end

%% autocorrelation

%define maximum lag
maxLag = round(20 / samplingPeriod);
maxLagSis = round(10 / samplingPeriod);

%initialization
for iLabel = 1 : 3
    eval(['sisterDistAutocorr' label{iLabel,1} ' = [];'])
    eval(['sisterVelAutocorr' label{iLabel,1} ' = [];'])
    eval(['angleNormalAutocorr' label{iLabel,1} ' = [];'])
    eval(['angularVelAutocorr' label{iLabel,1} ' = [];'])
    
    eval(['sisterDistIndAutocorr' label{iLabel,1} ' = NaN(maxLag+1,2,numMovies);'])
    eval(['sisterVelIndAutocorr' label{iLabel,1} ' = NaN(maxLag+1,2,numMovies);'])
    eval(['angleNormalIndAutocorr' label{iLabel,1} ' = NaN(maxLag+1,2,numMovies);'])
    eval(['angularVelIndAutocorr' label{iLabel,1} ' = NaN(maxLag+1,2,numMovies);'])
    
    eval(['sisterDistSisAutocorr' label{iLabel,1} ' = NaN(maxLagSis+1,2,iGlobal' label{iLabel,1} ');'])
    eval(['sisterVelSisAutocorr' label{iLabel,1} ' = NaN(maxLagSis+1,2,iGlobal' label{iLabel,1} ');'])
    eval(['angleNormalSisAutocorr' label{iLabel,1} ' = NaN(maxLagSis+1,2,iGlobal' label{iLabel,1} ');'])
    eval(['angularVelSisAutocorr' label{iLabel,1} ' = NaN(maxLagSis+1,2,iGlobal' label{iLabel,1} ');'])
end

%calculation
for iLabel = goodLabel

    % overall %
    
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
    
    % individual cells %
    
    %distance
    for iMovie = 1 : numMovies
        eval(['traj = sisterDist' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie));'])
        if ~isempty(traj)
            [tmpCorr,errFlag] = autoCorr(traj,maxLag);
            if ~errFlag
                eval(['sisterDistIndAutocorr' label{iLabel,1} '(:,:,iMovie) = tmpCorr;'])
            end
        end
    end

    %velocity
    for iMovie = 1 : numMovies
        eval(['traj = sisterVel' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie));'])
        if ~isempty(traj)
            [tmpCorr,errFlag] = autoCorr(traj,maxLag);
            if ~errFlag
                eval(['sisterVelIndAutocorr' label{iLabel,1} '(:,:,iMovie) = tmpCorr;'])
            end
        end
    end

    %angle with normal
    for iMovie = 1 : numMovies
        eval(['traj = angleNormal' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie));'])
        if ~isempty(traj)
            [tmpCorr,errFlag] = autoCorr(traj,maxLag);
            if ~errFlag
                eval(['angleNormalIndAutocorr' label{iLabel,1} '(:,:,iMovie) = tmpCorr;'])
            end
        end
    end

    %angular velocity
    for iMovie = 1 : numMovies
        eval(['traj = angularVel' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie));'])
        if ~isempty(traj)
            [tmpCorr,errFlag] = autoCorr(traj,maxLag);
            if ~errFlag
                eval(['angularVelIndAutocorr' label{iLabel,1} '(:,:,iMovie) = tmpCorr;'])
            end
        end
    end    
    
    % individual sisters %
    
    %distance
    for iSister = 1 : eval(['iGlobal' label{iLabel,1}])
        eval(['traj = sisterDist' label{iLabel,1} '(iSister).observations;'])
        if length(find(~isnan(traj(:,1))))>maxLagSis+10
            [tmpCorr,errFlag] = autoCorr(traj,maxLagSis);
            if ~errFlag
                eval(['sisterDistSisAutocorr' label{iLabel,1} '(:,:,iSister) = tmpCorr;'])
            end
        end
    end
    
    %velocity
    for iSister = 1 : eval(['iGlobal' label{iLabel,1}])
        eval(['traj = sisterVel' label{iLabel,1} '(iSister).observations;'])
        if length(find(~isnan(traj(:,1))))>maxLagSis+10
            [tmpCorr,errFlag] = autoCorr(traj,maxLagSis);
            if ~errFlag
                eval(['sisterVelSisAutocorr' label{iLabel,1} '(:,:,iSister) = tmpCorr;'])
            end
        end
    end
    
    %angle with normal
    for iSister = 1 : eval(['iGlobal' label{iLabel,1}])
        eval(['traj = angleNormal' label{iLabel,1} '(iSister).observations;'])
        if length(find(~isnan(traj(:,1))))>maxLagSis+10
            [tmpCorr,errFlag] = autoCorr(traj,maxLagSis);
            if ~errFlag
                eval(['angleNormalSisAutocorr' label{iLabel,1} '(:,:,iSister) = tmpCorr;'])
            end
        end
    end
    
    %angular velocity
    for iSister = 1 : eval(['iGlobal' label{iLabel,1}])
        eval(['traj = angularVel' label{iLabel,1} '(iSister).observations;'])
        if length(find(~isnan(traj(:,1))))>maxLagSis+10
            [tmpCorr,errFlag] = autoCorr(traj,maxLagSis);
            if ~errFlag
                eval(['angularVelSisAutocorr' label{iLabel,1} '(:,:,iSister) = tmpCorr;'])
            end
        end
    end
    
end

%% ARMA

%initialization
for iLabel = 1 : 3
    eval(['sisterDistIndArma' label{iLabel,1} ' = repmat(struct(''results'',[]),numMovies,1);'])
    eval(['sisterVelIndArma' label{iLabel,1} ' = repmat(struct(''results'',[]),numMovies,1);'])
    eval(['angleNormalIndArma' label{iLabel,1} ' = repmat(struct(''results'',[]),numMovies,1);'])
    eval(['angularVelIndArma' label{iLabel,1} ' = repmat(struct(''results'',[]),numMovies,1);'])
end

%define model orders to test
modelOrder = [0 5; 0 5; -1 -1];

%calculation
for iLabel = goodLabel
    
    %     %call ARMA analysis function for sister distance
    %     for iMovie = 1 : numMovies
    %         eval(['traj = sisterDist' label{iLabel,1} '(movieStartIndx' ...
    %             label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie));'])
    %         for i=1:length(traj)
    %             traj(i).observations(:,2) = 0;
    %         end
    %         if ~isempty(traj)
    %             fitResults = armaxFitKalmanMEX(traj,[],modelOrder,'tl');
    %             eval(['sisterDistIndArma' label{iLabel,1} '(iMovie).results = fitResults;'])
    %         end
    %     end
    
    %     %call ARMA analysis function for sister velocity
    %     for iMovie = 1 : numMovies
    %         eval(['traj = sisterVel' label{iLabel,1} '(movieStartIndx' ...
    %             label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie));'])
    %         for i=1:length(traj)
    %             traj(i).observations(:,2) = 0;
    %         end
    %         if ~isempty(traj)
    %             fitResults = armaxFitKalmanMEX(traj,[],modelOrder,'tl');
    %             eval(['sisterVelIndArma' label{iLabel,1} '(iMovie).results = fitResults;'])
    %         end
    %     end

    %     %call ARMA analysis function for angle with normal to plane
    %     for iMovie = 1 : numMovies
    %         eval(['fitResults = armaxFitKalmanMEX(angleNormal' label{iLabel,1} ...
    %             '(movieStartIndx' label{iLabel,1} '(iMovie):movieEndIndx' ...
    %             label{iLabel,1} '(iMovie)),[],modelOrder,''tl'');'])
    %         eval(['angleNormalIndArma' label{iLabel,1} '(iMovie).results = fitResults;'])
    %     end

    %     %call ARMA analysis function for angular velocity
    %     for iMovie = 1 : numMovies
    %         eval(['fitResults = armaxFitKalmanMEX(angularVel' label{iLabel,1} ...
    %             '(movieStartIndx' label{iLabel,1} '(iMovie):movieEndIndx' ...
    %             label{iLabel,1} '(iMovie)),[],modelOrder,''tl'');'])
    %         eval(['angularVelIndArma' label{iLabel,1} '(iMovie).results = fitResults;'])
    %     end

end

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
    eval(['meanStdMin25P50P75PMax.all = struct(''distance'',sisterDistParam' label{iLabel,1} ','...
        '''posRateChangeDist'',sisterVelPosParam' label{iLabel,1} ','...
        '''negRateChangeDist'',sisterVelNegParam' label{iLabel,1} ','...
        '''angleWithNormal'',angleNormalParam' label{iLabel,1} ','...
        '''angularVel'',angularVelParam' label{iLabel,1} ');']);
    eval(['meanStdMin25P50P75PMax.indcell = struct(''distance'',sisterDistIndParam' label{iLabel,1} ','...
        '''posRateChangeDist'',sisterVelPosIndParam' label{iLabel,1} ','...
        '''negRateChangeDist'',sisterVelNegIndParam' label{iLabel,1} ','...
        '''angleWithNormal'',angleNormalIndParam' label{iLabel,1} ','...
        '''angularVel'',angularVelIndParam' label{iLabel,1} ');']);
    eval(['temporalTrend = struct(''distance'',sisterDistTrend' label{iLabel,1} ','...
        '''absoluteRateChangeDist'',sisterVelTrend' label{iLabel,1} ','...
        '''angleWithNormal'',angleNormalTrend' label{iLabel,1} ','...
        '''angularVel'',angularVelTrend' label{iLabel,1} ');']);
    eval(['autocorr.all = struct(''distance'',sisterDistAutocorr' label{iLabel,1} ','...
        '''rateChangeDist'',sisterVelAutocorr' label{iLabel,1} ','...
        '''angleWithNormal'',angleNormalAutocorr' label{iLabel,1} ','...
        '''angularVel'',angularVelAutocorr' label{iLabel,1} ');']);
    eval(['autocorr.indcell = struct(''distance'',sisterDistIndAutocorr' label{iLabel,1} ','...
        '''rateChangeDist'',sisterVelIndAutocorr' label{iLabel,1} ','...
        '''angleWithNormal'',angleNormalIndAutocorr' label{iLabel,1} ','...
        '''angularVel'',angularVelIndAutocorr' label{iLabel,1} ');']);
    eval(['autocorr.indsis = struct(''distance'',sisterDistSisAutocorr' label{iLabel,1} ','...
        '''rateChangeDist'',sisterVelSisAutocorr' label{iLabel,1} ','...
        '''angleWithNormal'',angleNormalSisAutocorr' label{iLabel,1} ','...
        '''angularVel'',angularVelSisAutocorr' label{iLabel,1} ');']);
    eval(['arma.ind = struct(''distance'',sisterDistIndArma' label{iLabel,1} ','...
        '''rateChangeDist'',sisterVelIndArma' label{iLabel,1} ','...
        '''angleWithNormal'',angleNormalIndArma' label{iLabel,1} ','...
        '''angularVel'',angularVelIndArma' label{iLabel,1} ');']);
    eval(['numSistersCat = iGlobal' label{iLabel,1} ';']);

    eval([label{iLabel,1} ' = struct('...
        '''numSisters'',numSistersCat,'...
        '''distribution'',distribution,'...
        '''meanStdMin25P50P75PMax'',meanStdMin25P50P75PMax,'...
        '''temporalTrend'',temporalTrend,'...
        '''autocorr'',autocorr,'...
        '''arma'',arma);'])

end

inputParam = struct('removeNetDisp',removeNetDisp,'randomize',randomize,'samplingPeriod',samplingPeriod);
sisterConnection = struct('Inlier',Inlier,'Unaligned',Unaligned,...
    'Lagging',Lagging,'inputParam',inputParam);

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
analysisStruct = makiMakeAnalysisPlatformIndependent(analysisStruct,jobType);
save(fullfile(dir2SaveRes,fileName),'analysisStruct');

%% plots

if verbose

    %get number of frames in each movie
    numFrames = NaN(numMovies,1);
    for iMovie = 1 : numMovies
        numFrames(iMovie) = floor(dataStruct(iMovie).dataProperties.movieSize(end)/samplingPeriod);
    end
    numFrames = min(numFrames);

    %% distance stuff %%
    
    for iLabel = goodLabel

        %open figure and write title
        figFileName = [fileName(1:end-4) '-Distances-' label{iLabel,1}];
        figHandle = figure('Name',figFileName,'NumberTitle','off');

        %plot a sample of trajectories

        %create subplot 1
        subplot(2,2,1);
        hold on

        %put all distances together in one matrix
        distanceMat = [];
        for iSis = 1 : eval(['length(sisterDist' label{iLabel,1} ')'])
            eval(['distanceMat = [distanceMat sisterDist' label{iLabel,1} '(iSis).observations(1:numFrames,1)];'])
        end

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

        %hold off subplot 1
        hold off

        %plot overall autocorrelation functions

        %create subplot 2
        subplot(2,2,2);
        hold on

        %plot the distance and velocity autocorrelations
        eval(['sisterDistAutocorr = sisterDistAutocorr' label{iLabel,1} ';']);
        plot((0:maxLag)*timeLapse,sisterDistAutocorr(:,1),'k','marker','.');
        eval(['sisterVelAutocorr = sisterVelAutocorr' label{iLabel,1} ';']);
        plot((0:maxLag)*timeLapse,sisterVelAutocorr(:,1),'r','marker','.');
        plot([0 maxLag]*timeLapse,[0 0],'k--');

        %set axes limit
        axis([0 maxLag*timeLapse min(0,1.1*min([sisterDistAutocorr(:,1);sisterVelAutocorr(:,1)])) 1.1]);

        %write axes labels
        xlabel('Lag (s)');
        ylabel('Autocorrelation - overall');

        %write legend
        text(1*timeLapse,0.9,sprintf([' Black: Sister separation \n Red: ' ...
            'Rate of change of sister separation']));

        %hold off subplot 2
        hold off
        
        %plot individual autocorrelation functions
        
        %create subplot 3
        subplot(2,2,3)
        hold on
        
        %plot the distance autocorrelations
        eval(['tmpCorr = squeeze(sisterDistIndAutocorr' label{iLabel,1} '(:,1,:));'])
        plot((0:maxLag)*timeLapse,tmpCorr,'marker','.');
        plot([0 maxLag]*timeLapse,[0 0],'k--');
                
        %set axes limit
        axis([0 maxLag*timeLapse min(0,1.1*min(tmpCorr(:))) 1.1]);
        
        %write axes labels
        xlabel('Lag (s)');
        ylabel('Sister separation autocorrelation - ind movies');

        %hold off subplot 3
        hold off
        
        %create subplot 4
        subplot(2,2,4)
        hold on
        
        %plot the velocity autocorrelations
        eval(['tmpCorr = squeeze(sisterVelIndAutocorr' label{iLabel,1} '(:,1,:));'])
        plot((0:maxLag)*timeLapse,tmpCorr,'marker','.');
        plot([0 maxLag]*timeLapse,[0 0],'k--');

        %set axes limit
        axis([0 maxLag*timeLapse min(0,1.1*min(tmpCorr(:))) 1.1]);
        
        %write axes labels
        xlabel('Lag (s)');
        ylabel('Rate change sister separation autocorrelation - ind movies');

        %hold off subplot 4
        hold off
        
        %save figure in file
        saveas(figHandle,fullfile(dir2SaveRes,figFileName),'fig');
        
    end %(for iLabel = goodLabel)

    %% angle stuff %%
    
    for iLabel = goodLabel

        %open figure and write title
        figFileName = [fileName(1:end-4) '-Angles-' label{iLabel,1}];
        figHandle = figure('Name',figFileName,'NumberTitle','off');

        %plot a sample of time series of angle with normal

        %create subplot 1
        subplot(2,2,1);
        hold on

        %put all angles with normal together in one matrix
        angleMat = [];
        for iSis = 1 : eval(['length(angleNormal' label{iLabel,1} ')'])
            eval(['angleMat = [angleMat angleNormal' label{iLabel,1} '(iSis).observations(1:numFrames,1)];'])
        end

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

        %hold off subplot 1
        hold off

        %plot autocorrelation function of angle with normal

        %create subplot 3
        subplot(2,2,3);
        hold on

        %plot the autocorrelation of angle with normal
        eval(['angleNormalAutocorr = angleNormalAutocorr' label{iLabel,1} ';']);
        plot((0:maxLag)*timeLapse,angleNormalAutocorr(:,1),'k','marker','.');
        plot([0 maxLag]*timeLapse,[0 0],'k--');

        %set axes limit
        axis([0 maxLag*timeLapse min(0,1.1*min(angleNormalAutocorr(:,1))) 1.1]);

        %write axes labels
        xlabel('Lag (s)');
        ylabel('Autocorrelation of angle with normal');

        %hold off subplot 3
        hold off

        %plot a sample of time series of angular velocity

        %create subplot 2
        subplot(2,2,2);
        hold on

        %put all angular velocities together in one matrix
        angleMat = [];
        for iSis = 1 : eval(['length(angularVel' label{iLabel,1} ')'])
            eval(['angleMat = [angleMat angularVel' label{iLabel,1} '(iSis).observations(1:numFrames-1,1)];'])
        end

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

        %hold off subplot 2
        hold off

        %plot autocorrelation function of angle with normal

        %create subplot 4
        subplot(2,2,4);
        hold on

        %plot the autocorrelation of angular velocity
        eval(['angularVelAutocorr = angularVelAutocorr' label{iLabel,1} ';']);
        plot((0:maxLag)*timeLapse,angularVelAutocorr(:,1),'k','marker','.');
        plot([0 maxLag]*timeLapse,[0 0],'k--');

        %set axes limit
        axis([0 maxLag*timeLapse min(0,1.1*min(angularVelAutocorr(:,1))) 1.1]);

        %write axes labels
        xlabel('Lag (s)');
        ylabel('Autocorrelation of angular velocity');

        %hold off subplot 4
        hold off

        %save figure in file
        saveas(figHandle,fullfile(dir2SaveRes,figFileName),'fig');
        
    end %(for iLabel = goodLabel)
        
    %% angle vs. distance %%

    for iLabel = goodLabel

        %open figure and write title
        figFileName = [fileName(1:end-4) '-AngleVsDistance-' label{iLabel,1}];
        figHandle = figure('Name',figFileName,'NumberTitle','off');

        %plot angle vs. distance as a scatter plot
        eval(['plot(allDistances' label{iLabel,1} ',allAngles' label{iLabel,1} ',''k.'')']);

        %write axes labels
        xlabel('Sister separation (\mum)');
        ylabel('Angle with normal (degrees)');
        
        %save figure in file
        saveas(figHandle,fullfile(dir2SaveRes,figFileName),'fig');
        
    end %(or iLabel = goodLabel)
    
    %% histograms %%
    
    for iLabel = goodLabel
        
        %open figure and write title
        figFileName = [fileName(1:end-4) '-Histograms-' label{iLabel,1}];
        figHandle = figure('Name',figFileName,'NumberTitle','off');

        %create subplot 1 for distances
        subplot(2,2,1);
        hold on
        
        %get number of necessary bins
        eval(['[n,x] = histogram(sisterDistDistr' label{iLabel,1} ');']);
        eval(['hist(sisterDistDistr' label{iLabel,1} ',length(n));']);

        %write axes labels
        xlabel('Sister separation (\mum)');
        ylabel('# of occurances');        
    
        hold off
    
        %create subplot 2 for velocities
        subplot(2,2,2);
        hold on
        
        %get number of necessary bins
        eval(['[n,x] = histogram(sisterVelDistr' label{iLabel,1} ');']);
        eval(['hist(sisterVelDistr' label{iLabel,1} ',length(n));']);

        %write axes labels
        xlabel('Rate change sister separation (nm/s)');
        ylabel('# of occurances');        
    
        hold off
    
        %create subplot 3 for angles
        subplot(2,2,3);
        hold on
        
        %get number of necessary bins
        if eval(['~isempty(angleNormalDistr' label{iLabel,1} ')'])
            eval(['[n,x] = histogram(angleNormalDistr' label{iLabel,1} ');']);
            eval(['hist(angleNormalDistr' label{iLabel,1} ',length(n));']);
        end

        %write axes labels
        xlabel('Sister angle with normal (deg)');
        ylabel('# of occurances');        
    
        hold off
    
        %create subplot 4 for angular velocity
        subplot(2,2,4);
        hold on
        
        %get number of necessary bins
        eval(['[n,x] = histogram(angularVelDistr' label{iLabel,1} ');']);
        eval(['hist(angularVelDistr' label{iLabel,1} ',length(n));']);

        %write axes labels
        xlabel('Sister angular velocity (deg/s)');
        ylabel('# of occurances');        
    
        hold off
    
        %save figure in file
        saveas(figHandle,fullfile(dir2SaveRes,figFileName),'fig');
        
    end
    
end

%% ~~~ the end ~~~ %%

%% unused temporal trend stuff

% %calculation
% for iLabel = goodLabel
% 
%     %get number of sisters in this category
%     eval(['numSistersCat = iGlobal' label{iLabel,1} ';']);
%     
%     %distance
%     trendTmp = NaN(numSistersCat,3);
%     for iSister = 1 : numSistersCat
%         eval(['allValues = sisterDist' label{iLabel,1} '(iSister).observations(:,1);'])
%         indxAvail = find(~isnan(allValues));
%         [lineParam,S] = polyfit(indxAvail*timeLapse,allValues(indxAvail),1);
%         varCovMat = (inv(S.R)*inv(S.R)')*S.normr^2/S.df;
%         slopeStd = sqrt(varCovMat(1));
%         testStat = lineParam(1)/slopeStd;
%         if testStat > 0
%             pValue = 1 - tcdf(testStat,S.df);
%         else
%             pValue = tcdf(testStat,S.df);            
%         end
%         trendTmp(iSister,:) = [lineParam(1) slopeStd pValue];
%     end
%     eval(['sisterDistTrend' label{iLabel} ' = trendTmp;']);
%     
%     %absolute velocity
%     trendTmp = NaN(numSistersCat,3);
%     for iSister = 1 : numSistersCat
%         eval(['allValues = sisterVel' label{iLabel,1} '(iSister).observations(:,1);'])
%         allValues = abs(allValues);
%         indxAvail = find(~isnan(allValues));
%         [lineParam,S] = polyfit(indxAvail*timeLapse,allValues(indxAvail),1);
%         varCovMat = (inv(S.R)*inv(S.R)')*S.normr^2/S.df;
%         slopeStd = sqrt(varCovMat(1));
%         testStat = lineParam(1)/slopeStd;
%         if testStat > 0
%             pValue = 1 - tcdf(testStat,S.df);
%         else
%             pValue = tcdf(testStat,S.df);            
%         end
%         trendTmp(iSister,:) = [lineParam(1) slopeStd pValue];
%     end
%     eval(['sisterVelTrend' label{iLabel} ' = trendTmp;']);
%     
%     %angle with normal
%     trendTmp = NaN(numSistersCat,3);
%     for iSister = 1 : numSistersCat
%         eval(['allValues = angleNormal' label{iLabel,1} '(iSister).observations(:,1);'])
%         indxAvail = find(~isnan(allValues));
%         [lineParam,S] = polyfit(indxAvail*timeLapse,allValues(indxAvail),1);
%         varCovMat = (inv(S.R)*inv(S.R)')*S.normr^2/S.df;
%         slopeStd = sqrt(varCovMat(1));
%         testStat = lineParam(1)/slopeStd;
%         if testStat > 0
%             pValue = 1 - tcdf(testStat,S.df);
%         else
%             pValue = tcdf(testStat,S.df);            
%         end
%         trendTmp(iSister,:) = [lineParam(1) slopeStd pValue];
%     end
%     eval(['angleNormalTrend' label{iLabel} ' = trendTmp;']);
%     
%     %angular velocity
%     trendTmp = NaN(numSistersCat,3);
%     for iSister = 1 : numSistersCat
%         eval(['allValues = angularVel' label{iLabel,1} '(iSister).observations(:,1);'])
%         indxAvail = find(~isnan(allValues));
%         [lineParam,S] = polyfit(indxAvail*timeLapse,allValues(indxAvail),1);
%         varCovMat = (inv(S.R)*inv(S.R)')*S.normr^2/S.df;
%         slopeStd = sqrt(varCovMat(1));
%         testStat = lineParam(1)/slopeStd;
%         if testStat > 0
%             pValue = 1 - tcdf(testStat,S.df);
%         else
%             pValue = tcdf(testStat,S.df);            
%         end
%         trendTmp(iSister,:) = [lineParam(1) slopeStd pValue];
%     end
%     eval(['angularVelTrend' label{iLabel} ' = trendTmp;']);
%     
% end

    %     %% temporal trend stuff %%
    %
    %     for iLabel = goodLabel
    %
    %         %open figure and write title
    %         figure('Name',[fileName(1:end-4) ' - Temporal trends - ' label{iLabel,1}],'NumberTitle','off');
    %
    %         %create subplot 1
    %         subplot(4,2,1);
    %         hold on;
    %
    %         %plot histogram of distance temporal trend
    %         eval(['trend2plot = sisterDistTrend' label{iLabel,1} ';']);
    %         if ~isempty(trend2plot)
    %             [x] = histogram(trend2plot(:,1));
    %             hist(trend2plot(:,1),length(x));
    %         end
    %
    %         %write axes labels and title
    %         xlabel('Sister separation slopes (\mum/s)');
    %         ylabel('# of occurance');
    %
    %         %hold off figure
    %         hold off
    %
    %         %create subplot 2
    %         subplot(4,2,2);
    %         hold on;
    %
    %         %plot histogram of significant trends only
    %         trend2plot = trend2plot(trend2plot(:,3)<0.05,1);
    %         if ~isempty(trend2plot)
    %             [x] = histogram(trend2plot);
    %             hist(trend2plot,length(x));
    %         end
    %
    %         %write axes labels and title
    %         xlabel('Significant slopes only');
    %         ylabel('# of occurance');
    %
    %         %hold off figure
    %         hold off
    %
    %         %create subplot 3
    %         subplot(4,2,3);
    %         hold on;
    %
    %         %plot histogram of velocity temporal trend
    %         eval(['trend2plot = sisterVelTrend' label{iLabel,1} ';']);
    %         if ~isempty(trend2plot)
    %             [x] = histogram(trend2plot(:,1));
    %             hist(trend2plot(:,1),length(x));
    %         end
    %
    %         %write axes labels and title
    %         xlabel('Absolute rate change sister separation slopes (nm/s/s)');
    %         ylabel('# of occurances');
    %
    %         %hold off figure
    %         hold off
    %
    %         %create subplot 4
    %         subplot(4,2,4);
    %         hold on;
    %
    %         %plot histogram of significant trends only
    %         trend2plot = trend2plot(trend2plot(:,3)<0.05,1);
    %         if ~isempty(trend2plot)
    %             [x] = histogram(trend2plot);
    %             hist(trend2plot,length(x));
    %         end
    %
    %         %write axes labels and title
    %         xlabel('Significant slopes only');
    %         ylabel('# of occurances');
    %
    %         %hold off figure
    %         hold off
    %
    %         %create subplot 5
    %         subplot(4,2,5);
    %         hold on;
    %
    %         %plot histogram of angle with normal temporal trend
    %         eval(['trend2plot = angleNormalTrend' label{iLabel,1} ';']);
    %         if ~isempty(trend2plot)
    %             [x] = histogram(trend2plot(:,1));
    %             hist(trend2plot(:,1),length(x));
    %         end
    %
    %         %write axes labels and title
    %         xlabel('Angle with normal slopes (degrees/s)');
    %         ylabel('# of occurances');
    %
    %         %hold off figure
    %         hold off
    %
    %         %create subplot 6
    %         subplot(4,2,6);
    %         hold on;
    %
    %         %plot histogram of significant trends only
    %         trend2plot = trend2plot(trend2plot(:,3)<0.05,1);
    %         if ~isempty(trend2plot)
    %             [x] = histogram(trend2plot);
    %             hist(trend2plot,length(x));
    %         end
    %
    %         %write axes labels and title
    %         xlabel('Significant slopes only');
    %         ylabel('# of occurances');
    %
    %         %hold off figure
    %         hold off
    %
    %         %create subplot 7
    %         subplot(4,2,7);
    %         hold on;
    %
    %         %plot histogram of angular velocity temporal trend
    %         eval(['trend2plot = angularVelTrend' label{iLabel,1} ';']);
    %         if ~isempty(trend2plot)
    %             [x] = histogram(trend2plot(:,1));
    %             hist(trend2plot(:,1),length(x));
    %         end
    %
    %         %write axes labels and title
    %         xlabel('Angular velocity slopes  (degrees/s/s)');
    %         ylabel('# of occurance');
    %
    %         %hold off figure
    %         hold off
    %
    %         %create subplot 8
    %         subplot(4,2,8);
    %         hold on;
    %
    %         %plot histogram of angular velocity temporal trend
    %         trend2plot = trend2plot(trend2plot(:,3)<0.05,1);
    %         if ~isempty(trend2plot)
    %             [x] = histogram(trend2plot);
    %             hist(trend2plot,length(x));
    %         end
    %
    %         %write axes labels and title
    %         xlabel('Significant slopes only');
    %         ylabel('# of occurance');
    %
    %         %hold off figure
    %         hold off
    %
    %     end %(or iLabel = goodLabel)
    

