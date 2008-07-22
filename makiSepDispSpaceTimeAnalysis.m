function analysisStruct = makiSepDispSpaceTimeAnalysis(jobType,analysisStruct,verbose)
%MAKISEPDISPSPACETIMEANALYSIS analyzes the spatio-temporal dynamics of sister separation and oscillation
%
%SYNOPSIS analysisStruct = makiSepDispSpaceTimeAnalysis(jobType,analysisStruct,verbose)
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
%           .sepDispSpaceTime:
%
%REMARKS Code is not applicable to anaphase movies/frames
%
%Khuloud Jaqaman, July 2008

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

%get time between frames
timeLapse = round(2*dataStruct(1).dataProperties.timeLapse)/2;

%% collect sister separation and displacement information

%define cell array of sister labels
label{1,1} = 'Inlier';
label{2,1} = 'Unaligned';
label{3,1} = 'Lagging';

%initialize global sister index and structures
for iLabel = 1 : 3

    eval(['iGlobal' label{iLabel,1} ' = 0;'])

    eval(['separationSis12' label{iLabel,1} '(1:numSistersTot,1) = struct(''observations'',[],''time'',[]);'])
    eval(['separationChangeSis12' label{iLabel,1} '(1:numSistersTot,1) = struct(''observations'',[],''time'',[]);'])

    eval(['sepPChangeInt' label{iLabel,1} '(1:numSistersTot,1) = struct(''observations'',[],''time'',[]);'])
    eval(['sepNChangeInt' label{iLabel,1} '(1:numSistersTot,1) = struct(''observations'',[],''time'',[]);'])

    eval(['centerPositionSis12' label{iLabel,1} '(1:numSistersTot,1) = struct(''observations'',[],''time'',[]);'])
    eval(['centerPositionChangeSis12' label{iLabel,1} '(1:numSistersTot,1) = struct(''observations'',[],''time'',[]);'])

    eval(['centerPChangeInt' label{iLabel,1} '(1:numSistersTot,1) = struct(''observations'',[],''time'',[]);'])
    eval(['centerNChangeInt' label{iLabel,1} '(1:numSistersTot,1) = struct(''observations'',[],''time'',[]);'])

    eval(['sister1Disp' label{iLabel,1} '(1:numSistersTot,1) = struct(''observations'',[],''time'',[]);'])
    eval(['sister2Disp' label{iLabel,1} '(1:numSistersTot,1) = struct(''observations'',[],''time'',[]);'])
   
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
    
    if numSisters(iMovie) > 0

        %copy fields out of dataStruct(iMovie)
        sisterList = dataStruct(iMovie).sisterList;
        tracks = dataStruct(iMovie).tracks;
        frameAlignment = dataStruct(iMovie).frameAlignment;
        updatedClass = dataStruct(iMovie).updatedClass;
        planeFit = dataStruct(iMovie).planeFit;
        numFramesMovie = length(updatedClass);

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

            %%sort sisters to left and right
            
            %calculate the average coordinate of each sister along the normal
            %to the metaphase plate
            meanCoord1Normal = nanmean(coords1(:,1));
            meanCoord2Normal = nanmean(coords2(:,1));

            %put the sister with the smaller average coordinate on the "left"
            %negative is smaller than positive, no matter the magnitude
            if meanCoord2Normal < meanCoord1Normal
                tmp = coords2;
                coords2 = coords1;
                coords1 = tmp;
            end

            %%sister separation & its change

            %calculate vector connecting sisters
            sisterVec = coords2(:,1:3)-coords1(:,1:3); %um
            sisterVecVar = coords1(:,4:6).^2 + coords2(:,4:6).^2; %um^2

            %calculate sister separation
            sisterDist = sqrt(sum(sisterVec.^2,2)); %um
            sisterDistStd = sqrt( sum( sisterVecVar .* sisterVec.^2 ,2) ) ...
                ./ sisterDist;  %um
            
            %calculate change in sister separation between frames
            sisterDistChange = diff(sisterDist); %um
            sisterDistChangeStd = sqrt( sum( [sisterDistStd(2:end) ...
                sisterDistStd(1:end-1)].^2 ,2) ); %um
            
            %calculate intervals of +ve and -ve separation change
           
            %indicate intervals of +ve change by +1 and intervals of -ve change by -1
            changeSign = sign(sisterDistChange);
            
            %find transition points
            changeSignDiff = diff(changeSign);
            changeInterval = find(changeSignDiff~=0);
            changeIntValue = changeSignDiff(changeInterval);
            
            %go over transition points and calculate +ve and -ve change
            %periods
            sepPChangeInt = [];
            sepNChangeInt = [];
            for iInt = 1 : length(changeInterval)-1
                if ~isnan(changeIntValue(iInt))
                    if ~isnan(changeIntValue(iInt+1))
                        if changeIntValue(iInt) > 0
                            intCurrent = changeInterval(iInt+1) - changeInterval(iInt);
                            if intCurrent == 1
                                diffValue = -sisterDistChange(changeInterval(iInt)+1);
                                diffValueStd = sisterDistChangeStd(changeInterval(iInt)+1);
                                diffPValue = normcdf(diffValue,0,diffValueStd);
                                if diffPValue < 0.05
                                    sepPChangeInt = [sepPChangeInt; intCurrent]; %frames
                                end
                            else
                                sepPChangeInt = [sepPChangeInt; intCurrent]; %frames
                            end
                        else
                            intCurrent = changeInterval(iInt+1) - changeInterval(iInt);
                            if intCurrent == 1
                                diffValue = sisterDistChange(changeInterval(iInt)+1);
                                diffValueStd = sisterDistChangeStd(changeInterval(iInt)+1);
                                diffPValue = normcdf(diffValue,0,diffValueStd);
                                if diffPValue < 0.05
                                    sepNChangeInt = [sepNChangeInt; intCurrent]; %frames
                                end
                            else
                                sepNChangeInt = [sepNChangeInt; intCurrent]; %frames
                            end
                        end
                    end
                end
            end
            sepPChangeInt = sepPChangeInt * timeLapse; %s
            sepNChangeInt = sepNChangeInt * timeLapse; %s
           
            %%center position & its change
            
            %calculate position of sister center of mass along the normal
            %to the metaphase plate
            centerPos = mean([coords1(:,1) coords2(:,1)],2); %um
            centerPosStd = 0.5 * sqrt( sum( [coords1(:,4) coords2(:,4)].^2 ,2) ); %um

            %calculate change in position of center of mass between frames
            centerPosChange = diff(centerPos); %um
            centerPosChangeStd = sqrt( sum( [centerPosStd(2:end) ...
                centerPosStd(1:end-1)].^2 ,2) ); %um

            %calculate intervals of +ve and -ve displacement

            %calculate the p-values of these changes, assuming that they follow
            %a normal disitrubution with mean 0 and std given by the second
            %column in comPosChangeCurrent
            pValueChange = normcdf(centerPosChange,0,centerPosChangeStd);
            indxChange = find(pValueChange > 0.5);
            pValueChange(indxChange) = 1 - pValueChange(indxChange);

            %use an alpha-value of 0.05 to separate significant from
            %insignificant changes
            % 1 = significant movement in positive direction
            %-1 = insignificant movement in positive direction
            % 2 = significant movement in negative direction
            %-2 = insignificant movement in negative direction
            %NaN = missing because neighboring point is missing
            alpha = 0.05;
            trajStates = NaN(size(centerPosChange,1),1);
            trajStates(centerPosChange(:,1) >= 0 & pValueChange < alpha/2) =  1;
            trajStates(centerPosChange(:,1) >= 0 & pValueChange > alpha/2) = -1;
            trajStates(centerPosChange(:,1) <  0 & pValueChange < alpha/2) =  2;
            trajStates(centerPosChange(:,1) <  0 & pValueChange > alpha/2) = -2;

            %find intervals which have significant position change (positive or
            %negative) or which have unknown position change (due to missing
            %data)
            goodIntervals = find( trajStates~=-1 & trajStates~=-2 );

            %go over these intervals and assign whatever's in between its
            %proper state
            %conversion strategy:
            %**if either before or after is NaN, between is assigned NaN;
            %**if both before and after are +1/+2, between is assigned +1/+2;
            %**if before is +1 and after is +2, the point in between with
            %  maximum value is taken as the divider: intervals before it are
            %  assigned +1 and intervals after it are assigned +2;
            %**if before is +2 and after is -2, the point in between with
            %  minimum value is taken as the divider: intervals before it are
            %  assigned +2 and intervals after it are assigned +1;

            %assign intervals before the first good interval and intervals
            %after the last good interval a NaN
            trajStates(1:goodIntervals(1)-1) = NaN;
            trajStates(goodIntervals(end)+1:end) = NaN;
            
            %go over intervals in between good intervals
            for iInt = 1 : length(goodIntervals)-1

                %get states before and after
                stateBef = trajStates(goodIntervals(iInt));
                stateAft = trajStates(goodIntervals(iInt+1));

                %determine what to assign in between based on the states before
                %and after
                if isnan(stateBef) || isnan(stateAft) %either is NaN

                    trajStates(goodIntervals(iInt)+1:goodIntervals(iInt+1)-1) = NaN;

                elseif stateBef == stateAft %both are the same

                    trajStates(goodIntervals(iInt)+1:goodIntervals(iInt+1)-1) = stateBef;

                elseif stateBef == 1 %transition point from pos. to neg. movement

                    %get position values in this undetermined interval
                    trajValues = centerPos(goodIntervals(iInt)+1:...
                        goodIntervals(iInt+1),1);

                    %find the point with maximum value
                    maxPoint = find(trajValues==max(trajValues));

                    %give all intervals before this point a state of +1
                    trajStates(goodIntervals(iInt)+1:goodIntervals(iInt)+maxPoint-1) = 1;

                    %give all intervals after his point a state of +2;
                    trajStates(goodIntervals(iInt)+maxPoint:goodIntervals(iInt+1)-1) = 2;

                else %transition point from neg. to pos. movement

                    %get position values in this undetermined interval
                    trajValues = centerPos(goodIntervals(iInt)+1:...
                        goodIntervals(iInt+1),1);

                    %find the point with minimum value
                    minPoint = find(trajValues==min(trajValues));

                    %give all intervals before this point a state of +2
                    trajStates(goodIntervals(iInt)+1:goodIntervals(iInt)+minPoint-1) = 2;

                    %give all intervals after this point a state of +1
                    trajStates(goodIntervals(iInt)+minPoint:goodIntervals(iInt+1)-1) = 1;
                    
                end

            end

            %indicate intervals of +ve change by +1 and intervals of -ve change by -1
            changeSign = trajStates;
            changeSign(changeSign==2) = -1;
            
            %find transition points
            changeSignDiff = diff(changeSign);
            changeInterval = find(changeSignDiff~=0);
            changeIntValue = changeSignDiff(changeInterval);
            
            %go over transition points and calculate +ve and -ve change
            %periods
            centerPChangeInt = [];
            centerNChangeInt = [];
            for iInt = 1 : length(changeInterval)-1
                if ~isnan(changeIntValue(iInt))
                    if ~isnan(changeIntValue(iInt+1))
                        if changeIntValue(iInt) > 0
                            intCurrent = changeInterval(iInt+1) - changeInterval(iInt);
                            centerPChangeInt = [centerPChangeInt; intCurrent]; %frames
                        else
                            intCurrent = changeInterval(iInt+1) - changeInterval(iInt);
                            centerNChangeInt = [centerNChangeInt; intCurrent]; %frames
                        end
                    end
                end
            end
            centerPChangeInt = centerPChangeInt * timeLapse; %s
            centerNChangeInt = centerNChangeInt * timeLapse; %s

            %%sister displacement between frames
            
            %calculate sister displacement between frames along the normal
            coords1Diff = coords1(2:end,1) - coords1(1:end-1,1); %um
            coords2Diff = coords2(2:end,1) - coords2(1:end-1,1); %um

            %store sister information
            iLabel = sisterType(iSister) + 1;

            %increase global index of sister type by 1
            eval(['iGlobal' label{iLabel,1} ' = iGlobal' label{iLabel,1} ' + 1;'])

            %store information
            eval(['separationSis12' label{iLabel,1} '(iGlobal' label{iLabel,1} ').observations = [sisterDist sisterDistStd];']) %um
            eval(['separationChangeSis12' label{iLabel,1} '(iGlobal' label{iLabel,1} ').observations = [sisterDistChange sisterDistChangeStd];']) %um
            eval(['sepPChangeInt' label{iLabel,1} '(iGlobal' label{iLabel,1} ').observations = sepPChangeInt;']) %s
            eval(['sepNChangeInt' label{iLabel,1} '(iGlobal' label{iLabel,1} ').observations = sepNChangeInt;']) %s
            eval(['centerPositionSis12' label{iLabel,1} '(iGlobal' label{iLabel,1} ').observations = [centerPos centerPosStd];']) %um
            eval(['centerPositionChangeSis12' label{iLabel,1} '(iGlobal' label{iLabel,1} ').observations = [centerPosChange centerPosChangeStd];']) %um
            eval(['centerPChangeInt' label{iLabel,1} '(iGlobal' label{iLabel,1} ').observations = centerPChangeInt;']) %s
            eval(['centerNChangeInt' label{iLabel,1} '(iGlobal' label{iLabel,1} ').observations = centerNChangeInt;']) %s
            eval(['sister1Disp' label{iLabel,1} '(iGlobal' label{iLabel,1} ').observations = coords1Diff;']) %um
            eval(['sister2Disp' label{iLabel,1} '(iGlobal' label{iLabel,1} ').observations = coords2Diff;']) %um
           
        end %(for iSister = 1 : numSisters(iMovie) )

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

    eval(['separationSis12' label{iLabel,1} ' = separationSis12' ...
        label{iLabel,1} '(1:iGlobal' label{iLabel,1} ');'])
    eval(['separationChangeSis12' label{iLabel,1} ' = separationChangeSis12' ...
        label{iLabel,1} '(1:iGlobal' label{iLabel,1} ');'])
    eval(['sepPChangeInt' label{iLabel,1} ' = sepPChangeInt' ...
        label{iLabel,1} '(1:iGlobal' label{iLabel,1} ');'])
    eval(['sepNChangeInt' label{iLabel,1} ' = sepNChangeInt' ...
        label{iLabel,1} '(1:iGlobal' label{iLabel,1} ');'])
    eval(['centerPositionSis12' label{iLabel,1} ' = centerPositionSis12' ...
        label{iLabel,1} '(1:iGlobal' label{iLabel,1} ');'])
    eval(['centerPositionChangeSis12' label{iLabel,1} ' = centerPositionChangeSis12' ...
        label{iLabel,1} '(1:iGlobal' label{iLabel,1} ');'])
    eval(['centerPChangeInt' label{iLabel,1} ' = centerPChangeInt' ...
        label{iLabel,1} '(1:iGlobal' label{iLabel,1} ');'])
    eval(['centerNChangeInt' label{iLabel,1} ' = centerNChangeInt' ...
        label{iLabel,1} '(1:iGlobal' label{iLabel,1} ');'])
    eval(['sister1Disp' label{iLabel,1} ' = sister1Disp' ...
        label{iLabel,1} '(1:iGlobal' label{iLabel,1} ');'])
    eval(['sister2Disp' label{iLabel,1} ' = sister2Disp' ...
        label{iLabel,1} '(1:iGlobal' label{iLabel,1} ');'])

end

%define maximum lag
maxLag = 20;

%% sister separation analysis

%initialization
for iLabel = 1 : 3
    
    eval(['separationDistr' label{iLabel,1} ' = [];'])
    eval(['sepChangeDistr' label{iLabel,1} ' = [];'])

    eval(['separationParam' label{iLabel,1} ' = [];'])
    eval(['sepChangeParam' label{iLabel,1} ' = [];'])

    eval(['separationIndParam' label{iLabel,1} ' = [];'])
    eval(['sepChangeIndParam' label{iLabel,1} ' = [];'])
    
    eval(['sepChangeAutocorr' label{iLabel,1} ' = NaN(maxLag+1,2);'])
    eval(['sepChangeIndAutocorr' label{iLabel,1} ' = NaN(maxLag+1,2,numMovies);'])
    
    eval(['sepPChangeIntervalDistr' label{iLabel,1} ' = [];'])
    eval(['sepNChangeIntervalDistr' label{iLabel,1} ' = [];'])
    
    eval(['sepPChangeIntervalParam' label{iLabel,1} ' = [];'])
    eval(['sepNChangeIntervalParam' label{iLabel,1} ' = [];'])
   
    eval(['sepPChangeIntervalIndParam' label{iLabel,1} ' = [];'])
    eval(['sepNChangeIntervalIndParam' label{iLabel,1} ' = [];'])
   
end

%calculation
for iLabel = goodLabel

    % sister separation distribution %

    %overall
    eval(['allValues = vertcat(separationSis12' label{iLabel,1} '.observations);']);
    if ~isempty(allValues)
        allValues = allValues(:,1);
        allValues = allValues(~isnan(allValues));
        eval(['separationDistr' label{iLabel,1} ' = allValues;']);
        if ~isempty(allValues)
            eval(['separationParam' label{iLabel,1} ' = [mean(allValues) ' ...
                'std(allValues) min(allValues) prctile(allValues,25) ' ...
                'prctile(allValues,50) prctile(allValues,75) max(allValues)];']);
        end
    end

    %individual cells
    for iMovie = 1 : numMovies
        eval(['allValues = vertcat(separationSis12' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie)).observations);']);
        if ~isempty(allValues)
            allValues = allValues(:,1);
            allValues = allValues(~isnan(allValues));
            if ~isempty(allValues)
                eval(['separationIndParam' label{iLabel,1} '(iMovie,:)' ...
                    ' = [mean(allValues) std(allValues) min(allValues) ' ...
                    'prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) '...
                    'max(allValues)];']);
            end
        end
    end

    % sister separation change distribution %
    
    %overall
    eval(['allValues = vertcat(separationChangeSis12' label{iLabel,1} '.observations);']);
    if ~isempty(allValues)
        allValues = allValues(:,1);
        allValues = allValues(~isnan(allValues));
        eval(['sepChangeDistr' label{iLabel,1} ' = allValues;']);
        if ~isempty(allValues)
            eval(['sepChangeParam' label{iLabel,1} ' = [mean(allValues) ' ...
                'std(allValues) min(allValues) prctile(allValues,25) ' ...
                'prctile(allValues,50) prctile(allValues,75) max(allValues)];']);
        end
    end

    %individual cells
    for iMovie = 1 : numMovies
        eval(['allValues = vertcat(separationChangeSis12' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie)).observations);']);
        if ~isempty(allValues)
            allValues = allValues(:,1);
            allValues = allValues(~isnan(allValues));
            if ~isempty(allValues)
                eval(['sepChangeIndParam' label{iLabel,1} '(iMovie,:)' ...
                    ' = [mean(allValues) std(allValues) min(allValues) ' ...
                    'prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) '...
                    'max(allValues)];']);
            end
        end
    end

    % separation change autocorrelation %
    
    %overall
    eval(['[tmpCorr,errFlag] = autoCorr(separationChangeSis12' label{iLabel,1} ',maxLag);'])
    if ~errFlag
        eval(['sepChangeAutocorr' label{iLabel,1} ' = tmpCorr;'])
    end

    %individual cells
    for iMovie = 1 : numMovies
        eval(['traj = separationChangeSis12' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie));'])
        if ~isempty(traj)
            [tmpCorr,errFlag] = autoCorr(traj,maxLag);
            if ~errFlag
                eval(['sepChangeIndAutocorr' label{iLabel,1} '(:,:,iMovie) = tmpCorr;'])
            end
        end
    end

    % positive separation change intervals %
    
    %overall
    eval(['allValues = vertcat(sepPChangeInt' label{iLabel,1} '.observations);']);
    if ~isempty(allValues)
        allValues = allValues(:,1);
        allValues = allValues(~isnan(allValues));
        eval(['sepPChangeIntervalDistr' label{iLabel,1} ' = allValues;']);
        if ~isempty(allValues)
            eval(['sepPChangeIntervalParam' label{iLabel,1} ' = [mean(allValues) ' ...
                'std(allValues) min(allValues) prctile(allValues,25) ' ...
                'prctile(allValues,50) prctile(allValues,75) max(allValues)];']);
        end
    end

    %individual cells
    for iMovie = 1 : numMovies
        eval(['allValues = vertcat(sepPChangeInt' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie)).observations);']);
        if ~isempty(allValues)
            allValues = allValues(:,1);
            allValues = allValues(~isnan(allValues));
            if ~isempty(allValues)
                eval(['sepPChangeIntervalIndParam' label{iLabel,1} '(iMovie,:)' ...
                    ' = [mean(allValues) std(allValues) min(allValues) ' ...
                    'prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) '...
                    'max(allValues)];']);
            end
        end
    end

    % negative separation change intervals %
    
    %overall
    eval(['allValues = vertcat(sepNChangeInt' label{iLabel,1} '.observations);']);
    if ~isempty(allValues)
        allValues = allValues(:,1);
        allValues = allValues(~isnan(allValues));
        eval(['sepNChangeIntervalDistr' label{iLabel,1} ' = allValues;']);
        if ~isempty(allValues)
            eval(['sepNChangeIntervalParam' label{iLabel,1} ' = [mean(allValues) ' ...
                'std(allValues) min(allValues) prctile(allValues,25) ' ...
                'prctile(allValues,50) prctile(allValues,75) max(allValues)];']);
        end
    end

    %individual cells
    for iMovie = 1 : numMovies
        eval(['allValues = vertcat(sepNChangeInt' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie)).observations);']);
        if ~isempty(allValues)
            allValues = allValues(:,1);
            allValues = allValues(~isnan(allValues));
            if ~isempty(allValues)
                eval(['sepNChangeIntervalIndParam' label{iLabel,1} '(iMovie,:)' ...
                    ' = [mean(allValues) std(allValues) min(allValues) ' ...
                    'prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) '...
                    'max(allValues)];']);
            end
        end
    end
    
end %(for iLabel = goodLabel)

%% center of mass motion analysis

%initialization
for iLabel = 1 : 3
    
    eval(['centerPosDistr' label{iLabel,1} ' = [];'])
    eval(['centerPosChangeDistr' label{iLabel,1} ' = [];'])

    eval(['centerPosParam' label{iLabel,1} ' = [];'])
    eval(['centerPosChangeParam' label{iLabel,1} ' = [];'])

    eval(['centerPosIndParam' label{iLabel,1} ' = [];'])
    eval(['centerPosChangeIndParam' label{iLabel,1} ' = [];'])
    
    eval(['centerPosChangeAutocorr' label{iLabel,1} ' = NaN(maxLag+1,2);'])
    eval(['centerPosChangeIndAutocorr' label{iLabel,1} ' = NaN(maxLag+1,2,numMovies);'])
    
    eval(['centerPosPChangeIntervalDistr' label{iLabel,1} ' = [];'])
    eval(['centerPosNChangeIntervalDistr' label{iLabel,1} ' = [];'])
    
    eval(['centerPosPChangeIntervalParam' label{iLabel,1} ' = [];'])
    eval(['centerPosNChangeIntervalParam' label{iLabel,1} ' = [];'])
   
    eval(['centerPosPChangeIntervalIndParam' label{iLabel,1} ' = [];'])
    eval(['centerPosNChangeIntervalIndParam' label{iLabel,1} ' = [];'])
   
end

%calculation
for iLabel = goodLabel

    % center position distribution %

    %overall
    eval(['allValues = vertcat(centerPositionSis12' label{iLabel,1} '.observations);']);
    if ~isempty(allValues)
        allValues = allValues(:,1);
        allValues = allValues(~isnan(allValues));
        eval(['centerPosDistr' label{iLabel,1} ' = allValues;']);
        if ~isempty(allValues)
            eval(['centerPosParam' label{iLabel,1} ' = [mean(allValues) ' ...
                'std(allValues) min(allValues) prctile(allValues,25) ' ...
                'prctile(allValues,50) prctile(allValues,75) max(allValues)];']);
        end
    end

    %individual cells
    for iMovie = 1 : numMovies
        eval(['allValues = vertcat(centerPositionSis12' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie)).observations);']);
        if ~isempty(allValues)
            allValues = allValues(:,1);
            allValues = allValues(~isnan(allValues));
            if ~isempty(allValues)
                eval(['centerPosIndParam' label{iLabel,1} '(iMovie,:)' ...
                    ' = [mean(allValues) std(allValues) min(allValues) ' ...
                    'prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) '...
                    'max(allValues)];']);
            end
        end
    end

    % center position change distribution %
    
    %overall
    eval(['allValues = vertcat(centerPositionChangeSis12' label{iLabel,1} '.observations);']);
    if ~isempty(allValues)
        allValues = allValues(:,1);
        allValues = allValues(~isnan(allValues));
        eval(['centerPosChangeDistr' label{iLabel,1} ' = allValues;']);
        if ~isempty(allValues)
            eval(['centerPosChangeParam' label{iLabel,1} ' = [mean(allValues) ' ...
                'std(allValues) min(allValues) prctile(allValues,25) ' ...
                'prctile(allValues,50) prctile(allValues,75) max(allValues)];']);
        end
    end

    %individual cells
    for iMovie = 1 : numMovies
        eval(['allValues = vertcat(centerPositionChangeSis12' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie)).observations);']);
        if ~isempty(allValues)
            allValues = allValues(:,1);
            allValues = allValues(~isnan(allValues));
            if ~isempty(allValues)
                eval(['centerPosChangeIndParam' label{iLabel,1} '(iMovie,:)' ...
                    ' = [mean(allValues) std(allValues) min(allValues) ' ...
                    'prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) '...
                    'max(allValues)];']);
            end
        end
    end

    % center position change autocorrelation %
    
    %overall
    eval(['[tmpCorr,errFlag] = autoCorr(centerPositionChangeSis12' label{iLabel,1} ',maxLag);'])
    if ~errFlag
        eval(['centerPosChangeAutocorr' label{iLabel,1} ' = tmpCorr;'])
    end

    %individual cells
    for iMovie = 1 : numMovies
        eval(['traj = centerPositionChangeSis12' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie));'])
        if ~isempty(traj)
            [tmpCorr,errFlag] = autoCorr(traj,maxLag);
            if ~errFlag
                eval(['centerPosChangeIndAutocorr' label{iLabel,1} '(:,:,iMovie) = tmpCorr;'])
            end
        end
    end

    % positive center position change intervals %
    
    %overall
    eval(['allValues = vertcat(centerPChangeInt' label{iLabel,1} '.observations);']);
    if ~isempty(allValues)
        allValues = allValues(:,1);
        allValues = allValues(~isnan(allValues));
        eval(['centerPosPChangeIntervalDistr' label{iLabel,1} ' = allValues;']);
        if ~isempty(allValues)
            eval(['centerPosPChangeIntervalParam' label{iLabel,1} ' = [mean(allValues) ' ...
                'std(allValues) min(allValues) prctile(allValues,25) ' ...
                'prctile(allValues,50) prctile(allValues,75) max(allValues)];']);
        end
    end

    %individual cells
    for iMovie = 1 : numMovies
        eval(['allValues = vertcat(centerPChangeInt' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie)).observations);']);
        if ~isempty(allValues)
            allValues = allValues(:,1);
            allValues = allValues(~isnan(allValues));
            if ~isempty(allValues)
                eval(['centerPosPChangeIntervalIndParam' label{iLabel,1} '(iMovie,:)' ...
                    ' = [mean(allValues) std(allValues) min(allValues) ' ...
                    'prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) '...
                    'max(allValues)];']);
            end
        end
    end

    % negative separation change intervals %
    
    %overall
    eval(['allValues = vertcat(centerNChangeInt' label{iLabel,1} '.observations);']);
    if ~isempty(allValues)
        allValues = allValues(:,1);
        allValues = allValues(~isnan(allValues));
        eval(['centerPosNChangeIntervalDistr' label{iLabel,1} ' = allValues;']);
        if ~isempty(allValues)
            eval(['centerPosNChangeIntervalParam' label{iLabel,1} ' = [mean(allValues) ' ...
                'std(allValues) min(allValues) prctile(allValues,25) ' ...
                'prctile(allValues,50) prctile(allValues,75) max(allValues)];']);
        end
    end

    %individual cells
    for iMovie = 1 : numMovies
        eval(['allValues = vertcat(centerNChangeInt' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie)).observations);']);
        if ~isempty(allValues)
            allValues = allValues(:,1);
            allValues = allValues(~isnan(allValues));
            if ~isempty(allValues)
                eval(['centerPosNChangeIntervalIndParam' label{iLabel,1} '(iMovie,:)' ...
                    ' = [mean(allValues) std(allValues) min(allValues) ' ...
                    'prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) '...
                    'max(allValues)];']);
            end
        end
    end

end %(for iLabel = goodLabel)

%% sister displacement coupling analysis

%initialization
for iLabel = 1 : 3
    
    eval(['dispCrosscorr' label{iLabel,1} ' = NaN(1,2);'])
    eval(['dispIndCrosscorr' label{iLabel,1} ' = NaN(numMovies,2);'])
   
end

for iLabel = goodLabel

    %overall
    eval(['[tmpCorr,errFlag] = crossCorr(sister1Disp' label{iLabel,1} ...
        ',sister2Disp' label{iLabel,1} ',0);'])
    if ~errFlag
        eval(['dispCrosscorr' label{iLabel,1} ' = tmpCorr;'])
    end
    
    %individual cells
    for iMovie = 1 : numMovies
        eval(['traj1 = sister1Disp' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie));'])
        eval(['traj2 = sister2Disp' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie));'])
        if ~isempty(traj1)
            [tmpCorr,errFlag] = crossCorr(traj1,traj2,0);
            if ~errFlag
                eval(['dispIndCrosscorr' label{iLabel,1} '(iMovie,:) = tmpCorr;'])
            end
        end
    end

end %(for iLabel = goodLabel)

%% output to analysisStruct

for iLabel = 1 : 3

    eval(['distribution = struct('...
        '''sisterSeparation'',separationDistr' label{iLabel,1} ','...
        '''sisterSepChange'',sepChangeDistr' label{iLabel,1} ','...
        '''sepPChangeInterval'',sepPChangeIntervalDistr' label{iLabel,1} ','...
        '''sepNChangeInterval'',sepNChangeIntervalDistr' label{iLabel,1} ','...
        '''centerPosition'',centerPosDistr' label{iLabel,1} ','...
        '''centerPosChange'',centerPosChangeDistr' label{iLabel,1} ','...
        '''centerPosPChangeInterval'',centerPosPChangeIntervalDistr' label{iLabel,1} ','...
        '''centerPosNChangeInterval'',centerPosNChangeIntervalDistr' label{iLabel,1} ');']);
    
    eval(['meanStdMin25P50P75PMax.all = struct('...
        '''sisterSeparation'',separationParam' label{iLabel,1} ','...
        '''sisterSepChange'',sepChangeParam' label{iLabel,1} ','...
        '''sepPChangeInterval'',sepPChangeIntervalParam' label{iLabel,1} ','...
        '''sepNChangeInterval'',sepNChangeIntervalParam' label{iLabel,1} ','...
        '''centerPosition'',centerPosParam' label{iLabel,1} ','...
        '''centerPosChange'',centerPosChangeParam' label{iLabel,1} ','...
        '''centerPosPChangeInterval'',centerPosPChangeIntervalParam' label{iLabel,1} ','...
        '''centerPosNChangeInterval'',centerPosNChangeIntervalParam' label{iLabel,1} ');']);
    
    eval(['meanStdMin25P50P75PMax.indcell = struct('...
        '''sisterSeparation'',separationIndParam' label{iLabel,1} ','...
        '''sisterSepChange'',sepChangeIndParam' label{iLabel,1} ','...
        '''sepPChangeInterval'',sepPChangeIntervalIndParam' label{iLabel,1} ','...
        '''sepNChangeInterval'',sepNChangeIntervalIndParam' label{iLabel,1} ','...
        '''centerPosition'',centerPosIndParam' label{iLabel,1} ','...
        '''centerPosChange'',centerPosChangeIndParam' label{iLabel,1} ','...
        '''centerPosPChangeInterval'',centerPosPChangeIntervalIndParam' label{iLabel,1} ','...
        '''centerPosNChangeInterval'',centerPosNChangeIntervalIndParam' label{iLabel,1} ');']);
    
    eval(['autocorr.all = struct('...
        '''sisterSepChange'',sepChangeAutocorr' label{iLabel,1} ','...
        '''centerPosChange'',centerPosChangeAutocorr' label{iLabel,1} ');']);

    eval(['autocorr.indcell = struct('...
        '''sisterSepChange'',sepChangeIndAutocorr' label{iLabel,1} ','...
        '''centerPosChange'',centerPosChangeIndAutocorr' label{iLabel,1} ');']);

    eval(['crosscorr.all = struct('...
        '''sisterDisplacement'',dispCrosscorr' label{iLabel,1} ');']);

    eval(['crosscorr.indcell = struct('...
        '''sisterDisplacement'',dispIndCrosscorr' label{iLabel,1} ');']);

    eval(['numSistersCat = iGlobal' label{iLabel,1} ';']);

    eval([label{iLabel,1} ' = struct('...
        '''numSisters'',numSistersCat,'...
        '''distribution'',distribution,'...
        '''meanStdMin25P50P75PMax'',meanStdMin25P50P75PMax,'...
        '''autocorr'',autocorr,'...
        '''crosscorr'',crosscorr);'])

end %(for iLabel = 1 : 3)

inputParam = [];
sepDispSpaceTime = struct('Inlier',Inlier,'Unaligned',Unaligned,...
    'Lagging',Lagging,'inputParam',inputParam);

%check whether current analysisStruct already has the sepDispSpaceTime field
fieldExists = isfield(analysisStruct,'sepDispSpaceTime');

%store results in analysisStruct
analysisStruct.sepDispSpaceTime = sepDispSpaceTime;

%if sepDispSpaceTime field already existed, add 1 to the version number in
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
    
    %%sister separation%%

    for iLabel = goodLabel

        %open figure and write title
        figFileName = [fileName(1:end-4) '-SisterSeparation-' label{iLabel,1}];
        figHandle = figure('Name',figFileName,'NumberTitle','off');

        %create subplot 1
        subplot(3,2,1);
        hold on

        %histogram of sister separation
        eval(['tmpVar = separationDistr' label{iLabel,1} ';'])
        
        if ~isempty(tmpVar)

            n = histogram(tmpVar,[],0);
            histogram(tmpVar,[],0);

            %write axes labels
            xlabel('Sister separation (\mum)');
            ylabel('# of occurances');

            %add average and std
            eval(['tmpVar = separationParam' label{iLabel,1} '(1:2);'])
            text(tmpVar(1)+tmpVar(2),max(n),sprintf('mean+-std:\n%4.2f +- %4.2f',tmpVar(1),tmpVar(2)));
            
        end

        %hold off subplot 1
        hold off

        %create subplot 2
        subplot(3,2,2);
        hold on

        %histogram of sister separation change
        eval(['tmpVar = sepChangeDistr' label{iLabel,1} ';'])
        
        if ~isempty(tmpVar)

            n = histogram(tmpVar,[],0);
            histogram(tmpVar,[],0);

            %write axes labels
            xlabel('Frame-to-frame sister separation change (\mum)');
            ylabel('# of occurances');

            %add average and std
            eval(['tmpVar = sepChangeParam' label{iLabel,1} '(1:2);'])
            text(tmpVar(1)+tmpVar(2),max(n),sprintf('mean+-std:\n%4.2f +- %4.2f',tmpVar(1),tmpVar(2)));
            
        end

        %hold off subplot 2
        hold off
        
        %create subplot 3
        subplot(3,2,3);
        hold on

        %histogram of intervals of positive change in sister separation
        eval(['tmpVar = sepPChangeIntervalDistr' label{iLabel,1} ';']);
        
        if ~isempty(tmpVar)

            n = hist(tmpVar,(timeLapse:timeLapse:max(tmpVar)+timeLapse));
            hist(tmpVar,(timeLapse:timeLapse:max(tmpVar)+timeLapse));

            %write axes labels
            xlabel('Positive change intervals (s)');
            ylabel('# of occurances');

            %add average and std
            eval(['tmpVar = sepPChangeIntervalParam' label{iLabel,1} '(1:2);'])
            text(tmpVar(1)+tmpVar(2),max(n),sprintf('mean+-std:\n%4.2f +- %4.2f',tmpVar(1),tmpVar(2)));
            
        end

        %hold off subplot 3
        hold off
        
        %create subplot 4
        subplot(3,2,4);
        hold on

        %histogram of intervals of negative change in sister separation
        eval(['tmpVar = sepNChangeIntervalDistr' label{iLabel,1} ';']);
        
        if ~isempty(tmpVar)

            n = hist(tmpVar,(timeLapse:timeLapse:max(tmpVar)+timeLapse));
            hist(tmpVar,(timeLapse:timeLapse:max(tmpVar)+timeLapse));

            %write axes labels
            xlabel('Negative change intervals (s)');
            ylabel('# of occurances');

            %add average and std
            eval(['tmpVar = sepNChangeIntervalParam' label{iLabel,1} '(1:2);'])
            text(tmpVar(1)+tmpVar(2),max(n),sprintf('mean+-std:\n%4.2f +- %4.2f',tmpVar(1),tmpVar(2)));
            
        end

        %hold off subplot 4
        hold off
        
        %create subplot 5
        subplot(3,2,5);
        hold on

        %autocorrelation of separation change - overall
        eval(['plot((0:maxLag)*timeLapse,sepChangeAutocorr' label{iLabel,1} '(:,1),''marker'',''.'');']);
        plot([0 maxLag]*timeLapse,[0 0],'k--')
        
        %write axes labels
        xlabel('Time lag (s)');
        ylabel('Separation change autocorrelation - overall');

        %hold off subplot 5
        hold off
        
        %create subplot 6
        subplot(3,2,6);
        hold on

        %autocorrelation of separation change - individual cells
        eval(['tmpVar = squeeze(sepChangeIndAutocorr' label{iLabel,1} '(:,1,:));'])
        plot((0:maxLag)*timeLapse,tmpVar,'marker','.');
        plot([0 maxLag]*timeLapse,[0 0],'k--')

        %write axes labels
        xlabel('Time lag (s)');
        ylabel('Separation change autocorrelation - individual');

        %hold off subplot 6
        hold off
        
        %save figure
        saveas(figHandle,fullfile(dir2SaveRes,figFileName),'fig');

    end %(for iLabel = goodLabel)

    %%center position%%

    for iLabel = goodLabel

        %open figure and write title
        figFileName = [fileName(1:end-4) '-SisterCenterPosition-' label{iLabel,1}];
        figHandle = figure('Name',figFileName,'NumberTitle','off');

        %create subplot 1
        subplot(3,2,1);
        hold on
        
        %histogram of center position
        eval(['tmpVar = centerPosDistr' label{iLabel,1} ';'])

        if ~isempty(tmpVar)

            n = histogram(tmpVar,[],0);
            histogram(tmpVar,[],0);

            %write axes labels
            xlabel('Sister center position (\mum)');
            ylabel('# of occurances');

            %add average and std
            eval(['tmpVar = centerPosParam' label{iLabel,1} '(1:2);'])
            text(tmpVar(1)+tmpVar(2),max(n),sprintf('mean+-std:\n%4.2f +- %4.2f',tmpVar(1),tmpVar(2)));

        end

        %hold off subplot 1
        hold off

        %create subplot 2
        subplot(3,2,2);
        hold on

        %histogram of sister center displacement
        eval(['tmpVar = centerPosChangeDistr' label{iLabel,1} ';'])
        
        if ~isempty(tmpVar)

            n = histogram(tmpVar,[],0);
            histogram(tmpVar,[],0);

            %write axes labels
            xlabel('Frame-to-frame sister center displacement (\mum)');
            ylabel('# of occurances');

            %add average and std
            eval(['tmpVar = centerPosChangeParam' label{iLabel,1} '(1:2);'])
            text(tmpVar(1)+tmpVar(2),max(n),sprintf('mean+-std:\n%4.2f +- %4.2f',tmpVar(1),tmpVar(2)));

        end

        %hold off subplot 2
        hold off
        
        %create subplot 3
        subplot(3,2,3);
        hold on

        %histogram of intervals of positive change in center position
        eval(['tmpVar = centerPosPChangeIntervalDistr' label{iLabel,1} ';']);

        if ~isempty(tmpVar)

            n = hist(tmpVar,(timeLapse:timeLapse:max(tmpVar)+timeLapse));
            hist(tmpVar,(timeLapse:timeLapse:max(tmpVar)+timeLapse));

            %write axes labels
            xlabel('Positive displacement intervals (s)');
            ylabel('# of occurances');

            %add average and std
            eval(['tmpVar = centerPosPChangeIntervalParam' label{iLabel,1} '(1:2);'])
            text(tmpVar(1)+tmpVar(2),max(n),sprintf('mean+-std:\n%4.2f +- %4.2f',tmpVar(1),tmpVar(2)));

        end

        %hold off subplot 3
        hold off
        
        %create subplot 4
        subplot(3,2,4);
        hold on

        %histogram of intervals of negative change in center position
        eval(['tmpVar = centerPosNChangeIntervalDistr' label{iLabel,1} ';']);
        
        if ~isempty(tmpVar)

            n = hist(tmpVar,(timeLapse:timeLapse:max(tmpVar)+timeLapse));
            hist(tmpVar,(timeLapse:timeLapse:max(tmpVar)+timeLapse));

            %write axes labels
            xlabel('Negative displacement intervals (s)');
            ylabel('# of occurances');

            %add average and std
            eval(['tmpVar = centerPosNChangeIntervalParam' label{iLabel,1} '(1:2);'])
            text(tmpVar(1)+tmpVar(2),max(n),sprintf('mean+-std:\n%4.2f +- %4.2f',tmpVar(1),tmpVar(2)));
            
        end

        %hold off subplot 4
        hold off
        
        %create subplot 5
        subplot(3,2,5);
        hold on

        %autocorrelation of center position change - overall
        eval(['plot((0:maxLag)*timeLapse,centerPosChangeAutocorr' label{iLabel,1} '(:,1),''marker'',''.'');']);
        plot([0 maxLag]*timeLapse,[0 0],'k--')
        
        %write axes labels
        xlabel('Time lag (s)');
        ylabel('Sister center displacement autocorrelation - overall');

        %add cross-correlation information
        eval(['tmpVar = dispCrosscorr' label{iLabel,1} '(1:2);'])
        text(maxLag*timeLapse/2,0.7,sprintf('crosscorr:\n%4.2f +- %4.2f',tmpVar(1),tmpVar(2)));

        %hold off subplot 5
        hold off
        
        %create subplot 6
        subplot(3,2,6);
        hold on

        %autocorrelation of center displacement - individual cells
        eval(['tmpVar = squeeze(centerPosChangeIndAutocorr' label{iLabel,1} '(:,1,:));'])
        plot((0:maxLag)*timeLapse,tmpVar,'marker','.');
        plot([0 maxLag]*timeLapse,[0 0],'k--')

        %write axes labels
        xlabel('Time lag (s)');
        ylabel('Sister center displacement autocorrelation - individual');

        %hold off subplot 6
        hold off
        
        %save figure
        saveas(figHandle,fullfile(dir2SaveRes,figFileName),'fig');

    end %(for iLabel = goodLabel)

end




