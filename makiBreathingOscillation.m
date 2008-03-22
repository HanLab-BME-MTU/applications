function analysisStruct = makiBreathingOscillation(jobType,analysisStruct,verbose)
%MAKIBREATHINGOSCILLATION looks for coupling between sister breathing and sister oscillation normal to the metaphase plate
%
%SYNOPSIS analysisStruct = makiBreathingOscillation(jobType,analysisStruct,verbose)
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
%           .breathingOscillationCoupling: 
%
%REMARKS Code is not applicable to anaphase movies/frames
%
%Khuloud Jaqaman, March 2008

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

%get time between frames
timeLapse = round(dataStruct(1).dataProperties.timeLapse);

%% collect sister breathing and oscillation information

%define cell array of sister labels
label{1,1} = 'Inlier';
label{2,1} = 'Unaligned';
label{3,1} = 'Lagging';

%initialize global sister index and structures
for iLabel = 1 : 3

    eval(['iGlobal' label{iLabel,1} ' = 0;'])

    eval(['separationSis12' label{iLabel,1} '(1:numSistersTot,1) = struct(''observations'',[],''time'',[]);'])
    eval(['separationChangeSis12' label{iLabel,1} '(1:numSistersTot,1) = struct(''observations'',[],''time'',[]);'])
    eval(['angleNormSis12' label{iLabel,1} '(1:numSistersTot,1) = struct(''observations'',[],''time'',[]);'])
    eval(['angularDispSis12' label{iLabel,1} '(1:numSistersTot,1) = struct(''observations'',[],''time'',[]);'])

    eval(['centerPositionSis12' label{iLabel,1} '(1:numSistersTot,1) = struct(''observations'',[],''time'',[]);'])
    eval(['centerPositionChangeSis12' label{iLabel,1} '(1:numSistersTot,1) = struct(''observations'',[],''time'',[]);'])

    eval(['dispMagSis1' label{iLabel,1} '(1:numSistersTot,1) = struct(''observations'',[]);'])
    eval(['dispMagSis2' label{iLabel,1} '(1:numSistersTot,1) = struct(''observations'',[]);'])
    eval(['dispProjSis1' label{iLabel,1} '(1:numSistersTot,1) = struct(''observations'',[]);'])
    eval(['dispProjSis2' label{iLabel,1} '(1:numSistersTot,1) = struct(''observations'',[]);'])
    eval(['dispAngleSis1' label{iLabel,1} '(1:numSistersTot,1) = struct(''observations'',[]);'])
    eval(['dispAngleSis2' label{iLabel,1} '(1:numSistersTot,1) = struct(''observations'',[]);'])

    eval(['dotProdCrossProdSign12' label{iLabel,1} '(1:numSistersTot,1) = struct(''observations'',[]);'])
    
end

%go over all movies
for iMovie = 1 : numMovies

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

            %% sort sisters to left and right
            
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

            %% sister separation & its change

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
            
            %% sister angle with normal and angular velocity
            
            %normalize vector connecting sisters
            sisterVec = sisterVec ./ repmat(sisterDist(:,1),1,3); %unitless
            
            %calculate angle with normal
            angleWithNorm = NaN(numFrames,1);
            for iFrame = framesWithPlane
                angleWithNorm(iFrame) = acos(abs(sisterVec(iFrame,1))) * 180 / pi; %degrees
            end

            %calculate angle between consecutive frames
            angularDisp = NaN(numFrames-1,1);
            for iFrame = 1 : numFrames - 1
                angularDisp(iFrame) = acos(abs(sisterVec(iFrame,:)*sisterVec(iFrame+1,:)')) * 180 / pi; %degrees
            end

            %% center position & its change
            
            %calculate position of sister center of mass along the normal
            %to the metaphase plate
            centerPos = mean([coords1(:,1) coords2(:,1)],2); %um
            centerPosStd = 0.5 * sqrt( sum( [coords1(:,4) coords2(:,4)].^2 ,2) ); %um

            %calculate change in position of center of mass between frames
            centerPosChange = diff(centerPos); %um
            centerPosChangeStd = sqrt( sum( [centerPosStd(2:end) ...
                centerPosStd(1:end-1)].^2 ,2) ); %um
            
            %%sister displacement between frames
            
            %calculate sister displacement between frames
            coords1Diff = [coords1(2:end,1:3)-coords1(1:end-1,1:3) ...
                sqrt(coords1(2:end,4:6).^2+coords1(1:end-1,4:6).^2)]; %um
            coords2Diff = [coords2(2:end,1:3)-coords2(1:end-1,1:3) ...
                sqrt(coords2(2:end,4:6).^2+coords2(1:end-1,4:6).^2)]; %um

            %calculate projection of displacements on vector connecting sisters
            projection1 = diag( coords1Diff(:,1:3) * sisterVec(1:end-1,1:3)' ); %um
            projection2 = diag( coords2Diff(:,1:3) * sisterVec(1:end-1,1:3)' ); %um

            %calculate displacement magnitude
            %add minus sign to those displacements whose projection is
            %negative
            dispMag1 = sign(projection1) .* sqrt(sum(coords1Diff(:,1:3).^2,2)); %um
            dispMag2 = sign(projection2) .* sqrt(sum(coords2Diff(:,1:3).^2,2)); %um
            
            %calculate angle between displacements and vector connecting sisters
            angle1 = acos(projection1 ./ sqrt(sum(coords1Diff(:,1:3).^2,2))) * 180 / pi; %degrees
            angle2 = acos(projection2 ./ sqrt(sum(coords2Diff(:,1:3).^2,2))) * 180 / pi; %degrees

            %use the cross product to deduce whether the sisters
            %"displaced" or "tumbled"
            crossProd1 = cross(coords1Diff(:,1:3),sisterVec(1:end-1,1:3),2);
            crossProd2 = cross(coords2Diff(:,1:3),sisterVec(1:end-1,1:3),2);
            cosAngleCrossProd12 = diag(crossProd1 * crossProd2');
            crossProdSign = sign(cosAngleCrossProd12); %unitless
            
            %% time

            %store time vector
            timeVec = [(0:numFrames-1)'*timeLapse zeros(numFrames,1)]; %s
            
            %store sister information
            iLabel = sisterType(iSister) + 1;

            %increase global index of sister type by 1
            eval(['iGlobal' label{iLabel,1} ' = iGlobal' label{iLabel,1} ' + 1;'])

            %store information
            eval(['separationSis12' label{iLabel,1} '(iGlobal' label{iLabel,1} ').observations = [sisterDist sisterDistStd];']) %um
            eval(['separationChangeSis12' label{iLabel,1} '(iGlobal' label{iLabel,1} ').observations = [sisterDistChange sisterDistChangeStd];']) %um
            eval(['angleNormSis12' label{iLabel,1} '(iGlobal' label{iLabel,1} ').observations = angleWithNorm;']) %um
            eval(['angularDispSis12' label{iLabel,1} '(iGlobal' label{iLabel,1} ').observations = angularDisp;']) %um
            eval(['centerPositionSis12' label{iLabel,1} '(iGlobal' label{iLabel,1} ').observations = [centerPos centerPosStd];']) %um
            eval(['centerPositionChangeSis12' label{iLabel,1} '(iGlobal' label{iLabel,1} ').observations = [centerPosChange centerPosChangeStd];']) %um
            eval(['dispMagSis1' label{iLabel,1} '(iGlobal' label{iLabel,1} ').observations = dispMag1;'])
            eval(['dispMagSis2' label{iLabel,1} '(iGlobal' label{iLabel,1} ').observations = dispMag2;'])
            eval(['dispProjSis1' label{iLabel,1} '(iGlobal' label{iLabel,1} ').observations = projection1;'])
            eval(['dispProjSis2' label{iLabel,1} '(iGlobal' label{iLabel,1} ').observations = projection2;'])
            eval(['dispAngleSis1' label{iLabel,1} '(iGlobal' label{iLabel,1} ').observations = angle1;'])
            eval(['dispAngleSis2' label{iLabel,1} '(iGlobal' label{iLabel,1} ').observations = angle2;'])            
            eval(['dotProdCrossProdSign12' label{iLabel,1} '(iGlobal' label{iLabel,1} ').observations = crossProdSign;'])            
            
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

    eval(['separationSis12' label{iLabel,1} ' = separationSis12' label{iLabel,1} '(1:iGlobal' label{iLabel,1} ');'])
    eval(['separationChangeSis12' label{iLabel,1} ' = separationChangeSis12' label{iLabel,1} '(1:iGlobal' label{iLabel,1} ');'])
    eval(['angleNormSis12' label{iLabel,1} ' = angleNormSis12' label{iLabel,1} '(1:iGlobal' label{iLabel,1} ');'])
    eval(['angularDispSis12' label{iLabel,1} ' = angularDispSis12' label{iLabel,1} '(1:iGlobal' label{iLabel,1} ');'])
    eval(['centerPositionSis12' label{iLabel,1} ' = centerPositionSis12' label{iLabel,1} '(1:iGlobal' label{iLabel,1} ');'])
    eval(['centerPositionChangeSis12' label{iLabel,1} ' = centerPositionChangeSis12' label{iLabel,1} '(1:iGlobal' label{iLabel,1} ');'])
    eval(['dispMagSis1' label{iLabel,1} ' = dispMagSis1' label{iLabel,1} '(1:iGlobal' label{iLabel,1} ');'])
    eval(['dispMagSis2' label{iLabel,1} ' = dispMagSis2' label{iLabel,1} '(1:iGlobal' label{iLabel,1} ');'])
    eval(['dispProjSis1' label{iLabel,1} ' = dispProjSis1' label{iLabel,1} '(1:iGlobal' label{iLabel,1} ');'])
    eval(['dispProjSis2' label{iLabel,1} ' = dispProjSis2' label{iLabel,1} '(1:iGlobal' label{iLabel,1} ');'])
    eval(['dispAngleSis1' label{iLabel,1} ' = dispAngleSis1' label{iLabel,1} '(1:iGlobal' label{iLabel,1} ');'])
    eval(['dispAngleSis2' label{iLabel,1} ' = dispAngleSis2' label{iLabel,1} '(1:iGlobal' label{iLabel,1} ');'])
    eval(['dotProdCrossProdSign12' label{iLabel,1} ' = dotProdCrossProdSign12' label{iLabel,1} '(1:iGlobal' label{iLabel,1} ');'])
    
end
    
%% sister oscillation state analysis

%initialization
for iLabel = 1 : 3
    eval(['pointState' label{iLabel,1} ' = repmat(struct(''observations'',[]),iGlobal' ...
        label{iLabel,1} ',1);'])
end

%calculation
for iLabel = goodLabel

    %go over all center of mass trajectories
    for iSister = 1 : eval(['iGlobal' label{iLabel,1}])
        
        %get the sister pair's center of mass positions
        eval(['comPosCurrent = centerPositionSis12' label{iLabel} ...
            '(iSister).observations;'])
        
        %get the sister pair's center of mass position changes
        eval(['comPosChangeCurrent = centerPositionChangeSis12' label{iLabel} ...
            '(iSister).observations;'])
        
        %calculate the p-values of these changes, assuming that they follow
        %a normal disitrubution with mean 0 and std given by the second
        %column in comPosChangeCurrent
        pValueChange = normcdf(comPosChangeCurrent(:,1),0,...
            comPosChangeCurrent(:,2));
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
        trajStates = NaN(size(comPosChangeCurrent,1),1);
        trajStates(comPosChangeCurrent(:,1) >= 0 & pValueChange < alpha/2) =  1;
        trajStates(comPosChangeCurrent(:,1) >= 0 & pValueChange > alpha/2) = -1;
        trajStates(comPosChangeCurrent(:,1) <  0 & pValueChange < alpha/2) =  2;
        trajStates(comPosChangeCurrent(:,1) <  0 & pValueChange > alpha/2) = -2;
                
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
                trajValues = comPosCurrent(goodIntervals(iInt)+1:...
                    goodIntervals(iInt+1),1);
                
                %find the point with maximum value
                maxPoint = find(trajValues==max(trajValues));
                
                %give all intervals before this point a state of +1
                trajStates(goodIntervals(iInt)+1:goodIntervals(iInt)+maxPoint-1) = 1;
                
                %give all intervals after his point a state of +2;
                trajStates(goodIntervals(iInt)+maxPoint:goodIntervals(iInt+1)-1) = 2;
                
            else %transition point from neg. to pos. movement
                
                %get position values in this undetermined interval
                trajValues = comPosCurrent(goodIntervals(iInt)+1:...
                    goodIntervals(iInt+1),1);
                
                %find the point with minimum value
                minPoint = find(trajValues==min(trajValues));
                
                %give all intervals before this point a state of +2
                trajStates(goodIntervals(iInt)+1:goodIntervals(iInt)+minPoint-1) = 2;
                
                %give all intervals after this point a state of +1
                trajStates(goodIntervals(iInt)+minPoint:goodIntervals(iInt+1)-1) = 1;
                
            end
            
        end
        
        %Step 1: find transition points
        
        %take the first difference of trajStates to find transition points
        %+ve differences indicate transition from growth to shrinkage
        %-ve differences indicate transition from shrinkage to growth
        trajStatesDiff = diff(trajStates);
        pos2negTransPoints = find( trajStatesDiff > 0 ) + 1;
        neg2posTransPoints = find( trajStatesDiff < 0 ) + 1;
        
        %Step 2: find the points at the center of a state

        %e.g., if two transition points are at 4 and 10, the center
        %will be at 7; if they are at 4 and 11, there will be two center
        %points at 7 and 8; center points will only be chosen for a state
        %lasting at least 5 points, without missing points in between
        
        %find missing data points
        missDataPoints = find(isnan(comPosCurrent(:,1)));
        
        %concatenate the two transition point vectors and the vector of
        %missing data points
        transPointsVector = [pos2negTransPoints; neg2posTransPoints; ...
            missDataPoints];
        transPointsType = [ones(size(pos2negTransPoints)); ...
            -ones(size(neg2posTransPoints)); NaN(size(missDataPoints))];
        
        %sort transPointsVector in ascending order
        %sort transPointsType in the same way
        [transPointsVector,indxSort] = sort(transPointsVector);
        transPointsType = transPointsType(indxSort);
        
        %subtract consecutive elements in transPointsType to decide on
        %types of movement states
        % 2  = movement in the positive direction
        %-2  = movement in the negative direction
        %NaN = a transition point neighboring a missing data point
        movementType = diff(transPointsType);
        
        %subtract consecutive elements in transPointsVector to determine
        %the duration of each movement state
        movementDuration = diff(transPointsVector) + 1;
        
        %find movement states lasting for at least 5 consecutive points, so
        %that we find their center points
        goodIntervPos = find(movementType ==  2 & movementDuration >= 5);
        goodIntervNeg = find(movementType == -2 & movementDuration >= 5);
        
        %get the center points
        posCenterPoints = (transPointsVector(goodIntervPos) + ...
            transPointsVector(goodIntervPos+1))/2;
        posCenterPoints = sort(unique([floor(posCenterPoints); ...
            ceil(posCenterPoints)]));
        negCenterPoints = (transPointsVector(goodIntervNeg) + ...
            transPointsVector(goodIntervNeg+1))/2;
        negCenterPoints = sort(unique([floor(negCenterPoints); ...
            ceil(negCenterPoints)]));
        
        %Step 3: find general points

        %these are non-center points in intervals bounded by transition
        %points and without NaNs in between
        %they have to be one point away from any transition point and
        %center point
        %thay come from intervals that are at least 9 points long
        
        %find movement states lasting for at least 9 consecutive points, so
        %that we find their general points
        goodIntervPos = find(movementType ==  2 & movementDuration >= 9);
        goodIntervNeg = find(movementType == -2 & movementDuration >= 9);
        
        %extract all time points in positive movement intervals that are
        %not transition points or next to them
        posGenPoints = [];
        for iInt = goodIntervPos'
            posGenPoints = [posGenPoints; (transPointsVector(iInt)...
                +2:transPointsVector(iInt+1)-2)'];
        end
        posGenPoints = setdiff(posGenPoints,unique([posCenterPoints;...
            posCenterPoints-1;posCenterPoints+1]));
               
        %extract all time points in negative movement intervals that are
        %not transition points or center points
        negGenPoints = [];
        for iInt = goodIntervNeg'
            negGenPoints = [negGenPoints; (transPointsVector(iInt)...
                +2:transPointsVector(iInt+1)-2)'];
        end
        negGenPoints = setdiff(negGenPoints,unique([negCenterPoints;...
            negCenterPoints-1;negCenterPoints+1]));
               
        %store the state of each time point in the variable "pointState"
        %+1: point of transition from positive to negative movement
        %-1: point of transition from negative to positive movement
        %+2: center point in a positive movement interval
        %-2: center point in a negative movement interval
        %+3: non-center, non-transition point in a positive movement interval
        %-3: non-center, non-transition point in a negative movement interval
        % 0: point that does not fall in any of the above categories
        eval(['numTimePoints = size(separationSis12' label{iLabel,1} '(iSister).observations,1);'])
        eval(['pointState' label{iLabel,1} '(iSister).observations = zeros(numTimePoints,1);'])
        eval(['pointState' label{iLabel,1} '(iSister).observations(pos2negTransPoints) = 1;'])
        eval(['pointState' label{iLabel,1} '(iSister).observations(neg2posTransPoints) = -1;'])
        eval(['pointState' label{iLabel,1} '(iSister).observations(posCenterPoints) = 2;'])
        eval(['pointState' label{iLabel,1} '(iSister).observations(negCenterPoints) = -2;'])
        eval(['pointState' label{iLabel,1} '(iSister).observations(posGenPoints) = 3;'])
        eval(['pointState' label{iLabel,1} '(iSister).observations(negGenPoints) = -3;'])
                
    end

end

%% overall distributions and their parameters

%initialization
for iLabel = 1 : 3

    eval(['separationBeforeDistr' label{iLabel,1} ' = [];'])
    eval(['separationAfterDistr' label{iLabel,1} ' = [];'])
    eval(['separationChangeDistr' label{iLabel,1} ' = [];'])
    eval(['separationBeforeParam' label{iLabel,1} ' = [];'])
    eval(['separationAfterParam' label{iLabel,1} ' = [];'])
    eval(['separationChangeParam' label{iLabel,1} ' = [];'])
    
    eval(['angleNormBeforeDistr' label{iLabel,1} ' = [];'])
    eval(['angleNormAfterDistr' label{iLabel,1} ' = [];'])
    eval(['angularDispDistr' label{iLabel,1} ' = [];'])
    eval(['angleNormBeforeParam' label{iLabel,1} ' = [];'])
    eval(['angleNormAfterParam' label{iLabel,1} ' = [];'])
    eval(['angularDispParam' label{iLabel,1} ' = [];'])

    eval(['centerPositionBeforeDistr' label{iLabel,1} ' = [];'])
    eval(['centerPositionAfterDistr' label{iLabel,1} ' = [];'])
    eval(['centerPositionChangeDistr' label{iLabel,1} ' = [];'])
    eval(['centerPositionBeforeParam' label{iLabel,1} ' = [];'])
    eval(['centerPositionAfterParam' label{iLabel,1} ' = [];'])
    eval(['centerPositionChangeParam' label{iLabel,1} ' = [];'])

    eval(['displacement1Distr' label{iLabel,1} ' = [];'])
    eval(['dispProj1Distr' label{iLabel,1} ' = [];'])
    eval(['dispAngle1Distr' label{iLabel,1} ' = [];'])
    eval(['displacement1Param' label{iLabel,1} ' = [];'])
    eval(['dispProj1Param' label{iLabel,1} ' = [];'])
    eval(['dispAngle1Param' label{iLabel,1} ' = [];'])
    
    eval(['displacement2Distr' label{iLabel,1} ' = [];'])
    eval(['dispProj2Distr' label{iLabel,1} ' = [];'])
    eval(['dispAngle2Distr' label{iLabel,1} ' = [];'])
    eval(['displacement2Param' label{iLabel,1} ' = [];'])
    eval(['dispProj2Param' label{iLabel,1} ' = [];'])
    eval(['dispAngle2Param' label{iLabel,1} ' = [];'])
    
    eval(['dotProdCrossProdSignDistr' label{iLabel,1} ' = [];'])
    
    eval(['pointStateBeforeDistr' label{iLabel,1} ' = [];'])
    eval(['pointStateAfterDistr' label{iLabel,1} ' = [];'])    
    
end

for iLabel = goodLabel

    %sister separation change
    eval(['allValues = vertcat(separationChangeSis12' label{iLabel,1} '.observations);']);
    allValues = allValues(:,1);
    indxCommon = find(~isnan(allValues));
    allValues = allValues(indxCommon);
    eval(['separationChangeDistr' label{iLabel,1} ' = allValues;']);
    distrParam = [nanmean(allValues) nanstd(allValues) min(allValues) ...
        prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) ...
        max(allValues); ...
        nanmean(allValues(allValues>0)) nanstd(allValues(allValues>0)) ...
        min(allValues(allValues>0)) prctile(allValues(allValues>0),25) ...
        prctile(allValues(allValues>0),50) prctile(allValues(allValues>0),75) ...
        max(allValues(allValues>0)); ...
        nanmean(allValues(allValues<0)) nanstd(allValues(allValues<0)) ...
        min(allValues(allValues<0)) prctile(allValues(allValues<0),25) ...
        prctile(allValues(allValues<0),50) prctile(allValues(allValues<0),75) ...
        max(allValues(allValues<0))];
    eval(['separationChangeParam' label{iLabel,1} ' = distrParam;'])

    %sister separation, before and after
    eval(['allValues = separationSis12' label{iLabel,1} ';'])
    valuesBef = [];
    valuesAft = [];
    for i = 1 : length(allValues)
        valuesBef = [valuesBef; allValues(i).observations(1:end-1,1)];
        valuesAft = [valuesAft; allValues(i).observations(2:end,1)];
    end
    valuesBef = valuesBef(indxCommon);
    valuesAft = valuesAft(indxCommon);
    eval(['separationBeforeDistr' label{iLabel,1} ' = valuesBef;']);
    eval(['separationAfterDistr' label{iLabel,1} ' = valuesAft;']);
    distrParam = [nanmean(valuesBef) nanstd(valuesBef) min(valuesBef) ...
        prctile(valuesBef,25) prctile(valuesBef,50) prctile(valuesBef,75) ...
        max(valuesBef)];
    eval(['separationBeforeParam' label{iLabel,1} ' = distrParam;']);
    distrParam = [nanmean(valuesAft) nanstd(valuesAft) min(valuesAft) ...
        prctile(valuesAft,25) prctile(valuesAft,50) prctile(valuesAft,75) ...
        max(valuesAft)];
    eval(['separationAfterParam' label{iLabel,1} ' = distrParam;']);

    %angle with normal, before and after
    eval(['allValues = angleNormSis12' label{iLabel,1} ';'])
    valuesBef = [];
    valuesAft = [];
    for i = 1 : length(allValues)
        valuesBef = [valuesBef; allValues(i).observations(1:end-1,1)];
        valuesAft = [valuesAft; allValues(i).observations(2:end,1)];
    end
    valuesBef = valuesBef(indxCommon);
    valuesAft = valuesAft(indxCommon);
    eval(['angleNormBeforeDistr' label{iLabel,1} ' = valuesBef;']);
    eval(['angleNormAfterDistr' label{iLabel,1} ' = valuesAft;']);
    distrParam = [nanmean(valuesBef) nanstd(valuesBef) min(valuesBef) ...
        prctile(valuesBef,25) prctile(valuesBef,50) prctile(valuesBef,75) ...
        max(valuesBef)];
    eval(['angleNormBeforeParam' label{iLabel,1} ' = distrParam;']);
    distrParam = [nanmean(valuesAft) nanstd(valuesAft) min(valuesAft) ...
        prctile(valuesAft,25) prctile(valuesAft,50) prctile(valuesAft,75) ...
        max(valuesAft)];
    eval(['angleNormAfterParam' label{iLabel,1} ' = distrParam;']);

    %angular displacement
    eval(['allValues = vertcat(angularDispSis12' label{iLabel,1} '.observations);']);
    allValues = allValues(indxCommon,1);
    eval(['angularDispDistr' label{iLabel,1} ' = allValues;']);
    distrParam = [nanmean(allValues) nanstd(allValues) min(allValues) ...
        prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) ...
        max(allValues)];
    eval(['angularDispParam' label{iLabel,1} ' = distrParam;'])

    %center position, before and after
    eval(['allValues = centerPositionSis12' label{iLabel,1} ';'])
    valuesBef = [];
    valuesAft = [];
    for i = 1 : length(allValues)
        valuesBef = [valuesBef; allValues(i).observations(1:end-1,1)];
        valuesAft = [valuesAft; allValues(i).observations(2:end,1)];
    end
    valuesBef = valuesBef(indxCommon);
    valuesAft = valuesAft(indxCommon);
    eval(['centerPositionBeforeDistr' label{iLabel,1} ' = valuesBef;']);
    eval(['centerPositionAfterDistr' label{iLabel,1} ' = valuesAft;']);
    distrParam = [nanmean(valuesBef) nanstd(valuesBef) min(valuesBef) ...
        prctile(valuesBef,25) prctile(valuesBef,50) prctile(valuesBef,75) ...
        max(valuesBef)];
    eval(['centerPositionBeforeParam' label{iLabel,1} ' = distrParam;']);
    distrParam = [nanmean(valuesAft) nanstd(valuesAft) min(valuesAft) ...
        prctile(valuesAft,25) prctile(valuesAft,50) prctile(valuesAft,75) ...
        max(valuesAft)];
    eval(['centerPositionAfterParam' label{iLabel,1} ' = distrParam;']);

    %center position change
    eval(['allValues = vertcat(centerPositionChangeSis12' label{iLabel,1} '.observations);']);
    allValues = allValues(indxCommon,1);
    eval(['centerPositionChangeDistr' label{iLabel,1} ' = allValues;']);
    distrParam = [nanmean(allValues) nanstd(allValues) min(allValues) ...
        prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) ...
        max(allValues); ...
        nanmean(allValues(allValues>0)) nanstd(allValues(allValues>0)) ...
        min(allValues(allValues>0)) prctile(allValues(allValues>0),25) ...
        prctile(allValues(allValues>0),50) prctile(allValues(allValues>0),75) ...
        max(allValues(allValues>0)); ...
        nanmean(allValues(allValues<0)) nanstd(allValues(allValues<0)) ...
        min(allValues(allValues<0)) prctile(allValues(allValues<0),25) ...
        prctile(allValues(allValues<0),50) prctile(allValues(allValues<0),75) ...
        max(allValues(allValues<0))];
    eval(['centerPositionChangeParam' label{iLabel,1} ' = distrParam;'])

    %sister 1 displacement magnitude
    eval(['allValues = vertcat(dispMagSis1' label{iLabel,1} '.observations);']);
    allValues = allValues(indxCommon,1);
    eval(['displacement1Distr' label{iLabel,1} ' = allValues;']);
    distrParam = [nanmean(allValues) nanstd(allValues) min(allValues) ...
        prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) ...
        max(allValues); ...
        nanmean(allValues(allValues>0)) nanstd(allValues(allValues>0)) ...
        min(allValues(allValues>0)) prctile(allValues(allValues>0),25) ...
        prctile(allValues(allValues>0),50) prctile(allValues(allValues>0),75) ...
        max(allValues(allValues>0)); ...
        nanmean(allValues(allValues<0)) nanstd(allValues(allValues<0)) ...
        min(allValues(allValues<0)) prctile(allValues(allValues<0),25) ...
        prctile(allValues(allValues<0),50) prctile(allValues(allValues<0),75) ...
        max(allValues(allValues<0))];
    eval(['displacement1Param' label{iLabel,1} ' = distrParam;'])

    %sister 1 displacement projection
    eval(['allValues = vertcat(dispProjSis1' label{iLabel,1} '.observations);']);
    allValues = allValues(indxCommon,1);
    eval(['dispProj1Distr' label{iLabel,1} ' = allValues;']);
    distrParam = [nanmean(allValues) nanstd(allValues) min(allValues) ...
        prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) ...
        max(allValues); ...
        nanmean(allValues(allValues>0)) nanstd(allValues(allValues>0)) ...
        min(allValues(allValues>0)) prctile(allValues(allValues>0),25) ...
        prctile(allValues(allValues>0),50) prctile(allValues(allValues>0),75) ...
        max(allValues(allValues>0)); ...
        nanmean(allValues(allValues<0)) nanstd(allValues(allValues<0)) ...
        min(allValues(allValues<0)) prctile(allValues(allValues<0),25) ...
        prctile(allValues(allValues<0),50) prctile(allValues(allValues<0),75) ...
        max(allValues(allValues<0))];
    eval(['dispProj1Param' label{iLabel,1} ' = distrParam;'])

    %sister 1 displacement angle
    eval(['allValues = vertcat(dispAngleSis1' label{iLabel,1} '.observations);']);
    allValues = allValues(indxCommon,1);
    eval(['dispAngle1Distr' label{iLabel,1} ' = allValues;']);
    distrParam = [nanmean(allValues) nanstd(allValues) min(allValues) ...
        prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) ...
        max(allValues); ...
        nanmean(allValues(allValues>90)) nanstd(allValues(allValues>90)) ...
        min(allValues(allValues>90)) prctile(allValues(allValues>90),25) ...
        prctile(allValues(allValues>90),50) prctile(allValues(allValues>90),75) ...
        max(allValues(allValues>90)); ...
        nanmean(allValues(allValues<90)) nanstd(allValues(allValues<90)) ...
        min(allValues(allValues<90)) prctile(allValues(allValues<90),25) ...
        prctile(allValues(allValues<90),50) prctile(allValues(allValues<90),75) ...
        max(allValues(allValues<90))];
    eval(['dispAngle1Param' label{iLabel,1} ' = distrParam;'])

    %sister 2 displacement magnitude
    eval(['allValues = vertcat(dispMagSis2' label{iLabel,1} '.observations);']);
    allValues = allValues(indxCommon,1);
    eval(['displacement2Distr' label{iLabel,1} ' = allValues;']);
    distrParam = [nanmean(allValues) nanstd(allValues) min(allValues) ...
        prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) ...
        max(allValues); ...
        nanmean(allValues(allValues>0)) nanstd(allValues(allValues>0)) ...
        min(allValues(allValues>0)) prctile(allValues(allValues>0),25) ...
        prctile(allValues(allValues>0),50) prctile(allValues(allValues>0),75) ...
        max(allValues(allValues>0)); ...
        nanmean(allValues(allValues<0)) nanstd(allValues(allValues<0)) ...
        min(allValues(allValues<0)) prctile(allValues(allValues<0),25) ...
        prctile(allValues(allValues<0),50) prctile(allValues(allValues<0),75) ...
        max(allValues(allValues<0))];
    eval(['displacement2Param' label{iLabel,1} ' = distrParam;'])

    %sister 2 displacement projection
    eval(['allValues = vertcat(dispProjSis2' label{iLabel,1} '.observations);']);
    allValues = allValues(indxCommon,1);
    eval(['dispProj2Distr' label{iLabel,1} ' = allValues;']);
    distrParam = [nanmean(allValues) nanstd(allValues) min(allValues) ...
        prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) ...
        max(allValues); ...
        nanmean(allValues(allValues>0)) nanstd(allValues(allValues>0)) ...
        min(allValues(allValues>0)) prctile(allValues(allValues>0),25) ...
        prctile(allValues(allValues>0),50) prctile(allValues(allValues>0),75) ...
        max(allValues(allValues>0)); ...
        nanmean(allValues(allValues<0)) nanstd(allValues(allValues<0)) ...
        min(allValues(allValues<0)) prctile(allValues(allValues<0),25) ...
        prctile(allValues(allValues<0),50) prctile(allValues(allValues<0),75) ...
        max(allValues(allValues<0))];
    eval(['dispProj2Param' label{iLabel,1} ' = distrParam;'])

    %sister 2 displacement angle
    eval(['allValues = vertcat(dispAngleSis2' label{iLabel,1} '.observations);']);
    allValues = allValues(indxCommon,1);
    eval(['dispAngle2Distr' label{iLabel,1} ' = allValues;']);
    distrParam = [nanmean(allValues) nanstd(allValues) min(allValues) ...
        prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) ...
        max(allValues); ...
        nanmean(allValues(allValues>90)) nanstd(allValues(allValues>90)) ...
        min(allValues(allValues>90)) prctile(allValues(allValues>90),25) ...
        prctile(allValues(allValues>90),50) prctile(allValues(allValues>90),75) ...
        max(allValues(allValues>90)); ...
        nanmean(allValues(allValues<90)) nanstd(allValues(allValues<90)) ...
        min(allValues(allValues<90)) prctile(allValues(allValues<90),25) ...
        prctile(allValues(allValues<90),50) prctile(allValues(allValues<90),75) ...
        max(allValues(allValues<90))];
    eval(['dispAngle2Param' label{iLabel,1} ' = distrParam;'])

    %sign of dot product of cross products
    eval(['allValues = vertcat(dotProdCrossProdSign12' label{iLabel,1} '.observations);']);
    allValues = allValues(indxCommon,1);
    eval(['dotProdCrossProdSignDistr' label{iLabel,1} ' = allValues;']);

    %point state, before and after
    eval(['allValues = pointState' label{iLabel,1} ';'])
    valuesBef = [];
    valuesAft = [];
    for i = 1 : length(allValues)
        valuesBef = [valuesBef; allValues(i).observations(1:end-1,1)];
        valuesAft = [valuesAft; allValues(i).observations(2:end,1)];
    end
    valuesBef = valuesBef(indxCommon);
    valuesAft = valuesAft(indxCommon);
    eval(['pointStateBeforeDistr' label{iLabel,1} ' = valuesBef;']);
    eval(['pointStateAfterDistr' label{iLabel,1} ' = valuesAft;']);
    
end

%% indices of instances in various categories

%initialization
for iLabel = 1 : 3
    eval(['indx1R2R' label{iLabel,1} ' = [];'])
    eval(['indx1R2L' label{iLabel,1} ' = [];'])
    eval(['indx1L2R' label{iLabel,1} ' = [];'])
    eval(['indx1L2L' label{iLabel,1} ' = [];'])
    eval(['indxTransBef' label{iLabel,1} ' = [];'])
    eval(['indxTransAft' label{iLabel,1} ' = [];'])
    eval(['indxCenterBef' label{iLabel,1} ' = [];'])
    eval(['indxCenterAft' label{iLabel,1} ' = [];'])
    eval(['indxGeneralBef' label{iLabel,1} ' = [];'])
    eval(['indxGeneralAft' label{iLabel,1} ' = [];'])
    eval(['indxTumble' label{iLabel,1} ' = [];'])
    eval(['indxNoTumble' label{iLabel,1} ' = [];'])
    eval(['indxAll' label{iLabel,1} ' = [];'])
    eval(['indxClassBef' label{iLabel,1} ' = [];'])
    eval(['indxClassAft' label{iLabel,1} ' = [];'])
end

%calculation
for iLabel = goodLabel

    %read out variables for this category of sisters
    eval(['dispProj1Distr = dispProj1Distr' label{iLabel,1} ';'])
    eval(['dispProj2Distr = dispProj2Distr' label{iLabel,1} ';'])
    eval(['pointStateBeforeDistr = pointStateBeforeDistr' label{iLabel,1} ';'])
    eval(['pointStateAfterDistr = pointStateAfterDistr' label{iLabel,1} ';'])
    eval(['dotProdCrossProdSignDistr = dotProdCrossProdSignDistr' label{iLabel,1} ';'])

    %get intervals where sisters move ...
    %both to the right
    eval(['indx1R2R' label{iLabel,1} ' = find(dispProj1Distr > 0 & dispProj2Distr > 0);'])
    %the first to the right and the second to the left
    eval(['indx1R2L' label{iLabel,1} ' = find(dispProj1Distr > 0 & dispProj2Distr < 0);'])
    %the first to the left and the second to the right
    eval(['indx1L2R' label{iLabel,1} ' = find(dispProj1Distr < 0 & dispProj2Distr > 0);'])
    %both to the left
    eval(['indx1L2L' label{iLabel,1} ' = find(dispProj1Distr < 0 & dispProj2Distr < 0);'])
    
    %put all of these indices together
    eval(['indxAll' label{iLabel,1} ' = sort([indx1R2R' label{iLabel,1} ...
        '; indx1R2L' label{iLabel,1} '; indx1L2R' label{iLabel,1} ...
        '; indx1L2L' label{iLabel,1} ']);'])

    %get instances where points are ...
    %transition points
    eval(['indxTransBef' label{iLabel,1} ' = find(abs(pointStateBeforeDistr)==1);'])
    eval(['indxTransAft' label{iLabel,1} ' = find(abs(pointStateAfterDistr)==1);'])
    %center points
    eval(['indxCenterBef' label{iLabel,1} ' = find(abs(pointStateBeforeDistr)==2);'])
    eval(['indxCenterAft' label{iLabel,1} ' = find(abs(pointStateAfterDistr)==2);'])
    %general points
    eval(['indxGeneralBef' label{iLabel,1} ' = find(abs(pointStateBeforeDistr)==3);'])
    eval(['indxGeneralAft' label{iLabel,1} ' = find(abs(pointStateAfterDistr)==3);'])
    
    %put these indices together to get all instances where point state is
    %classified
    eval(['indxClassBef' label{iLabel,1} ' = sort(unique([indxTransBef' ...
        label{iLabel,1} '; indxCenterBef' label{iLabel,1} ...
        '; indxGeneralBef' label{iLabel,1} ']));'])
    eval(['indxClassAft' label{iLabel,1} ' = sort(unique([indxTransAft' ...
        label{iLabel,1} '; indxCenterAft' label{iLabel,1} ...
        '; indxGeneralAft' label{iLabel,1} ']));'])
    
    %get intervals where sisters tumble and do not tumble
    eval(['indxTumble' label{iLabel,1} ' = find(dotProdCrossProdSignDistr==-1);'])
    eval(['indxNoTumble' label{iLabel,1} ' = find(dotProdCrossProdSignDistr==1);'])

end

%% conditional probabilities

%initialization
for iLabel = 1 : 3
    eval(['probBothRight' label{iLabel,1} ' = [];'])
    eval(['probBothLeft' label{iLabel,1} ' = [];'])
    eval(['probApproach' label{iLabel,1} ' = [];'])
    eval(['probSeparate' label{iLabel,1} ' = [];'])
    eval(['probTransBef' label{iLabel,1} ' = [];'])
    eval(['probCenterBef' label{iLabel,1} ' = [];'])
    eval(['probGeneralBef' label{iLabel,1} ' = [];'])
    eval(['probTransAft' label{iLabel,1} ' = [];'])
    eval(['probCenterAft' label{iLabel,1} ' = [];'])
    eval(['probGeneralAft' label{iLabel,1} ' = [];'])
    eval(['probTumble' label{iLabel,1} ' = [];'])
end

%calculation
for iLabel = goodLabel
   
    %get indices
    eval(['indx1R2R = indx1R2R' label{iLabel,1} ';'])
    eval(['indx1R2L = indx1R2L' label{iLabel,1} ';'])
    eval(['indx1L2R = indx1L2R' label{iLabel,1} ';'])
    eval(['indx1L2L = indx1L2L' label{iLabel,1} ';'])
    eval(['indxTransBef = indxTransBef' label{iLabel,1} ';'])
    eval(['indxTransAft = indxTransAft' label{iLabel,1} ';'])
    eval(['indxCenterBef = indxCenterBef' label{iLabel,1} ';'])
    eval(['indxCenterAft = indxCenterAft' label{iLabel,1} ';'])
    eval(['indxGeneralBef = indxGeneralBef' label{iLabel,1} ';'])
    eval(['indxGeneralAft = indxGeneralAft' label{iLabel,1} ';'])
    eval(['indxTumble = indxTumble' label{iLabel,1} ';'])

    %get number of instances in each category
    numBothRight = length(indx1R2R);
    numBothLeft = length(indx1L2L);
    numApproach = length(indx1R2L);
    numSeparate = length(indx1L2R);
    numTransBef = length(indxTransBef);
    numTransAft = length(indxTransAft);
    numCenterBef = length(indxCenterBef);
    numCenterAft = length(indxCenterAft);
    numGeneralBef = length(indxGeneralBef);
    numGeneralAft = length(indxGeneralAft);
    numTumble = length(indxTumble);
    
    %get total number of instances
    numTotInstances = sum([numBothRight numBothLeft numApproach numSeparate]);
    
    %get number of instances of joint point states AND sister motion
    %state before
    numBothRightAndTransBef   = length(intersect(indx1R2R,indxTransBef));
    numBothRightAndCenterBef  = length(intersect(indx1R2R,indxCenterBef));
    numBothRightAndGeneralBef = length(intersect(indx1R2R,indxGeneralBef));
    numBothLeftAndTransBef    = length(intersect(indx1L2L,indxTransBef));
    numBothLeftAndCenterBef   = length(intersect(indx1L2L,indxCenterBef));
    numBothLeftAndGeneralBef  = length(intersect(indx1L2L,indxGeneralBef));
    numApproachAndTransBef    = length(intersect(indx1R2L,indxTransBef));
    numApproachAndCenterBef   = length(intersect(indx1R2L,indxCenterBef));
    numApproachAndGeneralBef  = length(intersect(indx1R2L,indxGeneralBef));
    numSeparateAndTransBef    = length(intersect(indx1L2R,indxTransBef));
    numSeparateAndCenterBef   = length(intersect(indx1L2R,indxCenterBef));
    numSeparateAndGeneralBef  = length(intersect(indx1L2R,indxGeneralBef));
    %state after
    numBothRightAndTransAft   = length(intersect(indx1R2R,indxTransAft));
    numBothRightAndCenterAft  = length(intersect(indx1R2R,indxCenterAft));
    numBothRightAndGeneralAft = length(intersect(indx1R2R,indxGeneralAft));
    numBothLeftAndTransAft    = length(intersect(indx1L2L,indxTransAft));
    numBothLeftAndCenterAft   = length(intersect(indx1L2L,indxCenterAft));
    numBothLeftAndGeneralAft  = length(intersect(indx1L2L,indxGeneralAft));
    numApproachAndTransAft    = length(intersect(indx1R2L,indxTransAft));
    numApproachAndCenterAft   = length(intersect(indx1R2L,indxCenterAft));
    numApproachAndGeneralAft  = length(intersect(indx1R2L,indxGeneralAft));
    numSeparateAndTransAft    = length(intersect(indx1L2R,indxTransAft));
    numSeparateAndCenterAft   = length(intersect(indx1L2R,indxCenterAft));
    numSeparateAndGeneralAft  = length(intersect(indx1L2R,indxGeneralAft));

    %get number of instances of joint sister motion AND tumbling
    numBothRightAndTumble = length(intersect(indx1R2R,indxTumble));
    numBothLeftAndTumble  = length(intersect(indx1L2L,indxTumble));
    numApproachAndTumble  = length(intersect(indx1R2L,indxTumble));
    numSeparateAndTumble  = length(intersect(indx1L2R,indxTumble));
    
    %get number of instances of joint point state AND tumbling
    %state before
    numTransBefAndTumble   = length(intersect(indxTransBef,indxTumble));
    numCenterBefAndTumble  = length(intersect(indxCenterBef,indxTumble));
    numGeneralBefAndTumble = length(intersect(indxGeneralBef,indxTumble));
    %state after
    numTransAftAndTumble   = length(intersect(indxTransAft,indxTumble));
    numCenterAftAndTumble  = length(intersect(indxCenterAft,indxTumble));
    numGeneralAftAndTumble = length(intersect(indxGeneralAft,indxTumble));
    
    %calculate conditional probabilities of sister motion IF point states
    %also IF tumbling
    %both to the right
    probBothRight = [numBothRight/numTotInstances NaN NaN NaN NaN ...
        numBothRightAndTransBef/numTransBef numBothRightAndCenterBef/numCenterBef ...
        numBothRightAndGeneralBef/numGeneralBef numBothRightAndTransAft/numTransAft ...
        numBothRightAndCenterAft/numCenterAft numBothRightAndGeneralAft/numGeneralAft ...
        numBothRightAndTumble/numTumble];
    %both to the left
    probBothLeft = [numBothLeft/numTotInstances NaN NaN NaN NaN ...
        numBothLeftAndTransBef/numTransBef numBothLeftAndCenterBef/numCenterBef ...
        numBothLeftAndGeneralBef/numGeneralBef numBothLeftAndTransAft/numTransAft ...
        numBothLeftAndCenterAft/numCenterAft numBothLeftAndGeneralAft/numGeneralAft ...
        numBothLeftAndTumble/numTumble];
    %approaching
    probApproach = [numApproach/numTotInstances NaN NaN NaN NaN ...
        numApproachAndTransBef/numTransBef numApproachAndCenterBef/numCenterBef ...
        numApproachAndGeneralBef/numGeneralBef numApproachAndTransAft/numTransAft ...
        numApproachAndCenterAft/numCenterAft numApproachAndGeneralAft/numGeneralAft ...
        numApproachAndTumble/numTumble];
    %separating
    probSeparate = [numSeparate/numTotInstances NaN NaN NaN NaN ...
        numSeparateAndTransBef/numTransBef numSeparateAndCenterBef/numCenterBef ...
        numSeparateAndGeneralBef/numGeneralBef numSeparateAndTransAft/numTransAft ...
        numSeparateAndCenterAft/numCenterAft numSeparateAndGeneralAft/numGeneralAft ...
        numSeparateAndTumble/numTumble];
    
    %calculate conditional probabilities of point state IF sister motion
    %also IF tumbling
    %transition point
    probTransBef = [numTransBef/numTotInstances numBothRightAndTransBef/numBothRight ...
        numBothLeftAndTransBef/numBothLeft numApproachAndTransBef/numApproach ...
        numSeparateAndTransBef/numSeparate NaN NaN NaN NaN NaN NaN ...
        numTransBefAndTumble/numTumble];
    probTransAft = [numTransAft/numTotInstances  numBothRightAndTransAft/numBothRight ...
        numBothLeftAndTransAft/numBothLeft numApproachAndTransAft/numApproach ...
        numSeparateAndTransAft/numSeparate NaN NaN NaN NaN NaN NaN ...
        numTransAftAndTumble/numTumble];
    %center point
    probCenterBef = [numCenterBef/numTotInstances numBothRightAndCenterBef/numBothRight ...
        numBothLeftAndCenterBef/numBothLeft numApproachAndCenterBef/numApproach ...
        numSeparateAndCenterBef/numSeparate NaN NaN NaN NaN NaN NaN ...
        numCenterBefAndTumble/numTumble];
    probCenterAft = [numCenterAft/numTotInstances  numBothRightAndCenterAft/numBothRight ...
        numBothLeftAndCenterAft/numBothLeft numApproachAndCenterAft/numApproach ...
        numSeparateAndCenterAft/numSeparate NaN NaN NaN NaN NaN NaN ...
        numCenterAftAndTumble/numTumble];
    %general point
    probGeneralBef = [numGeneralBef/numTotInstances numBothRightAndGeneralBef/numBothRight ...
        numBothLeftAndGeneralBef/numBothLeft numApproachAndGeneralBef/numApproach ...
        numSeparateAndGeneralBef/numSeparate NaN NaN NaN NaN NaN NaN ...
        numGeneralBefAndTumble/numTumble];
    probGeneralAft = [numGeneralAft/numTotInstances  numBothRightAndGeneralAft/numBothRight ...
        numBothLeftAndGeneralAft/numBothLeft numApproachAndGeneralAft/numApproach ...
        numSeparateAndGeneralAft/numSeparate NaN NaN NaN NaN NaN NaN ...
        numGeneralAftAndTumble/numTumble];
    
    %calculate conditional probability of tumbling IF sister motion
    %also IF point state
    probTumble = [numTumble/numTotInstances numBothRightAndTumble/numBothRight ...
        numBothLeftAndTumble/numBothLeft numApproachAndTumble/numApproach ...
        numSeparateAndTumble/numSeparate numTransBefAndTumble/numTransBef ...
        numCenterBefAndTumble/numCenterBef numGeneralBefAndTumble/numGeneralBef ...
        numTransAftAndTumble/numTransAft numCenterAftAndTumble/numCenterAft ...
        numGeneralAftAndTumble/numGeneralAft NaN];

    %store motion probabilities
    eval(['probBothRight' label{iLabel,1} ' = probBothRight;'])
    eval(['probBothLeft' label{iLabel,1} ' = probBothLeft;'])
    eval(['probApproach' label{iLabel,1} ' = probApproach;'])
    eval(['probSeparate' label{iLabel,1} ' = probSeparate;'])
    
    %store point state probabilities
    eval(['probTransBef' label{iLabel,1} ' = probTransBef;'])
    eval(['probCenterBef' label{iLabel,1} ' = probCenterBef;'])
    eval(['probGeneralBef' label{iLabel,1} ' = probGeneralBef;'])
    eval(['probTransAft' label{iLabel,1} ' = probTransAft;'])
    eval(['probCenterAft' label{iLabel,1} ' = probCenterAft;'])
    eval(['probGeneralAft' label{iLabel,1} ' = probGeneralAft;'])
    
    %store tumbling probability
    eval(['probTumble' label{iLabel,1} ' = probTumble;'])
    
end

%% sister coupling

%initialization
for iLabel = 1 : 3
    eval(['gammaAll' label{iLabel,1} ' = zeros(17,0);'])
end

%calculation of cross-correlation between sister displacements at lag 0
for iLabel = goodLabel

    %overall
    eval(['gammaAll = crossCorr(dispProj1Distr' label{iLabel,1} ...
        ',dispProj2Distr' label{iLabel,1} ',0);'])
    
    %both to the right
    eval(['gamma1R2R = crossCorr(dispProj1Distr' label{iLabel,1} ...
        '(indx1R2R' label{iLabel,1} '),dispProj2Distr' label{iLabel,1} ...
        '(indx1R2R' label{iLabel,1} '),0);'])
    
    %both to the left
    eval(['gamma1L2L = crossCorr(dispProj1Distr' label{iLabel,1} ...
        '(indx1L2L' label{iLabel,1} '),dispProj2Distr' label{iLabel,1} ...
        '(indx1L2L' label{iLabel,1} '),0);'])
    
    %first to the right, second to the left
    eval(['gamma1R2L = crossCorr(dispProj1Distr' label{iLabel,1} ...
        '(indx1R2L' label{iLabel,1} '),dispProj2Distr' label{iLabel,1} ...
        '(indx1R2L' label{iLabel,1} '),0);'])
    
    %first to the left, second to the right
    eval(['gamma1L2R = crossCorr(dispProj1Distr' label{iLabel,1} ...
        '(indx1L2R' label{iLabel,1} '),dispProj2Distr' label{iLabel,1} ...
        '(indx1L2R' label{iLabel,1} '),0);'])
    
    %both to the right or both to the left
    eval(['gammaSame = crossCorr(dispProj1Distr' label{iLabel,1} ...
        '([indx1R2R' label{iLabel,1} ';indx1L2L' label{iLabel,1} ...
        ']),dispProj2Distr' label{iLabel,1} '([indx1R2R' label{iLabel,1} ...
        ';indx1L2L' label{iLabel,1} ']),0);'])
    
    %one to the right and the other to the left (no matter which one)
    eval(['gammaDiff = crossCorr(dispProj1Distr' label{iLabel,1} ...
        '([indx1R2L' label{iLabel,1} ';indx1L2R' label{iLabel,1} ...
        ']),dispProj2Distr' label{iLabel,1} '([indx1R2L' label{iLabel,1} ...
        ';indx1L2R' label{iLabel,1} ']),0);'])
    
    %overall for points whose state is classified "before"
    eval(['gammaClassBef = crossCorr(dispProj1Distr' label{iLabel,1} ...
        '(indxClassBef' label{iLabel,1} '),dispProj2Distr' label{iLabel,1} ...
        '(indxClassBef' label{iLabel,1} '),0);'])
    
    %overall for points whose state is classified "after"
    eval(['gammaClassAft = crossCorr(dispProj1Distr' label{iLabel,1} ...
        '(indxClassAft' label{iLabel,1} '),dispProj2Distr' label{iLabel,1} ...
        '(indxClassAft' label{iLabel,1} '),0);'])
    
    %transition point before
    eval(['gammaTransBef = crossCorr(dispProj1Distr' label{iLabel,1} ...
        '(indxTransBef' label{iLabel,1} '),dispProj2Distr' label{iLabel,1} ...
        '(indxTransBef' label{iLabel,1} '),0);'])
    
    %transition point after
    eval(['gammaTransAft = crossCorr(dispProj1Distr' label{iLabel,1} ...
        '(indxTransAft' label{iLabel,1} '),dispProj2Distr' label{iLabel,1} ...
        '(indxTransAft' label{iLabel,1} '),0);'])
    
    %center point before
    eval(['gammaCenterBef = crossCorr(dispProj1Distr' label{iLabel,1} ...
        '(indxCenterBef' label{iLabel,1} '),dispProj2Distr' label{iLabel,1} ...
        '(indxCenterBef' label{iLabel,1} '),0);'])
    
    %center point after
    eval(['gammaCenterAft = crossCorr(dispProj1Distr' label{iLabel,1} ...
        '(indxCenterAft' label{iLabel,1} '),dispProj2Distr' label{iLabel,1} ...
        '(indxCenterAft' label{iLabel,1} '),0);'])
    
    %general point before
    eval(['gammaGeneralBef = crossCorr(dispProj1Distr' label{iLabel,1} ...
        '(indxGeneralBef' label{iLabel,1} '),dispProj2Distr' label{iLabel,1} ...
        '(indxGeneralBef' label{iLabel,1} '),0);'])
    
    %general point after
    eval(['gammaGeneralAft = crossCorr(dispProj1Distr' label{iLabel,1} ...
        '(indxGeneralAft' label{iLabel,1} '),dispProj2Distr' label{iLabel,1} ...
        '(indxGeneralAft' label{iLabel,1} '),0);'])
    
    %tumbling
    eval(['gammaTumble = crossCorr(dispProj1Distr' label{iLabel,1} ...
        '(indxTumble' label{iLabel,1} '),dispProj2Distr' label{iLabel,1} ...
        '(indxTumble' label{iLabel,1} '),0);'])
        
    %no tumbling
    eval(['gammaNoTumble = crossCorr(dispProj1Distr' label{iLabel,1} ...
        '(indxNoTumble' label{iLabel,1} '),dispProj2Distr' label{iLabel,1} ...
        '(indxNoTumble' label{iLabel,1} '),0);'])

    %store results
    eval(['gammaAll' label{iLabel,1} ' = [' ...
        'gammaAll length(indxAll' label{iLabel,1} '); ' ...
        'gamma1R2R length(indx1R2R' label{iLabel,1} '); ' ...
        'gamma1L2L length(indx1L2L' label{iLabel,1} '); ' ...
        'gamma1R2L length(indx1R2L' label{iLabel,1} '); ' ...
        'gamma1L2R length(indx1L2R' label{iLabel,1} '); ' ...
        'gammaSame length([indx1R2R' label{iLabel,1} ';indx1L2L' label{iLabel,1} ']); ' ...
        'gammaDiff length([indx1R2L' label{iLabel,1} ';indx1L2R' label{iLabel,1} ']); ' ...
        'gammaClassBef length(indxClassBef' label{iLabel,1} '); ' ...
        'gammaClassAft length(indxClassAft' label{iLabel,1} '); ' ...
        'gammaTransBef length(indxTransBef' label{iLabel,1} '); ' ...
        'gammaTransAft length(indxTransAft' label{iLabel,1} '); ' ...
        'gammaCenterBef length(indxCenterBef' label{iLabel,1} '); ' ...
        'gammaCenterAft length(indxCenterAft' label{iLabel,1} '); ' ...
        'gammaGeneralBef length(indxGeneralBef' label{iLabel,1} '); ' ...
        'gammaGeneralAft length(indxGeneralAft' label{iLabel,1} '); '...
        'gammaTumble length(indxTumble' label{iLabel,1} '); ' ...
        'gammaNoTumble length(indxNoTumble' label{iLabel,1} ')];'])

end

%% characteristics based on point state, sister motion and tumbling

%initialization
for iLabel = 1 : 3
    eval(['transPointChar' label{iLabel,1} ' = [];'])
    eval(['centerPointChar' label{iLabel,1} ' = [];'])
    eval(['generalPointChar' label{iLabel,1} ' = [];'])
end

%calculation
for iLabel = goodLabel

    %read out variables
    eval(['separationBeforeDistr = separationBeforeDistr' label{iLabel,1} ';'])
    eval(['separationAfterDistr = separationAfterDistr' label{iLabel,1} ';'])
    eval(['separationChangeDistr = separationChangeDistr' label{iLabel,1} ';'])
    eval(['angleNormBeforeDistr = angleNormBeforeDistr' label{iLabel,1} ';'])
    eval(['angleNormAfterDistr = angleNormAfterDistr' label{iLabel,1} ';'])
    eval(['angularDispDistr = angularDispDistr' label{iLabel,1} ';'])
    eval(['centerPositionChangeDistr = centerPositionChangeDistr' label{iLabel,1} ';'])
    eval(['dispProj1Distr = dispProj1Distr' label{iLabel,1} ';'])
    eval(['dispAngle1Distr = dispAngle1Distr' label{iLabel,1} ';'])
    eval(['dispProj2Distr = dispProj2Distr' label{iLabel,1} ';'])
    eval(['dispAngle2Distr = dispAngle2Distr' label{iLabel,1} ';'])
    
    %put all characteristics together
    characterBefore = [separationBeforeDistr separationChangeDistr ...
        angleNormBeforeDistr angularDispDistr centerPositionChangeDistr ...
        dispProj1Distr dispAngle1Distr dispProj2Distr dispAngle2Distr];
    characterAfter = [separationAfterDistr separationChangeDistr ...
        angleNormAfterDistr angularDispDistr centerPositionChangeDistr ...
        dispProj1Distr dispAngle1Distr dispProj2Distr dispAngle2Distr];

    %get indices
    eval(['indxAll = indxAll' label{iLabel,1} ';'])
    eval(['indx1R2R = indx1R2R' label{iLabel,1} ';'])
    eval(['indx1R2L = indx1R2L' label{iLabel,1} ';'])
    eval(['indx1L2R = indx1L2R' label{iLabel,1} ';'])
    eval(['indx1L2L = indx1L2L' label{iLabel,1} ';'])
    eval(['indxClassBef = indxClassBef' label{iLabel,1} ';'])
    eval(['indxClassAft = indxClassAft' label{iLabel,1} ';'])
    eval(['indxTransBef = indxTransBef' label{iLabel,1} ';'])
    eval(['indxTransAft = indxTransAft' label{iLabel,1} ';'])
    eval(['indxCenterBef = indxCenterBef' label{iLabel,1} ';'])
    eval(['indxCenterAft = indxCenterAft' label{iLabel,1} ';'])
    eval(['indxGeneralBef = indxGeneralBef' label{iLabel,1} ';'])
    eval(['indxGeneralAft = indxGeneralAft' label{iLabel,1} ';'])
    eval(['indxTumble = indxTumble' label{iLabel,1} ';'])
    eval(['indxNoTumble = indxNoTumble' label{iLabel,1} ';'])
    
    %get indices of transition points in the various motion and tumbling
    %categories
    %before
    indxTransBefBothRight = intersect(indxTransBef,indx1R2R);
    indxTransBefBothLeft  = intersect(indxTransBef,indx1L2L);
    indxTransBefApproach  = intersect(indxTransBef,indx1R2L);
    indxTransBefSeparate  = intersect(indxTransBef,indx1L2R);
    indxTransBefSame      = intersect(indxTransBef,[indx1R2R; indx1L2L]);
    indxTransBefOpp       = intersect(indxTransBef,[indx1R2L; indx1L2R]);
    indxTransBefTumble    = intersect(indxTransBef,indxTumble);
    indxTransBefNoTumble  = intersect(indxTransBef,indxNoTumble);
    %after
    indxTransAftBothRight = intersect(indxTransAft,indx1R2R);
    indxTransAftBothLeft  = intersect(indxTransAft,indx1L2L);
    indxTransAftApproach  = intersect(indxTransAft,indx1R2L);
    indxTransAftSeparate  = intersect(indxTransAft,indx1L2R);
    indxTransAftSame      = intersect(indxTransAft,[indx1R2R; indx1L2L]);
    indxTransAftOpp       = intersect(indxTransAft,[indx1R2L; indx1L2R]);
    indxTransAftTumble    = intersect(indxTransAft,indxTumble);
    indxTransAftNoTumble  = intersect(indxTransAft,indxNoTumble);
    
    %get indices of center points in the various motion and tumbling
    %categories
    %before
    indxCenterBefBothRight = intersect(indxCenterBef,indx1R2R);
    indxCenterBefBothLeft  = intersect(indxCenterBef,indx1L2L);
    indxCenterBefApproach  = intersect(indxCenterBef,indx1R2L);
    indxCenterBefSeparate  = intersect(indxCenterBef,indx1L2R);
    indxCenterBefSame      = intersect(indxCenterBef,[indx1R2R; indx1L2L]);
    indxCenterBefOpp       = intersect(indxCenterBef,[indx1R2L; indx1L2R]);
    indxCenterBefTumble    = intersect(indxCenterBef,indxTumble);
    indxCenterBefNoTumble  = intersect(indxCenterBef,indxNoTumble);
    %after
    indxCenterAftBothRight = intersect(indxCenterAft,indx1R2R);
    indxCenterAftBothLeft  = intersect(indxCenterAft,indx1L2L);
    indxCenterAftApproach  = intersect(indxCenterAft,indx1R2L);
    indxCenterAftSeparate  = intersect(indxCenterAft,indx1L2R);
    indxCenterAftSame      = intersect(indxCenterAft,[indx1R2R; indx1L2L]);
    indxCenterAftOpp       = intersect(indxCenterAft,[indx1R2L; indx1L2R]);
    indxCenterAftTumble    = intersect(indxCenterAft,indxTumble);
    indxCenterAftNoTumble  = intersect(indxCenterAft,indxNoTumble);
    
    %get indices of general points in the various motion and tumbling
    %categories
    %before
    indxGeneralBefBothRight = intersect(indxGeneralBef,indx1R2R);
    indxGeneralBefBothLeft  = intersect(indxGeneralBef,indx1L2L);
    indxGeneralBefApproach  = intersect(indxGeneralBef,indx1R2L);
    indxGeneralBefSeparate  = intersect(indxGeneralBef,indx1L2R);
    indxGeneralBefSame      = intersect(indxGeneralBef,[indx1R2R; indx1L2L]);
    indxGeneralBefOpp       = intersect(indxGeneralBef,[indx1R2L; indx1L2R]);
    indxGeneralBefTumble    = intersect(indxGeneralBef,indxTumble);
    indxGeneralBefNoTumble  = intersect(indxGeneralBef,indxNoTumble);
    %after
    indxGeneralAftBothRight = intersect(indxGeneralAft,indx1R2R);
    indxGeneralAftBothLeft  = intersect(indxGeneralAft,indx1L2L);
    indxGeneralAftApproach  = intersect(indxGeneralAft,indx1R2L);
    indxGeneralAftSeparate  = intersect(indxGeneralAft,indx1L2R);
    indxGeneralAftSame      = intersect(indxGeneralAft,[indx1R2R; indx1L2L]);
    indxGeneralAftOpp       = intersect(indxGeneralAft,[indx1R2L; indx1L2R]);
    indxGeneralAftTumble    = intersect(indxGeneralAft,indxTumble);
    indxGeneralAftNoTumble  = intersect(indxGeneralAft,indxNoTumble);
    
    %extract characteristics ...
    
    varName = [];
    varName{1} = 'trans';
    varName{2} = 'center';
    varName{3} = 'general';
    
    motionState = [];
    motionState{1} = 'All';
    motionState{2} = 'BothRight';
    motionState{3} = 'BothLeft';
    motionState{4} = 'Approach';
    motionState{5} = 'Separate';
    motionState{6} = 'Same';
    motionState{7} = 'Opp';
    motionState{8} = 'Tumble';
    motionState{9} = 'NoTumble';
    
    pointState = [];
    pointState{1} = 'TransBef';
    pointState{2} = 'CenterBef';
    pointState{3} = 'GeneralBef';
    
    for iState = 1 : 3

        %overall
        eval(['character2use = characterBefore(indx' pointState{iState} ',:);'])
        mean2use = mean(character2use);
        std2use = std(character2use);
        skew2use = skewness(character2use,0);
        kurt2use = kurtosis(character2use,0);
        entry2use.distribution = character2use;
        entry2use.distrParam = [mean2use; std2use; skew2use; kurt2use];
        eval([varName{iState} 'PointChar' label{iLabel,1} '.after.' motionState{1} ...
            ' = entry2use;'])
        
        %different motion states
        for iMotion = 2 : 9

            eval(['character2use = characterBefore(indx' pointState{iState} ...
                motionState{iMotion} ',:);'])
            mean2use = mean(character2use);
            std2use = std(character2use);
            skew2use = skewness(character2use,0);
            kurt2use = kurtosis(character2use,0);
            entry2use.distribution = character2use;
            entry2use.distrParam = [mean2use; std2use; skew2use; kurt2use];
            eval([varName{iState} 'PointChar' label{iLabel,1} '.after.' ...
                motionState{iMotion} ' = entry2use;'])

        end

    end
    
    pointState = [];
    pointState{1} = 'TransAft';
    pointState{2} = 'CenterAft';
    pointState{3} = 'GeneralAft';

    for iState = 1 : 3

        %overall
        eval(['character2use = characterBefore(indx' pointState{iState} ',:);'])
        mean2use = mean(character2use);
        std2use = std(character2use);
        skew2use = skewness(character2use,0);
        kurt2use = kurtosis(character2use,0);
        entry2use.distribution = character2use;
        entry2use.distrParam = [mean2use; std2use; skew2use; kurt2use];
        eval([varName{iState} 'PointChar' label{iLabel,1} '.before.' motionState{1} ...
            ' = entry2use;'])
        
        %different motion states
        for iMotion = 2 : 9

            eval(['character2use = characterBefore(indx' pointState{iState} ...
                motionState{iMotion} ',:);'])
            mean2use = mean(character2use);
            std2use = std(character2use);
            skew2use = skewness(character2use,0);
            kurt2use = kurtosis(character2use,0);
            entry2use.distribution = character2use;
            entry2use.distrParam = [mean2use; std2use; skew2use; kurt2use];
            eval([varName{iState} 'PointChar' label{iLabel,1} '.before.' ...
                motionState{iMotion} ' = entry2use;'])

        end

    end
    
end

%% output to analysisStruct

for iLabel = 1 : 3

    %store results in one structure
    eval(['numSistersCat = iGlobal' label{iLabel,1} ';'])
    eval(['distribution = struct('...
        '''separationBefore'',separationBeforeDistr' label{iLabel,1} ','...
        '''separationAfter'',separationAfterDistr' label{iLabel,1} ','...
        '''separationChange'',separationChangeDistr' label{iLabel,1} ','...
        '''angleNormBefore'',angleNormBeforeDistr' label{iLabel,1} ','...
        '''angleNormAfter'',angleNormAfterDistr' label{iLabel,1} ','...
        '''angularDisp'',angularDispDistr' label{iLabel,1} ','...
        '''centerPositionBefore'',centerPositionBeforeDistr' label{iLabel,1} ','...
        '''centerPositionAfter'',centerPositionAfterDistr' label{iLabel,1} ','...
        '''centerPositionChange'',centerPositionChangeDistr' label{iLabel,1} ','...
        '''displacement1'',displacement1Distr' label{iLabel,1} ','...
        '''dispProj1'',dispProj1Distr' label{iLabel,1} ','...
        '''dispAngle1'',dispAngle1Distr' label{iLabel,1} ','...
        '''displacement2'',displacement2Distr' label{iLabel,1} ','...
        '''dispProj2'',dispProj2Distr' label{iLabel,1} ','...
        '''dispAngle2'',dispAngle2Distr' label{iLabel,1} ','...
        '''dotProdCrossProdSign'',dotProdCrossProdSignDistr' label{iLabel,1} ','...
        '''pointStateBefore'',pointStateBeforeDistr' label{iLabel,1} ','...
        '''pointStateAfter'',pointStateAfterDistr' label{iLabel,1} ');'])
    eval(['distrParam = struct('...
        '''separationBefore'',separationBeforeParam' label{iLabel,1} ','...
        '''separationAfter'',separationAfterParam' label{iLabel,1} ','...
        '''separationChange'',separationChangeParam' label{iLabel,1} ','...
        '''angleNormBefore'',angleNormBeforeParam' label{iLabel,1} ','...
        '''angleNormAfter'',angleNormAfterParam' label{iLabel,1} ','...
        '''angularDisp'',angularDispParam' label{iLabel,1} ','...
        '''centerPositionBefore'',centerPositionBeforeParam' label{iLabel,1} ','...
        '''centerPositionAfter'',centerPositionAfterParam' label{iLabel,1} ','...
        '''centerPositionChange'',centerPositionChangeParam' label{iLabel,1} ','...
        '''displacement1'',displacement1Param' label{iLabel,1} ','...
        '''dispProj1'',dispProj1Param' label{iLabel,1} ','...
        '''dispAngle1'',dispAngle1Param' label{iLabel,1} ','...
        '''displacement2'',displacement2Param' label{iLabel,1} ','...
        '''dispProj2'',dispProj2Param' label{iLabel,1} ','...
        '''dispAngle2'',dispAngle2Param' label{iLabel,1} ');'])
    eval(['probability = struct('...
        '''bothMoveRight'',probBothRight' label{iLabel,1} ','...
        '''bothMoveLeft'',probBothLeft' label{iLabel,1} ','...
        '''bothApproach'',probApproach' label{iLabel,1} ','...
        '''bothSeparate'',probSeparate' label{iLabel,1} ','...
        '''transPointBef'',probTransBef' label{iLabel,1} ','...
        '''centerPointBef'',probCenterBef' label{iLabel,1} ','...
        '''generalPointBef'',probGeneralBef' label{iLabel,1} ','...
        '''transPointAft'',probTransAft' label{iLabel,1} ','...
        '''centerPointAft'',probCenterAft' label{iLabel,1} ','...
        '''generalPointAft'',probGeneralAft' label{iLabel,1} ','...
        '''tumbling'',probTumble' label{iLabel,1} ');'])
    eval(['crosscorr0 = struct('...
        '''overall'',gammaAll' label{iLabel,1} '(1,:),'...
        '''rightRight'',gammaAll' label{iLabel,1} '(2,:),'...
        '''leftLeft'',gammaAll' label{iLabel,1} '(3,:),'...
        '''rightLeft'',gammaAll' label{iLabel,1} '(4,:),'...
        '''leftRight'',gammaAll' label{iLabel,1} '(5,:),'...
        '''sameDir'',gammaAll' label{iLabel,1} '(6,:),'...
        '''OppDir'',gammaAll' label{iLabel,1} '(7,:),'...
        '''allPointClassBefore'',gammaAll' label{iLabel,1} '(8,:),'...
        '''allPointClassAfter'',gammaAll' label{iLabel,1} '(9,:),'...
        '''transitionBefore'',gammaAll' label{iLabel,1} '(10,:),'...
        '''transitionAfter'',gammaAll' label{iLabel,1} '(11,:),'...
        '''centerBefore'',gammaAll' label{iLabel,1} '(12,:),'...
        '''centerAfter'',gammaAll' label{iLabel,1} '(13,:),'...
        '''generalBefore'',gammaAll' label{iLabel,1} '(14,:),'...
        '''generalAfter'',gammaAll' label{iLabel,1} '(15,:),'...
        '''tumbling'',gammaAll' label{iLabel,1} '(16,:),'...
        '''noTumbling'',gammaAll' label{iLabel,1} '(17,:));'])
    eval(['characteristics = struct('...
        '''transitionPoint'',transPointChar' label{iLabel,1} ','...
        '''centerPoint'',centerPointChar' label{iLabel,1} ','...
        '''generalPoint'',generalPointChar' label{iLabel,1} ');'])
        
    eval([label{iLabel,1} ' = struct('...
        '''numSisters'',numSistersCat,'...
        '''distribution'',distribution,'...
        '''distrParam'',distrParam,'...
        '''probability'',probability,'...
        '''crosscorr0'',crosscorr0,'...
        '''characteristics'',characteristics);'])

end

breathingOscillation = struct('Inlier',Inlier,'Unaligned',Unaligned,...
    'Lagging',Lagging);

%check whether current analysisStruct already has the sisterConnection field
fieldExists = isfield(analysisStruct,'breathingOscillation');

%store results in analysisStruct
analysisStruct.breathingOscillation = breathingOscillation;

%if breathingOscillation field already existed, add 1 to the version number in
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

    %% probability stuff %%
    
    for iLabel = goodLabel

        %open figure and write title
        figure('Name',[fileName(1:end-4) ' - State probabilities - ' label{iLabel,1}],'NumberTitle','off');

        %subplot 1 - both move to right
        subplot(3,4,1);
        eval(['bar((0:11),probBothRight' label{iLabel,1} ',''FaceColor'',[0 1 1])';])
        axis([-1 12 0 1]);
        title('1: Right-Right');

        %subplot 2 - both move to left
        subplot(3,4,2);
        eval(['bar((0:11),probBothLeft' label{iLabel,1} ',''FaceColor'',[0 1 1])';])
        axis([-1 12 0 1]);
        title('2: Left-Left');

        %subplot 3 - both approach
        subplot(3,4,3);
        eval(['bar((0:11),probApproach' label{iLabel,1} ',''FaceColor'',[0 1 1])';])
        axis([-1 12 0 1]);
        title('3: Right-Left');

        %subplot 4 - both separate
        subplot(3,4,4);
        eval(['bar((0:11),probSeparate' label{iLabel,1} ',''FaceColor'',[0 1 1])';])
        axis([-1 12 0 1]);
        title('4: Left-Right');

        %subplot 5 - transition point "before"
        subplot(3,4,5);
        eval(['bar((0:11),probTransBef' label{iLabel,1} ',''FaceColor'',[0 1 1])';])
        axis([-1 12 0 1]);
        title('5: Trans. point before');

        %subplot 6 - center point "before"
        subplot(3,4,6);
        eval(['bar((0:11),probCenterBef' label{iLabel,1} ',''FaceColor'',[0 1 1])';])
        axis([-1 12 0 1]);
        title('6: Center point before');

        %subplot 7 - general point "before"
        subplot(3,4,7);
        eval(['bar((0:11),probGeneralBef' label{iLabel,1} ',''FaceColor'',[0 1 1])';])
        axis([-1 12 0 1]);
        title('7: Gen. point before');

        %subplot 8 - nothing

        %subplot 9 - transition point "after"
        subplot(3,4,9);
        eval(['bar((0:11),probTransAft' label{iLabel,1} ',''FaceColor'',[0 1 1])';])
        axis([-1 12 0 1]);
        title('8: Trans. point after');

        %subplot 10 - center point "after"
        subplot(3,4,10);
        eval(['bar((0:11),probCenterAft' label{iLabel,1} ',''FaceColor'',[0 1 1])';])
        axis([-1 12 0 1]);
        title('9: Center point after');

        %subplot 11 - general point "after"
        subplot(3,4,11);
        eval(['bar((0:11),probGeneralAft' label{iLabel,1} ',''FaceColor'',[0 1 1])';])
        axis([-1 12 0 1]);
        title('10: Gen. point after');

        %subplot 12 - tumbling
        subplot(3,4,12);
        eval(['bar((0:11),probTumble' label{iLabel,1} ',''FaceColor'',[0 1 1])';])
        axis([-1 12 0 1]);
        title('11: Tumbling');

    end
    
    %% cross-correlation stuff %%
    
    for iLabel = goodLabel

        %open figure and write title
        figure('Name',[fileName(1:end-4) ' - Sister coupling at lag 0 - ' label{iLabel,1}],'NumberTitle','off');

        %put x-axis labels
        axes(...
            'CameraUpVector',[0 1 0],...
            'XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17],...
            'XTickLabel',{'overall','r.r.','l.l.','r.l.','l.r.','r.r.+l.l.','r.l+l.r.','all b.','all a.','t.b.','t.a.','c.b.','c.a.','g.b.','g.a.','tumble','no tumble'});
        hold('all');

        %plot cross-correlations
        eval(['bar(gammaAll' label{iLabel,1} '(:,1),''FaceColor'',[0 1 1])'])

        %add error bars
        eval(['myErrorbar(gammaAll' label{iLabel,1} '(:,1),gammaAll' label{iLabel,1} '(:,2))'])
        
        %label axes
        xlabel('Category');
        ylabel('Cross-correlation at lag 0');

    end

    %% characteristics stuff %%

    for iLabel = 1

        %assign figure names
        figureName{1} = 'Sister separation';
        figureName{2} = 'Sister separation change';
        figureName{3} = 'Angle with normal';
        figureName{4} = 'Angular displacement';
        figureName{5} = 'Center position change';
        figureName{6} = 'Sister 1 displacement projection';
        figureName{7} = 'Sister 1 displacement angle';
        figureName{8} = 'Sister 2 displacement projection';
        figureName{9} = 'Sister 2 displacement angle';
        
        %assign point states
        pointState = [];
        pointState{1} = 'trans';
        pointState{2} = 'center';
        pointState{3} = 'general';
        
        for iColumn = 1 : 9

            %open figure and write title
            figure('Name',[fileName(1:end-4) ' - ' figureName{iColumn} ...
                ' - ' label{iLabel,1}],'NumberTitle','off');

            for iState = 1 : 3

                subplot(6,9,(iState-1)*3+1)
                eval(['histogram(' pointState{iState} 'PointChar' ...
                    label{iLabel,1} '.after.All.distribution(:,iColumn),[],[],0)'])
                xlabel('overall','Fontsize',8)
                h = gca; set(h,'Fontsize',6);
                subplot(6,9,(iState-1)*3+2)
                eval(['histogram(' pointState{iState} 'PointChar' ...
                    label{iLabel,1} '.after.Same.distribution(:,iColumn),[],[],0)'])
                title(['After ' pointState{iState} ' point'],'Fontsize',12)
                xlabel('r.r.+l.l.','Fontsize',8);
                h = gca; set(h,'Fontsize',6);
                subplot(6,9,(iState-1)*3+3)
                eval(['histogram(' pointState{iState} 'PointChar' ...
                    label{iLabel,1} '.after.Opp.distribution(:,iColumn),[],[],0)'])
                xlabel('r.l.+l.r.','Fontsize',8)
                h = gca; set(h,'Fontsize',6);
                
                subplot(6,9,(iState-1)*3+10)
                eval(['histogram(' pointState{iState} 'PointChar' ...
                    label{iLabel,1} '.after.BothRight.distribution(:,iColumn),[],[],0)'])
                xlabel('r.r.','Fontsize',8)
                h = gca; set(h,'Fontsize',6);
                subplot(6,9,(iState-1)*3+11)
                eval(['histogram(' pointState{iState} 'PointChar' ...
                    label{iLabel,1} '.after.BothLeft.distribution(:,iColumn),[],[],0)'])
                xlabel('l.l.','Fontsize',8)
                h = gca; set(h,'Fontsize',6);
                subplot(6,9,(iState-1)*3+12)
                eval(['histogram(' pointState{iState} 'PointChar' ...
                    label{iLabel,1} '.after.Approach.distribution(:,iColumn),[],[],0)'])
                xlabel('r.l.','Fontsize',8)
                h = gca; set(h,'Fontsize',6);
                
                subplot(6,9,(iState-1)*3+19)
                eval(['histogram(' pointState{iState} 'PointChar' ...
                    label{iLabel,1} '.after.Separate.distribution(:,iColumn),[],[],0)'])
                xlabel('l.r.','Fontsize',8)
                h = gca; set(h,'Fontsize',6);
                subplot(6,9,(iState-1)*3+20)
                eval(['histogram(' pointState{iState} 'PointChar' ...
                    label{iLabel,1} '.after.Tumble.distribution(:,iColumn),[],[],0)'])
                xlabel('tumble','Fontsize',8)
                h = gca; set(h,'Fontsize',6);
                subplot(6,9,(iState-1)*3+21)
                eval(['histogram(' pointState{iState} 'PointChar' ...
                    label{iLabel,1} '.after.NoTumble.distribution(:,iColumn),[],[],0)'])
                xlabel('no tumble','Fontsize',8)
                h = gca; set(h,'Fontsize',6);
                
                subplot(6,9,(iState-1)*3+28)
                eval(['histogram(' pointState{iState} 'PointChar' ...
                    label{iLabel,1} '.before.All.distribution(:,iColumn),[],[],0)'])
                xlabel('overall','Fontsize',8)
                h = gca; set(h,'Fontsize',6);
                subplot(6,9,(iState-1)*3+29)
                eval(['histogram(' pointState{iState} 'PointChar' ...
                    label{iLabel,1} '.before.Same.distribution(:,iColumn),[],[],0)'])
                title(['Before ' pointState{iState} ' point'],'Fontsize',12)
                xlabel('r.r.+l.l.','Fontsize',8);
                h = gca; set(h,'Fontsize',6);
                subplot(6,9,(iState-1)*3+30)
                eval(['histogram(' pointState{iState} 'PointChar' ...
                    label{iLabel,1} '.before.Opp.distribution(:,iColumn),[],[],0)'])
                xlabel('r.l.+l.r.','Fontsize',8)
                h = gca; set(h,'Fontsize',6);
                
                subplot(6,9,(iState-1)*3+37)
                eval(['histogram(' pointState{iState} 'PointChar' ...
                    label{iLabel,1} '.before.BothRight.distribution(:,iColumn),[],[],0)'])
                xlabel('r.r.','Fontsize',8)
                h = gca; set(h,'Fontsize',6);
                subplot(6,9,(iState-1)*3+38)
                eval(['histogram(' pointState{iState} 'PointChar' ...
                    label{iLabel,1} '.before.BothLeft.distribution(:,iColumn),[],[],0)'])
                xlabel('l.l.','Fontsize',8)
                h = gca; set(h,'Fontsize',6);
                subplot(6,9,(iState-1)*3+39)
                eval(['histogram(' pointState{iState} 'PointChar' ...
                    label{iLabel,1} '.before.Approach.distribution(:,iColumn),[],[],0)'])
                xlabel('r.l.','Fontsize',8)
                h = gca; set(h,'Fontsize',6);
                
                subplot(6,9,(iState-1)*3+46)
                eval(['histogram(' pointState{iState} 'PointChar' ...
                    label{iLabel,1} '.before.Separate.distribution(:,iColumn),[],[],0)'])
                xlabel('l.r.','Fontsize',8)
                h = gca; set(h,'Fontsize',6);
                subplot(6,9,(iState-1)*3+47)
                eval(['histogram(' pointState{iState} 'PointChar' ...
                    label{iLabel,1} '.before.Tumble.distribution(:,iColumn),[],[],0)'])
                xlabel('tumble','Fontsize',8)
                h = gca; set(h,'Fontsize',6);
                subplot(6,9,(iState-1)*3+48)
                eval(['histogram(' pointState{iState} 'PointChar' ...
                    label{iLabel,1} '.before.NoTumble.distribution(:,iColumn),[],[],0)'])
                xlabel('no tumble','Fontsize',8)
                h = gca; set(h,'Fontsize',6);

            end %(for iState = 1 : 3)
            
        end %(for iColumn = 1 : 7)

    end %(for iLabel = goodLabel)

end

%% ~~~ the end ~~~ %%

%% old stuff

% %% auto- and cross-correlations
% 
% %initialization
% for iLabel = 1 : 3
%     eval(['sisterSepAutocorr' label{iLabel,1} ' = [];'])
%     eval(['sisterCoMAutocorr' label{iLabel,1} ' = [];'])
%     eval(['sepCoMPosCrosscorr' label{iLabel,1} ' = [];'])
%     eval(['sepAbsCoMPosCrosscorr' label{iLabel,1} ' = [];'])
% end
% 
% %define maximum lag
% maxLag = 10;
% 
% %calculation
% for iLabel = goodLabel
% 
%     %autocorrelation of change in sister separation
%     eval(['sisterSepAutocorr' label{iLabel,1} ' = autoCorr(sisterSepChange' ...
%         label{iLabel,1} ',maxLag);'])
% 
%     %autocorrelation of center of mass displacement
%     eval(['sisterCoMAutocorr' label{iLabel,1} ' = autoCorr(sisterCoMPosChange' ...
%         label{iLabel,1} ',maxLag);'])
% 
%     %crosscorrelation between separation and abs(center of mass position
%     %change)
%     eval(['series1 = sisterCoMPosChange' label{iLabel,1} ';'])
%     for i=1:length(series1)
%         series1(i).observations(:,1) = series1(i).observations(:,1) - ...
%             nanmean(series1(i).observations(:,1));
%         series1(i).observations(:,1) = abs(series1(i).observations(:,1));
%     end
%     eval(['series2 = sisterSep' label{iLabel,1} ';'])
%     series2bef = [];
%     series2aft = [];
%     for i=1:length(series2)
%         series2(i).observations(:,1) = series2(i).observations(:,1) - ...
%             nanmean(series2(i).observations(:,1));
%         series2bef(i).observations = series2(i).observations(1:end-1,:);
%         series2aft(i).observations = series2(i).observations(2:end,:);
%     end
%     eval(['sepBefCoMPosChangeCrosscorr' label{iLabel,1} ...
%         ' = crossCorr(series1,series2bef,maxLag);'])
%     eval(['sepAftCoMPosChangeCrosscorr' label{iLabel,1} ...
%         ' = crossCorr(series1,series2aft,maxLag);'])
% 
% end

%     
%     %test whether the distributions change from one bin to the next ...
%     
%     
%     
% end

%     eval(['autocorr = struct(''sisterSepChange'',sisterSepAutocorr' label{iLabel,1} ...
%         ',''sisterCoMPosChange'',sisterCoMAutocorr' label{iLabel,1} ');']);
%     eval(['crosscorr = struct(''separationCoMPosition'',sepCoMPosCrosscorr' label{iLabel,1} ','...
%         '''separationAbsCoMPosition'',sepAbsCoMPosCrosscorr' label{iLabel,1} ');']);

%     %% crosscorrelation stuff %%
%     
%     for iLabel = goodLabel
% 
%         %open figure and write title
%         figure('Name',[fileName(1:end-4) ' - Breathing-oscillation crosscorrelation - ' label{iLabel,1}],'NumberTitle','off');
% 
%         %create subplot 1
%         subplot(2,2,1);
%         hold on;
% 
%         %plot sister separation change autocorrelation
%         eval(['autocorr1 = sisterSepAutocorr' label{iLabel,1} ';']);
%         plot((0:maxLag)*timeLapse,autocorr1(:,1),'k','marker','.');
% 
%         %set axes limit
%         minVal = min(autocorr1(:,1));
%         axis([0 maxLag*timeLapse min(0,1.1*minVal) 1.1]);
% 
%         %write axes labels
%         xlabel('Lag (s)');
%         ylabel('Autocorrelation of change in sister separation');
% 
%         %hold off figure
%         hold off
% 
%         %create subplot 2
%         subplot(2,2,2);
%         hold on;
% 
%         %plot center of mass position change autocorrelation
%         eval(['autocorr2 = sisterCoMAutocorr' label{iLabel,1} ';']);
%         plot((0:maxLag)*timeLapse,autocorr2(:,1),'k','marker','.');
% 
%         %set axes limit
%         minVal = min(autocorr2(:,1));
%         axis([0 maxLag*timeLapse min(0,1.1*minVal) 1.1]);
% 
%         %write axes labels
%         xlabel('Lag (s)');
%         ylabel('Autocorrelation of sister center displacement');
% 
%         %hold off figure
%         hold off
% 
%         %create subplot 3
%         subplot(2,2,3);
%         hold on;
% 
%         %plot crosscorrelation between sister separation and center of mass
%         %position
%         eval(['crosscorr1 = sepCoMPosCrosscorr' label{iLabel,1} ';']);
%         plot((-maxLag:maxLag)*timeLapse,crosscorr1(:,1),'k','marker','.');
% 
%         %set axes limit
%         minVal = min(crosscorr1(:,1));
%         maxVal = max(crosscorr1(:,1));
%         axis([-maxLag*timeLapse maxLag*timeLapse min(0,1.1*minVal) ...
%             max(0,1.1*maxVal)]);
% 
%         %write axes labels
%         xlabel('Lag (s)');
%         ylabel('crosscorrelation between sister separation and center position');
% 
%         %hold off figure
%         hold off
% 
%         %create subplot 3
%         subplot(2,2,4);
%         hold on;
% 
%         %plot crosscorrelation between sister separation and absolute
%         %values of center of mass position (where 0 is average center of
%         %mass position)
%         eval(['crosscorr2 = sepAbsCoMPosCrosscorr' label{iLabel,1} ';']);
%         plot((-maxLag:maxLag)*timeLapse,crosscorr2(:,1),'k','marker','.');
% 
%         %set axes limit
%         minVal = min(crosscorr2(:,1));
%         maxVal = max(crosscorr2(:,1));
%         axis([-maxLag*timeLapse maxLag*timeLapse min(0,1.1*minVal) ...
%             max(0,1.1*maxVal)]);
% 
%         %write axes labels
%         xlabel('Lag (s)');
%         ylabel('crosscorrelation between sister separation and absolute center position');
% 
%         %hold off figure
%         hold off
% 
%     end

%     %% separation and velocity vs. position stuff %%
% 
%     for iLabel = goodLabel
% 
%         %open figure and write title
%         figure('Name',[fileName(1:end-4) ' - Separation & velocity vs. position - ' label{iLabel,1}],'NumberTitle','off');
% 
%         %create subplot 1
%         subplot(2,3,1);
%         hold on;
% 
%         %plot sister separation vs. position, taking an equal number of
%         %observations per bin
%         eval(['values = comPosSisSep' label{iLabel,1} ';'])
%         numObsPerBin = vertcat(values.numObserve);
%         numObsPerBinSort = sort(numObsPerBin);
%         numObs2use = numObsPerBinSort(find(numObsPerBinSort>=50,1,'first'));
%         for i=1:length(numObsPerBin)
%             if numObsPerBin(i) >= numObs2use
%                 randIndx = randsample(numObsPerBin(i),numObs2use);
%                 val2plot = values(i).observations(randIndx,:);
%                 plot(val2plot(val2plot(:,3)==3,1),...
%                     val2plot(val2plot(:,3)==3,2),'b.');
%                 plot(val2plot(val2plot(:,3)==-3,1),...
%                     val2plot(val2plot(:,3)==-3,2),'k.');
%                 plot(val2plot(val2plot(:,3)==1,1),...
%                     val2plot(val2plot(:,3)==1,2),'r.');
%                 plot(val2plot(val2plot(:,3)==-1,1),...
%                     val2plot(val2plot(:,3)==-1,2),'g.');
%                 plot(val2plot(val2plot(:,3)==2,1),...
%                     val2plot(val2plot(:,3)==2,2),'m.');
%                 plot(val2plot(val2plot(:,3)==-2,1),...
%                     val2plot(val2plot(:,3)==-2,2),'c.');
%             else
%                 val2plot = values(i).observations;
%                 plot(val2plot(val2plot(:,3)==3,1),...
%                     val2plot(val2plot(:,3)==3,2),'b+');
%                 plot(val2plot(val2plot(:,3)==-3,1),...
%                     val2plot(val2plot(:,3)==-3,2),'k+');
%                 plot(val2plot(val2plot(:,3)==1,1),...
%                     val2plot(val2plot(:,3)==1,2),'r+');
%                 plot(val2plot(val2plot(:,3)==-1,1),...
%                     val2plot(val2plot(:,3)==-1,2),'g+');
%                 plot(val2plot(val2plot(:,3)==2,1),...
%                     val2plot(val2plot(:,3)==2,2),'m+');
%                 plot(val2plot(val2plot(:,3)==-2,1),...
%                     val2plot(val2plot(:,3)==-2,2),'c+');
%             end
%         end
%         
%         %write axes labels
%         xlabel('Center of mass position');
%         ylabel('Sister separation at particular position');
% 
%         %hold off figure
%         hold off
% 
%         %create subplot 2
%         subplot(2,3,2);
%         hold on;
% 
%         %plot sister "after" velocity vs. position, taking an equal number of
%         %observations per bin
%         eval(['values = comPosSisVelAft' label{iLabel,1} ';'])
%         numObsPerBin = vertcat(values.numObserve);
%         numObsPerBinSort = sort(numObsPerBin);
%         numObs2use = numObsPerBinSort(find(numObsPerBinSort>=50,1,'first'));
%         for i=1:length(numObsPerBin)
%             if numObsPerBin(i) >= numObs2use
%                 randIndx = randsample(numObsPerBin(i),numObs2use);
%                 val2plot = values(i).observations(randIndx,:);
%                 plot(val2plot(val2plot(:,3)==3,1),...
%                     val2plot(val2plot(:,3)==3,2),'b.');
%                 plot(val2plot(val2plot(:,3)==-3,1),...
%                     val2plot(val2plot(:,3)==-3,2),'k.');
%                 plot(val2plot(val2plot(:,3)==1,1),...
%                     val2plot(val2plot(:,3)==1,2),'r.');
%                 plot(val2plot(val2plot(:,3)==-1,1),...
%                     val2plot(val2plot(:,3)==-1,2),'g.');
%                 plot(val2plot(val2plot(:,3)==2,1),...
%                     val2plot(val2plot(:,3)==2,2),'m.');
%                 plot(val2plot(val2plot(:,3)==-2,1),...
%                     val2plot(val2plot(:,3)==-2,2),'c.');
%             else
%                 val2plot = values(i).observations;
%                 plot(val2plot(val2plot(:,3)==3,1),...
%                     val2plot(val2plot(:,3)==3,2),'b+');
%                 plot(val2plot(val2plot(:,3)==-3,1),...
%                     val2plot(val2plot(:,3)==-3,2),'k+');
%                 plot(val2plot(val2plot(:,3)==1,1),...
%                     val2plot(val2plot(:,3)==1,2),'r+');
%                 plot(val2plot(val2plot(:,3)==-1,1),...
%                     val2plot(val2plot(:,3)==-1,2),'g+');
%                 plot(val2plot(val2plot(:,3)==2,1),...
%                     val2plot(val2plot(:,3)==2,2),'m+');
%                 plot(val2plot(val2plot(:,3)==-2,1),...
%                     val2plot(val2plot(:,3)==-2,2),'c+');
%             end
%         end
% 
%         %write axes labels
%         xlabel('Center of mass position');
%         ylabel('Change in sister separation after particular position');
% 
%         %hold off figure
%         hold off
% 
%         %create subplot 3
%         subplot(2,3,3);
%         hold on;
% 
%         %plot sister "before" velocity vs. position, taking an equal number of
%         %observations per bin
%         eval(['values = comPosSisVelBef' label{iLabel,1} ';'])
%         numObsPerBin = vertcat(values.numObserve);
%         numObsPerBinSort = sort(numObsPerBin);
%         numObs2use = numObsPerBinSort(find(numObsPerBinSort>=50,1,'first'));
%         for i=1:length(numObsPerBin)
%             if numObsPerBin(i) >= numObs2use
%                 randIndx = randsample(numObsPerBin(i),numObs2use);
%                 val2plot = values(i).observations(randIndx,:);
%                 plot(val2plot(val2plot(:,3)==3,1),...
%                     val2plot(val2plot(:,3)==3,2),'b.');
%                 plot(val2plot(val2plot(:,3)==-3,1),...
%                     val2plot(val2plot(:,3)==-3,2),'k.');
%                 plot(val2plot(val2plot(:,3)==1,1),...
%                     val2plot(val2plot(:,3)==1,2),'r.');
%                 plot(val2plot(val2plot(:,3)==-1,1),...
%                     val2plot(val2plot(:,3)==-1,2),'g.');
%                 plot(val2plot(val2plot(:,3)==2,1),...
%                     val2plot(val2plot(:,3)==2,2),'m.');
%                 plot(val2plot(val2plot(:,3)==-2,1),...
%                     val2plot(val2plot(:,3)==-2,2),'c.');
%             else
%                 val2plot = values(i).observations;
%                 plot(val2plot(val2plot(:,3)==3,1),...
%                     val2plot(val2plot(:,3)==3,2),'b+');
%                 plot(val2plot(val2plot(:,3)==-3,1),...
%                     val2plot(val2plot(:,3)==-3,2),'k+');
%                 plot(val2plot(val2plot(:,3)==1,1),...
%                     val2plot(val2plot(:,3)==1,2),'r+');
%                 plot(val2plot(val2plot(:,3)==-1,1),...
%                     val2plot(val2plot(:,3)==-1,2),'g+');
%                 plot(val2plot(val2plot(:,3)==2,1),...
%                     val2plot(val2plot(:,3)==2,2),'m+');
%                 plot(val2plot(val2plot(:,3)==-2,1),...
%                     val2plot(val2plot(:,3)==-2,2),'c+');
%             end
%         end
% 
%         %write axes labels
%         xlabel('Center of mass position');
%         ylabel('Change in sister separation before particular position');
% 
%         %hold off figure
%         hold off
% 
%         %do some calculations for bar graphs below
%         fracPlus1  = NaN(length(numObsPerBin),1);
%         fracMinus1 = fracPlus1;
%         binMean = fracPlus1;
%         for i=1:length(numObsPerBin)
%             fracPlus1(i)  = length(find(values(i).observations(:,3)==1)) ...
%                 / numObsPerBin(i); %pos. to neg. transitions
%             fracMinus1(i) = length(find(values(i).observations(:,3)==-1)) ...
%                 / numObsPerBin(i); %neg. to pos. transitions
%             binMean(i) = mean(values(i).minMax); %mean bin value
%         end
% 
%         %create subplot 4
%         subplot(2,3,4);
%         hold on;
%         
%         %make bar graph of number of observations per bin
%         bar(binMean,numObsPerBin);
% 
%         %write axes labels
%         xlabel('Center of mass position');
%         ylabel('Number of observations');
%         
%         %hold off figure
%         hold off
% 
%         %create subplot 5
%         subplot(2,3,5);
%         hold on;
%         
%         %make bar graph of ratio of number of pos. to neg. transition
%         %points to total number of points
%         bar(binMean,fracPlus1);
% 
%         %write axes labels
%         xlabel('Center of mass position');
%         ylabel('Frac pos to neg transition points');
% 
%         %hold off figure
%         hold off
% 
%         %create subplot 6
%         subplot(2,3,6);
%         hold on;
%         
%         %make bar graph of ratio of number of neg. to pos. transition
%         %points to total number of points
%         bar(binMean,fracMinus1);
% 
%         %write axes labels
%         xlabel('Center of mass position');
%         ylabel('Frac neg to pos transition points');
% 
%         %hold off figure
%         hold off
% 
%     end
% 
% end
% 
% %% distributions per oscillation state
% 
% %initialization
% for iLabel = 1 : 3
%     eval(['sisterSepDistr' label{iLabel,1} ' = [];'])
%     eval(['sisterSepChangeAftDistr' label{iLabel,1} ' = [];'])
%     eval(['sisterSepChangeBefDistr' label{iLabel,1} ' = [];'])
%     eval(['comPosDistr' label{iLabel,1} ' = [];'])
%     eval(['comPosChangeAftDistr' label{iLabel,1} ' = [];'])
%     eval(['comPosChangeBefDistr' label{iLabel,1} ' = [];'])
% end
% 
% %specify point states to look for
% pointStateLook4 = [-3 -2 -1 1 2 3];
%
% %go over all available sister categories
% for iLabel = goodLabel
%
%     %initialize
%     sisterSepDistr = repmat(struct('observations',[]),7,1);
%     sisterSepChangeAftDistr = repmat(struct('observations',[]),7,1);
%     sisterSepChangeBefDistr = repmat(struct('observations',[]),7,1);
%     comPosDistr = repmat(struct('observations',[]),7,1);
%     comPosChangeAftDistr = repmat(struct('observations',[]),7,1);
%     comPosChangeBefDistr = repmat(struct('observations',[]),7,1);
%
%     %go over all sisters
%     for iSister = 1 : eval(['iGlobal' label{iLabel,1}])
%
%         %get sister information
%         eval(['pointState = pointState' label{iLabel,1} '(iSister).values;'])
%         eval(['sisterSep = sisterSep' label{iLabel,1} '(iSister).observations;'])
%         eval(['sisterSepChange = sisterSepChange' label{iLabel,1} '(iSister).observations;'])
%         eval(['comPos = sisterCoMPos' label{iLabel,1} '(iSister).observations;'])
%         eval(['comPosChange = sisterCoMPosChange' label{iLabel,1} '(iSister).observations;'])
%         
%         %go over point states and collect variable values
%         for iState = pointStateLook4
%             indxGood = find( pointState == iState );
%             sisterSepDistr(iState+4).observations = ...
%                 [sisterSepDistr(iState+4).observations; ...
%                 sisterSep(indxGood,:)];
%             sisterSepChangeAftDistr(iState+4).observations = ...
%                 [sisterSepChangeAftDistr(iState+4).observations; ...
%                 sisterSepChange(indxGood,:)];
%             sisterSepChangeBefDistr(iState+4).observations = ...
%                 [sisterSepChangeBefDistr(iState+4).observations; ...
%                 sisterSepChange(indxGood-1,:)];
%             comPosDistr(iState+4).observations = ...
%                 [comPosDistr(iState+4).observations; ...
%                 comPos(indxGood,:)];
%             comPosChangeAftDistr(iState+4).observations = ...
%                 [comPosChangeAftDistr(iState+4).observations; ...
%                 comPosChange(indxGood,:)];
%             comPosChangeBefDistr(iState+4).observations = ...
%                 [comPosChangeBefDistr(iState+4).observations; ...
%                 comPosChange(indxGood-1,:)];
%         end
% 
%     end
% 
%     %store the information for this category of sisters
%     eval(['sisterSepDistr' label{iLabel,1} ' = sisterSepDistr;'])
%     eval(['sisterSepChangeAftDistr' label{iLabel,1} ' = sisterSepChangeAftDistr;'])
%     eval(['sisterSepChangeBefDistr' label{iLabel,1} ' = sisterSepChangeBefDistr;'])
%     eval(['comPosDistr' label{iLabel,1} ' = comPosDistr;'])
%     eval(['comPosChangeAftDistr' label{iLabel,1} ' = comPosChangeAftDistr;'])
%     eval(['comPosChangeBefDistr' label{iLabel,1} ' = comPosChangeBefDistr;'])
% 
% end
% 
% %% separation & change vs. position & change
% 
% %initialization
% for iLabel = 1 : 3
%     eval(['comPosSisSep' label{iLabel,1} ' = [];'])
%     eval(['comPosSisVelAft' label{iLabel,1} ' = [];'])
%     eval(['comPosSisVelBef' label{iLabel,1} ' = [];'])
%     eval(['sisSepComVelAft' label{iLabel,1} ' = [];'])
%     eval(['sisSepComVelBef' label{iLabel,1} ' = [];'])
%     eval(['comVelSisVel' label{iLabel,1} ' = [];'])
%     eval(['comPosSisSepAft' label{iLabel,1} ' = [];'])
%     eval(['comVelSisVelAft' label{iLabel,1} ' = [];'])
% end
% 
% %calculation
% for iLabel = goodLabel
% 
%     %get this sister category's information
%     eval(['comPos = sisterCoMPos' label{iLabel,1} ';'])
%     eval(['comVel = sisterCoMPosChange' label{iLabel,1} ';'])
%     eval(['sisSep = sisterSep' label{iLabel,1} ';'])
%     eval(['sisVel = sisterSepChange' label{iLabel,1} ';'])
%     eval(['pointState = pointState' label{iLabel,1} ';'])
%     
%     %make vector of center positions, sister separations, and point states
%     comPosAll = vertcat(comPos.observations);
%     sisSepAll = vertcat(sisSep.observations);
%     pointStateAll = vertcat(pointState.values);
%     availPoints = find(~isnan(sisSepAll(:,1)));
%     comPosSisSep = [comPosAll(availPoints,1) sisSepAll(availPoints,1) ...
%         pointStateAll(availPoints,1)];
% 
%     %make vector of center positions, sister separation changes 'after',
%     %and point states
%     comPosAll = [];
%     sisVelAll = [];
%     pointStateAll = [];
%     for i = 1 : eval(['iGlobal' label{iLabel,1}])
%         comPosAll = [comPosAll; comPos(i).observations(1:end-1,1)];
%         sisVelAll = [sisVelAll; sisVel(i).observations(:,1)];
%         pointStateAll = [pointStateAll; pointState(i).values(1:end-1,1)];
%     end
%     availPoints = find(~isnan(sisVelAll(:,1)));
%     comPosSisVelAft = [comPosAll(availPoints,1) sisVelAll(availPoints,1) ...
%         pointStateAll(availPoints,1)];
% 
%     %make vector of center positions, sister separation changes 'before',
%     %and point states
%     comPosAll = [];
%     pointStateAll = [];
%     for i = 1 : eval(['iGlobal' label{iLabel,1}])
%         comPosAll = [comPosAll; comPos(i).observations(2:end,1)];
%         pointStateAll = [pointStateAll; pointState(i).values(2:end,1)];
%     end
%     comPosSisVelBef = [comPosAll(availPoints,1) sisVelAll(availPoints,1) ...
%         pointStateAll(availPoints,1)];
%     
%     %make vector of sister separations, center position changes 'after',
%     %and point states
%     sisSepAll = [];
%     comVelAll = [];
%     pointStateAll = [];
%     for i = 1 : eval(['iGlobal' label{iLabel,1}])
%         sisSepAll = [sisSepAll; sisSep(i).observations(1:end-1,1)];
%         comVelAll = [comVelAll; comVel(i).observations(:,1)];
%         pointStateAll = [pointStateAll; pointState(i).values(1:end-1,1)];
%     end
%     availPoints = find(~isnan(comVelAll(:,1)));
%     sisSepComVelAft = [sisSepAll(availPoints,1) comVelAll(availPoints,1) ...
%         pointStateAll(availPoints,1)];
%         
%     %make vector of sister separation, center position changes 'before',
%     %and point states
%     sisSepAll = [];
%     pointStateAll = [];
%     for i = 1 : eval(['iGlobal' label{iLabel,1}])
%         sisSepAll = [sisSepAll; sisSep(i).observations(2:end,1)];
%         pointStateAll = [pointStateAll; pointState(i).values(2:end,1)];
%     end
%     sisSepComVelBef = [sisSepAll(availPoints,1) comVelAll(availPoints,1) ...
%         pointStateAll(availPoints,1)];
%     
%     %make vector of sister separation changes, center position changes,
%     %and point states
%     comVelAll = [];
%     sisVelAll = [];
%     pointStateAll = [];
%     for i = 1 : eval(['iGlobal' label{iLabel,1}])
%         comVelAll = [comVelAll; comVel(i).observations(:,1)];
%         sisVelAll = [sisVelAll; sisVel(i).observations(:,1)];
%         pointStateAll = [pointStateAll; pointState(i).values(1:end-1,1)];
%     end
%     availPoints = find(~isnan(sisVelAll(:,1)));
%     comVelSisVel = [comVelAll(availPoints,1) sisVelAll(availPoints,1) ...
%         pointStateAll(availPoints,1)];
% 
%     %make vector of center positions, sister separations 'after',
%     %and point states
%     comPosAll = [];
%     sisSepAll = [];
%     pointStateAll = [];
%     for i = 1 : eval(['iGlobal' label{iLabel,1}])
%         comPosAll = [comPosAll; comPos(i).observations(1:end-1,1)];
%         sisSepAll = [sisSepAll; sisSep(i).observations(2:end,1)];
%         pointStateAll = [pointStateAll; pointState(i).values(1:end-1,1)];
%     end
%     availPoints = find( ~isnan(comPosAll(:,1)) & ~isnan(sisSepAll(:,1)) );
%     comPosSisSepAft = [comPosAll(availPoints,1) sisSepAll(availPoints,1) ...
%         pointStateAll(availPoints,1)];
%         
%     %make vector of center position changes, sister separation changes
%     %'after', and point states
%     comVelAll = [];
%     sisVelAll = [];
%     pointStateAll = [];
%     for i = 1 : eval(['iGlobal' label{iLabel,1}])
%         comVelAll = [comVelAll; comVel(i).observations(1:end-1,1)];
%         sisVelAll = [sisVelAll; sisVel(i).observations(2:end,1)];
%         pointStateAll = [pointStateAll; pointState(i).values(1:end-1,1)];
%     end
%     availPoints = find( ~isnan(comVelAll(:,1)) & ~isnan(sisVelAll(:,1)) );
%     comVelSisVelAft = [comVelAll(availPoints,1) sisVelAll(availPoints,1) ...
%         pointStateAll(availPoints,1)];
%         
%     %put observations into bins
% 
%     %find maximum absolute center position
%     maxComPos = max(abs(comPosSisSep(:,1)));
% 
%     %specify bin width and calculate number of bins accordingly
%     binWidth = 0.2; %microns
%     numBins = ceil( maxComPos / binWidth );
% 
%     %go over bins ...
%     for iBin = 1 - numBins : numBins
% 
%         %find bin edges
%         minBinPos = (iBin-1) * binWidth;
%         maxBinPos = iBin * binWidth;
% 
%         %for sister separation
% 
%         %find indices of center position observations within this bin
%         indxBin = find( comPosSisSep(:,1) >= minBinPos & ...
%             comPosSisSep(:,1) <= maxBinPos );
% 
%         %store data in a binned format
%         eval(['comPosSisSep' label{iLabel,1} '(iBin+numBins).minMax = [minBinPos maxBinPos];'])
%         eval(['comPosSisSep' label{iLabel,1} '(iBin+numBins).numObserve = length(indxBin);'])
%         eval(['comPosSisSep' label{iLabel,1} '(iBin+numBins).'...
%             'observations = comPosSisSep(indxBin,:);'])
% 
%         %for change in sister separation starting from a point
% 
%         %find indices of center position observations within this bin
%         indxBin = find( comPosSisVelAft(:,1) >= minBinPos & ...
%             comPosSisVelAft(:,1) <= maxBinPos );
% 
%         %store data in a binned format
%         eval(['comPosSisVelAft' label{iLabel,1} '(iBin+numBins).minMax = [minBinPos maxBinPos];'])
%         eval(['comPosSisVelAft' label{iLabel,1} '(iBin+numBins).numObserve = length(indxBin);'])
%         eval(['comPosSisVelAft' label{iLabel,1} '(iBin+numBins).'...
%             'observations = comPosSisVelAft(indxBin,:);'])
% 
%         %for change in sister separation ending in a point
% 
%         %find indices of center position observations within this bin
%         indxBin = find( comPosSisVelBef(:,1) >= minBinPos & ...
%             comPosSisVelBef(:,1) <= maxBinPos );
% 
%         %store data in a binned format
%         eval(['comPosSisVelBef' label{iLabel,1} '(iBin+numBins).minMax = [minBinPos maxBinPos];'])
%         eval(['comPosSisVelBef' label{iLabel,1} '(iBin+numBins).numObserve = length(indxBin);'])
%         eval(['comPosSisVelBef' label{iLabel,1} '(iBin+numBins).'...
%             'observations = comPosSisVelBef(indxBin,:);'])
% 
%     end
%     
%     %find maximum absolute sister separation
%     maxSisSep = max(abs(sisSepComVelAft(:,1)));
% 
%     %specify bin width and calculate number of bins accordingly
%     binWidth = 0.1; %microns
%     numBins = ceil( maxSisSep / binWidth );
% 
%     %go over bins ...
%     for iBin = 1 - numBins : numBins
% 
%         %find bin edges
%         minBinPos = (iBin-1) * binWidth;
%         maxBinPos = iBin * binWidth;
% 
%         %for change in center position starting from a point
% 
%         %find indices of sister separation observations within this bin
%         indxBin = find( sisSepComVelAft(:,1) >= minBinPos & ...
%             sisSepComVelAft(:,1) <= maxBinPos );
% 
%         %store data in a binned format
%         eval(['sisSepComVelAft' label{iLabel,1} '(iBin+numBins).minMax = [minBinPos maxBinPos];'])
%         eval(['sisSepComVelAft' label{iLabel,1} '(iBin+numBins).numObserve = length(indxBin);'])
%         eval(['sisSepComVelAft' label{iLabel,1} '(iBin+numBins).'...
%             'observations = sisSepComVelAft(indxBin,:);'])
% 
%         %for change in center position ending in a point
% 
%         %find indices of sister separation observations within this bin
%         indxBin = find( sisSepComVelBef(:,1) >= minBinPos & ...
%             sisSepComVelBef(:,1) <= maxBinPos );
% 
%         %store data in a binned format
%         eval(['sisSepComVelBef' label{iLabel,1} '(iBin+numBins).minMax = [minBinPos maxBinPos];'])
%         eval(['sisSepComVelBef' label{iLabel,1} '(iBin+numBins).numObserve = length(indxBin);'])
%         eval(['sisSepComVelBef' label{iLabel,1} '(iBin+numBins).'...
%             'observations = sisSepComVelBef(indxBin,:);'])
% 
%     end
% 
%     %find maximum absolute center position change
%     maxComVel = max(abs(comVelSisVel(:,1)));
% 
%     %specify bin width and calculate number of bins accordingly
%     binWidth = 0.1; %microns
%     numBins = ceil( maxComVel / binWidth );
% 
%     %go over bins ...
%     for iBin = 1 - numBins : numBins
% 
%         %find bin edges
%         minBinPos = (iBin-1) * binWidth;
%         maxBinPos = iBin * binWidth;
% 
%         %for change in sister separation
% 
%         %find indices of center position change observations within this bin
%         indxBin = find( comVelSisVel(:,1) >= minBinPos & ...
%             comVelSisVel(:,1) <= maxBinPos );
% 
%         %store data in a binned format
%         eval(['comVelSisVel' label{iLabel,1} '(iBin+numBins).minMax = [minBinPos maxBinPos];'])
%         eval(['comVelSisVel' label{iLabel,1} '(iBin+numBins).numObserve = length(indxBin);'])
%         eval(['comVelSisVel' label{iLabel,1} '(iBin+numBins).'...
%             'observations = comVelSisVel(indxBin,:);'])
% 
%     end
% 
% end

    %     eval(['distrPerState = struct(''sisterSep'',sisterSepDistr' label{iLabel,1} ...
    %         ',''sisterSepChangeAft'',sisterSepChangeAftDistr' label{iLabel,1} ...
    %         ',''sisterSepChangeBef'',sisterSepChangeBefDistr' label{iLabel,1} ...
    %         ',''centerPos'',comPosDistr' label{iLabel,1} ...
    %         ',''centerPosChangeAft'',comPosChangeAftDistr' label{iLabel,1} ...
    %         ',''centerPosChangeBef'',comPosChangeBefDistr' label{iLabel,1} ');']);
    %     eval(['separationVsPosition = struct(''comPosSisSep'',comPosSisSep' label{iLabel,1} ...
    %         ',''comPosSisVelAft'',comPosSisVelAft' label{iLabel,1} ...
    %         ',''comPosSisVelBef'',comPosSisVelBef' label{iLabel,1} ...
    %         ',''sisSepComVelAft'',sisSepComVelAft' label{iLabel,1} ...
    %         ',''sisSepComVelBef'',sisSepComVelBef' label{iLabel,1} ...
    %         ',''comVelSisVel'',comVelSisVel' label{iLabel,1} ');']);

        %         '''distrPerState'',distrPerState,'...
        %         '''separationVsPosition'',separationVsPosition);'])

            %% distributions per state stuff %%

%     for iLabel = goodLabel
% 
%         %assign figure titles
%         figureInfo(1).title = 'Sister separation';
%         figureInfo(2).title = 'Sister separation change ''after''';
%         figureInfo(3).title = 'Sister separation change ''before''';
%         figureInfo(4).title = 'Center position';
%         figureInfo(5).title = 'Center position change ''after''';
%         figureInfo(6).title = 'Center position change ''before''';
%         
%         %assign variable names
%         figureInfo(1).varName = 'sisterSepDistr';
%         figureInfo(2).varName = 'sisterSepChangeAftDistr';
%         figureInfo(3).varName = 'sisterSepChangeBefDistr';
%         figureInfo(4).varName = 'comPosDistr';
%         figureInfo(5).varName = 'comPosChangeAftDistr';
%         figureInfo(6).varName = 'comPosChangeBefDistr';
%         
%         for iFig = 1 : 6
% 
%             %open figure and write title
%             figure('Name',[fileName(1:end-4) ' - ' figureInfo(iFig).title ' - ' label{iLabel,1}],'NumberTitle','off');
% 
%             %fetch variable to plot
%             eval(['var2plot = ' figureInfo(iFig).varName label{iLabel,1} ';'])
% 
%             %subplot 1: non-center, non-transition points, +ve movement
%             subplot(2,3,1);
%             histogram(var2plot(7).observations(:,1),[],[],0);
%             xlabel('+ve, general');
%             ylabel('Counts');
% 
%             %subplot 2: center points, +ve movement
%             subplot(2,3,2);
%             histogram(var2plot(6).observations(:,1),[],[],0);
%             xlabel('+ve, center');
%             ylabel('Counts');
% 
%             %subplot 3: +ve to -ve transition points
%             subplot(2,3,3);
%             histogram(var2plot(5).observations(:,1),[],[],0);
%             xlabel('+ve to -ve transition');
%             ylabel('Counts');
% 
%             %subplot 4: non-center, non-transition points, -ve movement
%             subplot(2,3,4);
%             histogram(var2plot(1).observations(:,1),[],[],0);
%             xlabel('-ve, general');
%             ylabel('Counts');
% 
%             %subplot 5: center points, -ve movement
%             subplot(2,3,5);
%             histogram(var2plot(2).observations(:,1),[],[],0);
%             xlabel('-ve, center');
%             ylabel('Counts');
% 
%             %subplot 6: -ve to +ve transition points
%             subplot(2,3,6);
%             histogram(var2plot(3).observations(:,1),[],[],0);
%             xlabel('-ve to +ve transition');
%             ylabel('Counts');
%             
%         end
% 
%     end
