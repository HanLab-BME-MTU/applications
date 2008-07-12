function analysisStruct = makiSisterMotionCoupling(jobType,...
    analysisStruct,verbose,samplingPeriod)
%MAKISISTERMOTIONCOUPLING looks for coupling in the motion of sister kinetochores
%
%SYNOPSIS analysisStruct = makiSisterMotionCoupling(jobType,...
%    analysisStruct,verbose,samplingPeriod)
%
%INPUT  jobType       : string which can take the values:
%                       'TEST', 'HERCULES', 'DANUSER', 'MERALDI',
%                       'SWEDLOW' or 'MCAINSH'
%       analysisStruct: Structure with fields
%               .fileName: Name of file where analysisStruct is stored.
%               .filePath: Directory where analysisStruct is stored.
%               .movies  : nx2 cell array of the names of the n selected movies.
%                          First column: file name, second column: file path.
%                       Optional. If not input, GUI to load movies is launched.
%       verbose       : 1 to make plots, 0 otherwise. Optional. Default: 0.
%       samplingPeriod: 1 to keep sampling as is, 2 to downsample by taking
%                       every 2nd time point, 3 to downsample by taking
%                       every 3rd time point, etc.
%                       Optional. Default: 1.
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

if nargin < 4 || isempty(samplingPeriod)
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

%multiply number of sisters by samplingPeriod to get the effective number
%of sisters
numSistersTotEff = numSistersTot * samplingPeriod;

%% collect sister displacements and separations

%define cell array of sister labels
label{1,1} = 'Inlier';
label{2,1} = 'Unaligned';
label{3,1} = 'Lagging';

%initialize global sister index and structures
for iLabel = 1 : 3

    eval(['iGlobal' label{iLabel,1} ' = 0;'])

    eval(['projectionSis1' label{iLabel,1} '(1:numSistersTotEff,1) = struct(''observations'',[]);'])
    eval(['projectionSis2' label{iLabel,1} '(1:numSistersTotEff,1) = struct(''observations'',[]);'])
    eval(['angleSis1' label{iLabel,1} '(1:numSistersTotEff,1) = struct(''observations'',[]);'])
    eval(['angleSis2' label{iLabel,1} '(1:numSistersTotEff,1) = struct(''observations'',[]);'])
    eval(['dispMagSis1' label{iLabel,1} '(1:numSistersTotEff,1) = struct(''observations'',[]);'])
    eval(['dispMagSis2' label{iLabel,1} '(1:numSistersTotEff,1) = struct(''observations'',[]);'])
    eval(['separationSis12' label{iLabel,1} '(1:numSistersTotEff,1) = struct(''observations'',[]);'])
    eval(['sepChangeSis12' label{iLabel,1} '(1:numSistersTotEff,1) = struct(''observations'',[]);'])
    eval(['crossProdSign12' label{iLabel,1} '(1:numSistersTotEff,1) = struct(''observations'',[]);'])
        
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

        %get sister list
        sisterListTmp = dataStruct(iMovie).sisterList;
        
        %get tracks
        tracksTmp = dataStruct(iMovie).tracks;
        [tracksTmp,tracksIndxTmp] = convStruct2MatNoMS(tracksTmp);

        %sister type = 0 if all kinetochores are inliers
        %sister type = 1 if some kinetochores are unaligned
        %sister type = 2 if some kinetochores are lagging
        sisterType = zeros(numSisters(iMovie),1);
        sisterType(dataStruct(iMovie).updatedClass(1).sistersUnaligned) = 1;
        sisterType(dataStruct(iMovie).updatedClass(1).sistersLagging) = 2;

        %subsample as requested
        for iSample = 1 : samplingPeriod
            
            %get sister coordinates with proper subsampling
            for iSister = 1 : numSisters(iMovie)
                sisterList(iSister).distances = sisterListTmp(iSister).distances(iSample:samplingPeriod:end,:);
            end
            
            %get track coordinates with proper subsampling
            frames2keep = (iSample:samplingPeriod:size(tracksTmp,2)/8);
            columns2keep = [8*(frames2keep-1)+1; 8*(frames2keep-1)+2; ...
                8*(frames2keep-1)+3; 8*(frames2keep-1)+4; ...
                8*(frames2keep-1)+5; 8*(frames2keep-1)+6; ...
                8*(frames2keep-1)+7; 8*(frames2keep-1)+8];
            columns2keep = columns2keep(:)';
            tracks = tracksTmp(:,columns2keep);
            trackFeatIndx = tracksIndxTmp(:,frames2keep);
            
            %get track start and end times
            nonemptyIndx = find(any(~isnan(tracks),2));
            trackSEL = NaN(size(tracks,1),3);
            trackSEL(nonemptyIndx,:) = getTrackSEL(tracks(nonemptyIndx,:));

            %copy fields out of dataStruct(iMovie)
            frameAlignment = dataStruct(iMovie).frameAlignment(iSample:samplingPeriod:end);
            updatedClass = dataStruct(iMovie).updatedClass(iSample:samplingPeriod:end);
            numFramesMovie = length(updatedClass);

            %find frame where anaphase starts (if it starts at all)
            framePhase = vertcat(updatedClass.phase);
            firstFrameAna = find(framePhase=='a',1,'first');
            if isempty(firstFrameAna)
                firstFrameAna = numFramesMovie + 1;
            end

            %go over all sisters in movie
            for iSister = 1 : numSisters(iMovie)

                %find track indices
                tracksIndx = sisterListTmp(1).trackPairs(iSister,1:2);

                %determine frame where each track starts
                trackStart = [trackSEL(tracksIndx(1),1) trackSEL(tracksIndx(2),1)];

                %find number of frames and frames where pair "exists"
                goodFrames = ~isnan(sisterList(iSister).distances(:,1));
                numFrames = length(goodFrames);
                goodFrames = find(goodFrames);
                goodFrames = goodFrames(goodFrames < firstFrameAna);

                %find feature indices making up sisters
                sisterIndx1 = NaN(numFrames,1);
                sisterIndx2 = NaN(numFrames,1);
                sisterIndx1(goodFrames) = trackFeatIndx(tracksIndx(1),goodFrames);
                sisterIndx2(goodFrames) = trackFeatIndx(tracksIndx(2),goodFrames);

                %get aligned sister coordinates
                coords1 = NaN(numFrames,6);
                coords2 = NaN(numFrames,6);
                for iFrame = goodFrames'
                    coords1(iFrame,:) = frameAlignment(iFrame).alignedCoord(sisterIndx1(iFrame),:);
                    coords2(iFrame,:) = frameAlignment(iFrame).alignedCoord(sisterIndx2(iFrame),:);
                end

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

                %calculate vector between sisters and normalize it
                sisterVec = [coords2(:,1:3)-coords1(:,1:3) sqrt(coords1(:,4:6).^2+coords2(:,4:6).^2)];

                %calculate sister separation
                sisterSep = sqrt(sum(sisterVec(:,1:3).^2,2));
                sisterSep = [sisterSep sqrt( sum( (sisterVec(:,1:3) .* ...
                    sisterVec(:,4:6)).^2 ,2) )./sisterSep];

                %normalize vector connecting sisters
                sisterVec = sisterVec ./ repmat(sisterSep(:,1),1,6);

                %calculate change in sister separation
                sisterVel = [diff(sisterSep(:,1)) sqrt( sum( [sisterSep(1:end-1,2) ...
                    sisterSep(2:end,2)].^2 ,2) )];

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

                %calculate displacement magnitude
                %add minus sign to those displacements whose projection is
                %negative
                dispMag1 = sqrt(sum(coords1Diff(:,1:3).^2,2));
                dispMag1 = sign(projection1) .* dispMag1;
                dispMag2 = sqrt(sum(coords2Diff(:,1:3).^2,2));
                dispMag2 = sign(projection2) .* dispMag2;

                %calculate angle between displacements and vector connecting sisters
                angle1 = acos(projection1 ./ sqrt(sum(coords1Diff(:,1:3).^2,2))) * 180 / pi;
                angle2 = acos(projection2 ./ sqrt(sum(coords2Diff(:,1:3).^2,2))) * 180 / pi;

                %use the cross product to deduce whether the sisters
                %"displaced" or "tumbled"
                crossProd1 = cross(coords1Diff(:,1:3),sisterVec(1:end-1,1:3),2);
                crossProd2 = cross(coords2Diff(:,1:3),sisterVec(1:end-1,1:3),2);
                cosAngleCrossProd12 = diag(crossProd1 * crossProd2');
                crossProdSign = sign(cosAngleCrossProd12);

                %store sister information based on sister type
                iLabel = sisterType(iSister) + 1;

                %increase global index of sister type by 1
                eval(['iGlobal' label{iLabel,1} ' = iGlobal' label{iLabel,1} ' + 1;'])

                %store information
                eval(['projectionSis1' label{iLabel,1} '(iGlobal' label{iLabel,1} ...
                    ').observations = projection1;']) %um
                eval(['projectionSis2' label{iLabel,1} '(iGlobal' label{iLabel,1} ...
                    ').observations = projection2;']) %um
                eval(['angleSis1' label{iLabel,1} '(iGlobal' label{iLabel,1} ...
                    ').observations = angle1;']) %degrees
                eval(['angleSis2' label{iLabel,1} '(iGlobal' label{iLabel,1} ...
                    ').observations = angle2;']) %degrees
                eval(['dispMagSis1' label{iLabel,1} '(iGlobal' label{iLabel,1} ...
                    ').observations = dispMag1;']) %um
                eval(['dispMagSis2' label{iLabel,1} '(iGlobal' label{iLabel,1} ...
                    ').observations = dispMag2;']) %um
                eval(['separationSis12' label{iLabel,1} '(iGlobal' label{iLabel,1} ...
                    ').observations = sisterSep;']); %um
                eval(['sepChangeSis12' label{iLabel,1} '(iGlobal' label{iLabel,1} ...
                    ').observations = sisterVel;']); %um
                eval(['crossProdSign12' label{iLabel,1} '(iGlobal' label{iLabel,1} ...
                    ').observations = crossProdSign;']);

            end %(for iSister = 1 : numSisters(iMovie) )

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

    eval(['projectionSis1' label{iLabel,1} ' = projectionSis1' label{iLabel,1} ...
        '(1:iGlobal' label{iLabel,1} ');'])
    eval(['projectionSis2' label{iLabel,1} ' = projectionSis2' label{iLabel,1} ...
        '(1:iGlobal' label{iLabel,1} ');'])
    eval(['angleSis1' label{iLabel,1} ' = angleSis1' label{iLabel,1} ...
        '(1:iGlobal' label{iLabel,1} ');'])
    eval(['angleSis2' label{iLabel,1} ' = angleSis2' label{iLabel,1} ...
        '(1:iGlobal' label{iLabel,1} ');'])
    eval(['dispMagSis1' label{iLabel,1} ' = dispMagSis1' label{iLabel,1} ...
        '(1:iGlobal' label{iLabel,1} ');'])
    eval(['dispMagSis2' label{iLabel,1} ' = dispMagSis2' label{iLabel,1} ...
        '(1:iGlobal' label{iLabel,1} ');'])
    eval(['separationSis12' label{iLabel,1} ' = separationSis12' label{iLabel,1} ...
        '(1:iGlobal' label{iLabel,1} ');'])
    eval(['sepChangeSis12' label{iLabel,1} ' = sepChangeSis12' label{iLabel,1} ...
        '(1:iGlobal' label{iLabel,1} ');'])
    eval(['crossProdSign12' label{iLabel,1} ' = crossProdSign12' label{iLabel,1} ...
        '(1:iGlobal' label{iLabel,1} ');'])

end
    
%% distributions and some distribution parameters

%initialization
for iLabel = 1 : 3
    eval(['sister1ProjDistr' label{iLabel,1} ' = [];'])
    eval(['sister1ProjPosParam' label{iLabel,1} ' = [];'])
    eval(['sister1ProjNegParam' label{iLabel,1} ' = [];'])
    eval(['sister1DispDistr' label{iLabel,1} ' = [];'])
    eval(['sister1DispPosParam' label{iLabel,1} ' = [];'])
    eval(['sister1DispNegParam' label{iLabel,1} ' = [];'])
    eval(['sister1AngleDistr' label{iLabel,1} ' = [];'])
    eval(['sister1AngleParam' label{iLabel,1} ' = [];'])
    eval(['sister2ProjDistr' label{iLabel,1} ' = [];'])
    eval(['sister2ProjPosParam' label{iLabel,1} ' = [];'])
    eval(['sister2ProjNegParam' label{iLabel,1} ' = [];'])
    eval(['sister2DispDistr' label{iLabel,1} ' = [];'])
    eval(['sister2DispPosParam' label{iLabel,1} ' = [];'])
    eval(['sister2DispNegParam' label{iLabel,1} ' = [];'])
    eval(['sister2AngleDistr' label{iLabel,1} ' = [];'])
    eval(['sister2AngleParam' label{iLabel,1} ' = [];'])
    eval(['sisterSepBefDistr' label{iLabel,1} ' = [];'])
    eval(['sisterSepBefParam' label{iLabel,1} ' = [];'])
    eval(['sisterSepAftDistr' label{iLabel,1} ' = [];'])
    eval(['sisterSepAftParam' label{iLabel,1} ' = [];'])
    eval(['sisterSepChangeDistr' label{iLabel,1} ' = [];'])
    eval(['sisterSepChangePosParam' label{iLabel,1} ' = [];'])
    eval(['sisterSepChangeNegParam' label{iLabel,1} ' = [];'])
    eval(['crossProdSignDistr' label{iLabel,1} ' = [];'])
    
    eval(['sister1ProjPosIndParam' label{iLabel,1} ' = NaN(numMovies,7);'])
    eval(['sister1ProjNegIndParam' label{iLabel,1} ' = NaN(numMovies,7);'])
    eval(['sister1DispPosIndParam' label{iLabel,1} ' = NaN(numMovies,7);'])
    eval(['sister1DispNegIndParam' label{iLabel,1} ' = NaN(numMovies,7);'])
    eval(['sister1AngleIndParam' label{iLabel,1} ' = NaN(numMovies,7);'])
    eval(['sister2ProjPosIndParam' label{iLabel,1} ' = NaN(numMovies,7);'])
    eval(['sister2ProjNegIndParam' label{iLabel,1} ' = NaN(numMovies,7);'])
    eval(['sister2DispPosIndParam' label{iLabel,1} ' = NaN(numMovies,7);'])
    eval(['sister2DispNegIndParam' label{iLabel,1} ' = NaN(numMovies,7);'])
    eval(['sister2AngleIndParam' label{iLabel,1} ' = NaN(numMovies,7);'])
    eval(['sisterSepBefIndParam' label{iLabel,1} ' = NaN(numMovies,7);'])
    eval(['sisterSepAftIndParam' label{iLabel,1} ' = NaN(numMovies,7);'])
    eval(['sisterSepChangePosIndParam' label{iLabel,1} ' = NaN(numMovies,7);'])
    eval(['sisterSepChangeNegIndParam' label{iLabel,1} ' = NaN(numMovies,7);'])
end

%calculation
for iLabel = goodLabel

    % overall %
    
    %sister 1 projections
    eval(['allValues = vertcat(projectionSis1' label{iLabel,1} '.observations);']);
    allValues = allValues(:,1);
    eval(['sister1ProjDistr' label{iLabel,1} ' = allValues;']);
    eval(['sister1ProjPosParam' label{iLabel,1} ...
        ' = [nanmean(allValues(allValues>0)) nanstd(allValues(allValues>0)) ' ...
        'min(allValues(allValues>0)) prctile(allValues(allValues>0),25) '...
        'prctile(allValues(allValues>0),50) prctile(allValues(allValues>0),75) '...
        'max(allValues(allValues>0))];']);
    eval(['sister1ProjNegParam' label{iLabel,1} ...
        ' = [nanmean(allValues(allValues<0)) nanstd(allValues(allValues<0)) ' ...
        'min(allValues(allValues<0)) prctile(allValues(allValues<0),25) '...
        'prctile(allValues(allValues<0),50) prctile(allValues(allValues<0),75) '...
        'max(allValues(allValues<0))];']);

    %sister 1 displacements
    eval(['allValues = vertcat(dispMagSis1' label{iLabel,1} '.observations);']);
    allValues = allValues(:,1);
    eval(['sister1DispDistr' label{iLabel,1} ' = allValues;']);
    eval(['sister1DispPosParam' label{iLabel,1} ...
        ' = [nanmean(allValues(allValues>0)) nanstd(allValues(allValues>0)) ' ...
        'min(allValues(allValues>0)) prctile(allValues(allValues>0),25) '...
        'prctile(allValues(allValues>0),50) prctile(allValues(allValues>0),75) '...
        'max(allValues(allValues>0))];']);
    eval(['sister1DispNegParam' label{iLabel,1} ...
        ' = [nanmean(allValues(allValues<0)) nanstd(allValues(allValues<0)) ' ...
        'min(allValues(allValues<0)) prctile(allValues(allValues<0),25) '...
        'prctile(allValues(allValues<0),50) prctile(allValues(allValues<0),75) '...
        'max(allValues(allValues<0))];']);

    %sister 1 angles
    eval(['allValues = vertcat(angleSis1' label{iLabel,1} '.observations);']);
    allValues = allValues(:,1);
    eval(['sister1AngleDistr' label{iLabel,1} ' = allValues;']);
    eval(['sister1AngleParam' label{iLabel,1} ...
        ' = [nanmean(allValues) nanstd(allValues) min(allValues) ' ...
        'prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) '...
        'max(allValues)];']);
    
    %sister 2 projections
    eval(['allValues = vertcat(projectionSis2' label{iLabel,1} '.observations);']);
    allValues = allValues(:,1);
    eval(['sister2ProjDistr' label{iLabel,1} ' = allValues;']);
    eval(['sister2ProjPosParam' label{iLabel,1} ...
        ' = [nanmean(allValues(allValues>0)) nanstd(allValues(allValues>0)) ' ...
        'min(allValues(allValues>0)) prctile(allValues(allValues>0),25) '...
        'prctile(allValues(allValues>0),50) prctile(allValues(allValues>0),75) '...
        'max(allValues(allValues>0))];']);
    eval(['sister2ProjNegParam' label{iLabel,1} ...
        ' = [nanmean(allValues(allValues<0)) nanstd(allValues(allValues<0)) ' ...
        'min(allValues(allValues<0)) prctile(allValues(allValues<0),25) '...
        'prctile(allValues(allValues<0),50) prctile(allValues(allValues<0),75) '...
        'max(allValues(allValues<0))];']);

    %sister 2 displacements
    eval(['allValues = vertcat(dispMagSis2' label{iLabel,1} '.observations);']);
    allValues = allValues(:,1);
    eval(['sister2DispDistr' label{iLabel,1} ' = allValues;']);
    eval(['sister2DispPosParam' label{iLabel,1} ...
        ' = [nanmean(allValues(allValues>0)) nanstd(allValues(allValues>0)) ' ...
        'min(allValues(allValues>0)) prctile(allValues(allValues>0),25) '...
        'prctile(allValues(allValues>0),50) prctile(allValues(allValues>0),75) '...
        'max(allValues(allValues>0))];']);
    eval(['sister2DispNegParam' label{iLabel,1} ...
        ' = [nanmean(allValues(allValues<0)) nanstd(allValues(allValues<0)) ' ...
        'min(allValues(allValues<0)) prctile(allValues(allValues<0),25) '...
        'prctile(allValues(allValues<0),50) prctile(allValues(allValues<0),75) '...
        'max(allValues(allValues<0))];']);

    %sister 2 angles
    eval(['allValues = vertcat(angleSis2' label{iLabel,1} '.observations);'])
    allValues = allValues(:,1);
    eval(['sister2AngleDistr' label{iLabel,1} ' = allValues;']);
    eval(['sister2AngleParam' label{iLabel,1} ...
        ' = [nanmean(allValues) nanstd(allValues) min(allValues) ' ...
        'prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) '...
        'max(allValues)];']);
    
    %sister separation, before and after
    eval(['allValues = separationSis12' label{iLabel,1} ';'])
    valuesBef = [];
    valuesAft = [];
    for i = 1 : length(allValues)
        valuesBef = [valuesBef; allValues(i).observations(1:end-1,1)];
        valuesAft = [valuesAft; allValues(i).observations(2:end,1)];
    end
    eval(['sisterSepBefDistr' label{iLabel,1} ' = valuesBef;']);
    eval(['sisterSepAftDistr' label{iLabel,1} ' = valuesAft;']);
    eval(['sisterSepBefParam' label{iLabel,1} ...
        ' = [nanmean(valuesBef) nanstd(valuesBef) min(valuesBef) ' ...
        'prctile(valuesBef,25) prctile(valuesBef,50) prctile(valuesBef,75) '...
        'max(valuesBef)];']);
    eval(['sisterSepAftParam' label{iLabel,1} ...
        ' = [nanmean(valuesAft) nanstd(valuesAft) min(valuesAft) ' ...
        'prctile(valuesAft,25) prctile(valuesAft,50) prctile(valuesAft,75) '...
        'max(valuesAft)];']);

    %sister separation change
    eval(['allValues = vertcat(sepChangeSis12' label{iLabel,1} '.observations);']);
    allValues = allValues(:,1);
    eval(['sisterSepChangeDistr' label{iLabel,1} ' = allValues;']);
    eval(['sisterSepChangePosParam' label{iLabel,1} ...
        ' = [nanmean(allValues(allValues>0)) nanstd(allValues(allValues>0)) ' ...
        'min(allValues(allValues>0)) prctile(allValues(allValues>0),25) '...
        'prctile(allValues(allValues>0),50) prctile(allValues(allValues>0),75) '...
        'max(allValues(allValues>0))];']);
    eval(['sisterSepChangeNegParam' label{iLabel,1} ...
        ' = [nanmean(allValues(allValues<0)) nanstd(allValues(allValues<0)) ' ...
        'min(allValues(allValues<0)) prctile(allValues(allValues<0),25) '...
        'prctile(allValues(allValues<0),50) prctile(allValues(allValues<0),75) '...
        'max(allValues(allValues<0))];']);

    %sign of dot products of cross products
    eval(['allValues = vertcat(crossProdSign12' label{iLabel,1} '.observations);']);
    allValues = allValues(:,1);
    eval(['crossProdSignDistr' label{iLabel,1} ' = allValues;']);
    
    % individual cells %
    
    %sister 1 projections
    for iMovie = 1 : numMovies
        eval(['allValues = vertcat(projectionSis1' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie)).observations);']);
        if ~isempty(allValues)
            allValues = allValues(:,1);
            eval(['sister1ProjPosIndParam' label{iLabel,1} '(iMovie,:)' ...
                ' = [nanmean(allValues(allValues>0)) nanstd(allValues(allValues>0)) ' ...
                'min(allValues(allValues>0)) prctile(allValues(allValues>0),25) '...
                'prctile(allValues(allValues>0),50) prctile(allValues(allValues>0),75) '...
                'max(allValues(allValues>0))];']);
            eval(['sister1ProjNegIndParam' label{iLabel,1} '(iMovie,:)' ...
                ' = [nanmean(allValues(allValues<0)) nanstd(allValues(allValues<0)) ' ...
                'min(allValues(allValues<0)) prctile(allValues(allValues<0),25) '...
                'prctile(allValues(allValues<0),50) prctile(allValues(allValues<0),75) '...
                'max(allValues(allValues<0))];']);
        end
    end
    
    %sister 1 displacements
    for iMovie = 1 : numMovies
        eval(['allValues = vertcat(dispMagSis1' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie)).observations);']);
        if ~isempty(allValues)
            allValues = allValues(:,1);
            eval(['sister1DispPosIndParam' label{iLabel,1} '(iMovie,:)' ...
                ' = [nanmean(allValues(allValues>0)) nanstd(allValues(allValues>0)) ' ...
                'min(allValues(allValues>0)) prctile(allValues(allValues>0),25) '...
                'prctile(allValues(allValues>0),50) prctile(allValues(allValues>0),75) '...
                'max(allValues(allValues>0))];']);
            eval(['sister1DispNegIndParam' label{iLabel,1} '(iMovie,:)' ...
                ' = [nanmean(allValues(allValues<0)) nanstd(allValues(allValues<0)) ' ...
                'min(allValues(allValues<0)) prctile(allValues(allValues<0),25) '...
                'prctile(allValues(allValues<0),50) prctile(allValues(allValues<0),75) '...
                'max(allValues(allValues<0))];']);
        end
    end
    
    %sister 1 angles
    for iMovie = 1 : numMovies
        eval(['allValues = vertcat(angleSis1' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie)).observations);']);
        if ~isempty(allValues)
            allValues = allValues(:,1);
            eval(['sister1AngleIndParam' label{iLabel,1} '(iMovie,:)' ...
                ' = [nanmean(allValues) nanstd(allValues) min(allValues) ' ...
                'prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) '...
                'max(allValues)];']);
        end
    end

    %sister 2 projections
    for iMovie = 1 : numMovies
        eval(['allValues = vertcat(projectionSis2' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie)).observations);']);
        if ~isempty(allValues)
            allValues = allValues(:,1);
            eval(['sister2ProjPosIndParam' label{iLabel,1} '(iMovie,:)' ...
                ' = [nanmean(allValues(allValues>0)) nanstd(allValues(allValues>0)) ' ...
                'min(allValues(allValues>0)) prctile(allValues(allValues>0),25) '...
                'prctile(allValues(allValues>0),50) prctile(allValues(allValues>0),75) '...
                'max(allValues(allValues>0))];']);
            eval(['sister2ProjNegIndParam' label{iLabel,1} '(iMovie,:)' ...
                ' = [nanmean(allValues(allValues<0)) nanstd(allValues(allValues<0)) ' ...
                'min(allValues(allValues<0)) prctile(allValues(allValues<0),25) '...
                'prctile(allValues(allValues<0),50) prctile(allValues(allValues<0),75) '...
                'max(allValues(allValues<0))];']);
        end
    end

    %sister 2 displacements
    for iMovie = 1 : numMovies
        eval(['allValues = vertcat(dispMagSis2' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie)).observations);']);
        if ~isempty(allValues)
            allValues = allValues(:,1);
            eval(['sister2DispPosIndParam' label{iLabel,1} '(iMovie,:)' ...
                ' = [nanmean(allValues(allValues>0)) nanstd(allValues(allValues>0)) ' ...
                'min(allValues(allValues>0)) prctile(allValues(allValues>0),25) '...
                'prctile(allValues(allValues>0),50) prctile(allValues(allValues>0),75) '...
                'max(allValues(allValues>0))];']);
            eval(['sister2DispNegIndParam' label{iLabel,1} '(iMovie,:)' ...
                ' = [nanmean(allValues(allValues<0)) nanstd(allValues(allValues<0)) ' ...
                'min(allValues(allValues<0)) prctile(allValues(allValues<0),25) '...
                'prctile(allValues(allValues<0),50) prctile(allValues(allValues<0),75) '...
                'max(allValues(allValues<0))];']);
        end
    end
    
    %sister 2 angles
    for iMovie = 1 : numMovies
        eval(['allValues = vertcat(angleSis2' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie)).observations);']);
        if ~isempty(allValues)
            allValues = allValues(:,1);
            eval(['sister2AngleIndParam' label{iLabel,1} '(iMovie,:)' ...
                ' = [nanmean(allValues) nanstd(allValues) min(allValues) ' ...
                'prctile(allValues,25) prctile(allValues,50) prctile(allValues,75) '...
                'max(allValues)];']);
        end
    end
        
    %sister separation, before and after
    for iMovie = 1 : numMovies
        eval(['allValues = separationSis12' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie));']);
        valuesBef = [];
        valuesAft = [];
        for i = 1 : length(allValues)
            valuesBef = [valuesBef; allValues(i).observations(1:end-1,1)];
            valuesAft = [valuesAft; allValues(i).observations(2:end,1)];
        end
        if ~isempty(valuesBef)
            eval(['sisterSepBefIndParam' label{iLabel,1} '(iMovie,:)' ...
                ' = [nanmean(valuesBef) nanstd(valuesBef) min(valuesBef) ' ...
                'prctile(valuesBef,25) prctile(valuesBef,50) prctile(valuesBef,75) '...
                'max(valuesBef)];']);
        end
        if ~isempty(valuesAft)
            eval(['sisterSepAftIndParam' label{iLabel,1} '(iMovie,:)' ...
                ' = [nanmean(valuesAft) nanstd(valuesAft) min(valuesAft) ' ...
                'prctile(valuesAft,25) prctile(valuesAft,50) prctile(valuesAft,75) '...
                'max(valuesAft)];']);
        end
    end

    %sister separation change
    for iMovie = 1 : numMovies
        eval(['allValues = vertcat(sepChangeSis12' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie)).observations);']);
        if ~isempty(allValues)
            allValues = allValues(:,1);
            eval(['sisterSepChangePosIndParam' label{iLabel,1} '(iMovie,:)' ...
                ' = [nanmean(allValues(allValues>0)) nanstd(allValues(allValues>0)) ' ...
                'min(allValues(allValues>0)) prctile(allValues(allValues>0),25) '...
                'prctile(allValues(allValues>0),50) prctile(allValues(allValues>0),75) '...
                'max(allValues(allValues>0))];']);
            eval(['sisterSepChangeNegIndParam' label{iLabel,1} '(iMovie,:)' ...
                ' = [nanmean(allValues(allValues<0)) nanstd(allValues(allValues<0)) ' ...
                'min(allValues(allValues<0)) prctile(allValues(allValues<0),25) '...
                'prctile(allValues(allValues<0),50) prctile(allValues(allValues<0),75) '...
                'max(allValues(allValues<0))];']);
        end
    end

end

%% cross-correlation of motion

%define maximum lag
maxLag = round(20 / samplingPeriod);
maxLagSis = round(10 / samplingPeriod);

%initialization
for iLabel = 1 : 3
    eval(['projectionCrosscorr' label{iLabel,1} ' = [];'])
    eval(['angleCrosscorr' label{iLabel,1} ' = [];'])
    eval(['projection1Autocorr' label{iLabel,1} ' = [];'])
    eval(['angle1Autocorr' label{iLabel,1} ' = [];'])
    eval(['projection2Autocorr' label{iLabel,1} ' = [];'])
    eval(['angle2Autocorr' label{iLabel,1} ' = [];'])
    
    eval(['projectionIndCrosscorr' label{iLabel,1} ' = NaN(2*maxLag+1,2,numMovies);'])
    eval(['angleIndCrosscorr' label{iLabel,1} ' = NaN(2*maxLag+1,2,numMovies);'])
    eval(['projection1IndAutocorr' label{iLabel,1} ' = NaN(maxLag+1,2,numMovies);'])
    eval(['angle1IndAutocorr' label{iLabel,1} ' = NaN(maxLag+1,2,numMovies);'])
    eval(['projection2IndAutocorr' label{iLabel,1} ' = NaN(maxLag+1,2,numMovies);'])
    eval(['angle2IndAutocorr' label{iLabel,1} ' = NaN(maxLag+1,2,numMovies);'])
    
    eval(['projectionSisCrosscorr' label{iLabel,1} ' = NaN(2*maxLagSis+1,2,iGlobal' label{iLabel,1} ');'])
    eval(['angleSisCrosscorr' label{iLabel,1} ' = NaN(2*maxLagSis+1,2,iGlobal' label{iLabel,1} ');'])
    eval(['projection1SisAutocorr' label{iLabel,1} ' = NaN(maxLagSis+1,2,iGlobal' label{iLabel,1} ');'])
    eval(['angle1SisAutocorr' label{iLabel,1} ' = NaN(maxLagSis+1,2,iGlobal' label{iLabel,1} ');'])
    eval(['projection2SisAutocorr' label{iLabel,1} ' = NaN(maxLagSis+1,2,iGlobal' label{iLabel,1} ');'])
    eval(['angle2SisAutocorr' label{iLabel,1} ' = NaN(maxLagSis+1,2,iGlobal' label{iLabel,1} ');'])
end

%calculation
for iLabel = goodLabel
    
    % overall %

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

    % individual cells %

    %between sisters

    %projections
    for iMovie = 1 : numMovies
        eval(['traj1 = projectionSis1' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie));'])
        eval(['traj2 = projectionSis2' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie));'])
        if ~isempty(traj1) && ~isempty(traj2)
            [tmpCorr,errFlag] = crossCorr(traj1,traj2,maxLag);
            if ~errFlag
                eval(['projectionIndCrosscorr' label{iLabel,1} '(:,:,iMovie) = tmpCorr;'])
            end
        end
    end

    %angles
    for iMovie = 1 : numMovies
        eval(['traj1 = angleSis1' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie));'])
        eval(['traj2 = angleSis2' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie));'])
        if ~isempty(traj1) && ~isempty(traj2)
            [tmpCorr,errFlag] = crossCorr(traj1,traj2,maxLag);
            if ~errFlag
                eval(['angleIndCrosscorr' label{iLabel,1} '(:,:,iMovie) = tmpCorr;'])
            end
        end
    end

    %sister with itself

    %projections
    for iMovie = 1 : numMovies
        eval(['traj = projectionSis1' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie));'])
        if ~isempty(traj)
            [tmpCorr,errFlag] = autoCorr(traj,maxLag);
            if ~errFlag
                eval(['projection1IndAutocorr' label{iLabel,1} '(:,:,iMovie) = tmpCorr;'])
            end
        end
    end

    %angles
    for iMovie = 1 : numMovies
        eval(['traj = angleSis1' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie));'])
        if ~isempty(traj)
            [tmpCorr,errFlag] = autoCorr(traj,maxLag);
            if ~errFlag
                eval(['angle1IndAutocorr' label{iLabel,1} '(:,:,iMovie) = tmpCorr;'])
            end
        end
    end

    %projections
    for iMovie = 1 : numMovies
        eval(['traj = projectionSis2' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie));'])
        if ~isempty(traj)
            [tmpCorr,errFlag] = autoCorr(traj,maxLag);
            if ~errFlag
                eval(['projection2IndAutocorr' label{iLabel,1} '(:,:,iMovie) = tmpCorr;'])
            end
        end
    end

    %angles
    for iMovie = 1 : numMovies
        eval(['traj = angleSis2' label{iLabel,1} '(movieStartIndx' ...
            label{iLabel,1} '(iMovie):movieEndIndx' label{iLabel,1} '(iMovie));'])
        if ~isempty(traj)
            [tmpCorr,errFlag] = autoCorr(traj,maxLag);
            if ~errFlag
                eval(['angle2IndAutocorr' label{iLabel,1} '(:,:,iMovie) = tmpCorr;'])
            end
        end
    end

    % individual sisters %

    %between sisters
    
    %projections
    for iSister = 1 : eval(['iGlobal' label{iLabel,1}])
        eval(['traj1 = projectionSis1' label{iLabel,1} '(iSister).observations;'])
        eval(['traj2 = projectionSis2' label{iLabel,1} '(iSister).observations;'])
        if length(find(~isnan(traj1(:,1))))>maxLagSis+10
            [tmpCorr,errFlag] = crossCorr(traj1,traj2,maxLagSis);
            if ~errFlag
                eval(['projectionSisCrosscorr' label{iLabel,1} '(:,:,iSister) = tmpCorr;'])
            end
        end
    end

    %angles
    for iSister = 1 : eval(['iGlobal' label{iLabel,1}])
        eval(['traj1 = angleSis1' label{iLabel,1} '(iSister).observations;'])
        eval(['traj2 = angleSis2' label{iLabel,1} '(iSister).observations;'])
        if length(find(~isnan(traj1(:,1))))>maxLagSis+10
            [tmpCorr,errFlag] = crossCorr(traj1,traj2,maxLagSis);
            if ~errFlag
                eval(['angleSisCrosscorr' label{iLabel,1} '(:,:,iSister) = tmpCorr;'])
            end
        end
    end

    %sister with itself

    %projections
    for iSister = 1 : eval(['iGlobal' label{iLabel,1}])
        eval(['traj = projectionSis1' label{iLabel,1} '(iSister).observations;'])
        if length(find(~isnan(traj(:,1))))>maxLagSis+10
            [tmpCorr,errFlag] = autoCorr(traj,maxLagSis);
            if ~errFlag
                eval(['projection1SisAutocorr' label{iLabel,1} '(:,:,iSister) = tmpCorr;'])
            end
        end
    end

    %angles
    for iMovie = 1 : eval(['iGlobal' label{iLabel,1}])
        eval(['traj = angleSis1' label{iLabel,1} '(iSister).observations;'])
        if length(find(~isnan(traj(:,1))))>maxLagSis+10
            [tmpCorr,errFlag] = autoCorr(traj,maxLag);
            if ~errFlag
                eval(['angle1SisAutocorr' label{iLabel,1} '(:,:,iSister) = tmpCorr;'])
            end
        end
    end

    %projections
    for iSister = 1 : eval(['iGlobal' label{iLabel,1}])
        eval(['traj = projectionSis2' label{iLabel,1} '(iSister).observations;'])
        if length(find(~isnan(traj(:,1))))>maxLagSis+10
            [tmpCorr,errFlag] = autoCorr(traj,maxLagSis);
            if ~errFlag
                eval(['projection2SisAutocorr' label{iLabel,1} '(:,:,iSister) = tmpCorr;'])
            end
        end
    end

    %angles
    for iMovie = 1 : eval(['iGlobal' label{iLabel,1}])
        eval(['traj = angleSis2' label{iLabel,1} '(iSister).observations;'])
        if length(find(~isnan(traj(:,1))))>maxLagSis+10
            [tmpCorr,errFlag] = autoCorr(traj,maxLag);
            if ~errFlag
                eval(['angle2SisAutocorr' label{iLabel,1} '(:,:,iSister) = tmpCorr;'])
            end
        end
    end

end

%% ARMAX

%this part will be done later - ARMAX still needs some work

%% output to analysisStruct

for iLabel = 1 : 3

    %store results in one structure
    eval(['distribution = struct('...
        '''projectionSis1'',sister1ProjDistr' label{iLabel,1} ','...
        '''signedDispMagSis1'',sister1DispDistr' label{iLabel,1} ','...
        '''angleSis1'',sister1AngleDistr' label{iLabel,1} ','...
        '''projectionSis2'',sister2ProjDistr' label{iLabel,1} ','...
        '''signedDispMagSis2'',sister2DispDistr' label{iLabel,1} ','...
        '''angleSis2'',sister2AngleDistr' label{iLabel,1} ',' ...
        '''sisterSepBefore'',sisterSepBefDistr' label{iLabel,1} ','...
        '''sisterSepAfter'',sisterSepAftDistr' label{iLabel,1} ','...
        '''sisterSepChange'',sisterSepChangeDistr' label{iLabel,1} ',' ...
        '''signDotProdCrossProd'',crossProdSignDistr' label{iLabel,1} ');']);
    eval(['meanStdMin25P50P75PMax.all = struct('...
        '''projectionSis1Pos'',sister1ProjPosParam' label{iLabel,1} ','...
        '''projectionSis1Neg'',sister1ProjNegParam' label{iLabel,1} ','...
        '''dispMagSis1Pos'',sister1DispPosParam' label{iLabel,1} ','...
        '''dispMagSis1Neg'',sister1DispNegParam' label{iLabel,1} ','...
        '''angleSis1'',sister1AngleParam' label{iLabel,1} ','...
        '''projectionSis2Pos'',sister2ProjPosParam' label{iLabel,1} ','...
        '''projectionSis2Neg'',sister2ProjNegParam' label{iLabel,1} ','...
        '''dispMagSis2Pos'',sister2DispPosParam' label{iLabel,1} ','...
        '''dispMagSis2Neg'',sister2DispNegParam' label{iLabel,1} ','...
        '''angleSis2'',sister2AngleParam' label{iLabel,1} ',' ...
        '''sisterSepBefore'',sisterSepBefParam' label{iLabel,1} ','...
        '''sisterSepAfter'',sisterSepAftParam' label{iLabel,1} ','...
        '''sisterSepChangePos'',sisterSepChangePosParam' label{iLabel,1} ','...
        '''sisterSepChangeNeg'',sisterSepChangeNegParam' label{iLabel,1} ');']);
    eval(['meanStdMin25P50P75PMax.indcell = struct('...
        '''projectionSis1Pos'',sister1ProjPosIndParam' label{iLabel,1} ','...
        '''projectionSis1Neg'',sister1ProjNegIndParam' label{iLabel,1} ','...
        '''dispMagSis1Pos'',sister1DispPosIndParam' label{iLabel,1} ','...
        '''dispMagSis1Neg'',sister1DispNegIndParam' label{iLabel,1} ','...
        '''angleSis1'',sister1AngleIndParam' label{iLabel,1} ','...
        '''projectionSis2Pos'',sister2ProjPosIndParam' label{iLabel,1} ','...
        '''projectionSis2Neg'',sister2ProjNegIndParam' label{iLabel,1} ','...
        '''dispMagSis2Pos'',sister2DispPosIndParam' label{iLabel,1} ','...
        '''dispMagSis2Neg'',sister2DispNegIndParam' label{iLabel,1} ','...
        '''angleSis2'',sister2AngleIndParam' label{iLabel,1} ',' ...
        '''sisterSepBefore'',sisterSepBefIndParam' label{iLabel,1} ','...
        '''sisterSepAfter'',sisterSepAftIndParam' label{iLabel,1} ','...
        '''sisterSepChangePos'',sisterSepChangePosIndParam' label{iLabel,1} ','...
        '''sisterSepChangeNeg'',sisterSepChangeNegIndParam' label{iLabel,1} ');']);
    eval(['crosscorr.all = struct(''projections'',projectionCrosscorr' label{iLabel,1} ','...
        '''angles'',angleCrosscorr' label{iLabel,1} ');']);
    eval(['crosscorr.indcell = struct(''projections'',projectionIndCrosscorr' label{iLabel,1} ','...
        '''angles'',angleIndCrosscorr' label{iLabel,1} ');']);
    eval(['crosscorr.indsis = struct(''projections'',projectionSisCrosscorr' label{iLabel,1} ','...
        '''angles'',angleSisCrosscorr' label{iLabel,1} ');']);
    eval(['autocorr.all = struct(''projectionsSis1'',projection1Autocorr' label{iLabel,1} ...
        ',''anglesSis1'',angle1Autocorr' label{iLabel,1} ...
        ',''projectionsSis2'',projection2Autocorr' label{iLabel,1} ...
        ',''anglesSis2'',angle2Autocorr' label{iLabel,1} ');']);
    eval(['autocorr.indcell = struct(''projectionsSis1'',projection1IndAutocorr' label{iLabel,1} ...
        ',''anglesSis1'',angle1IndAutocorr' label{iLabel,1} ...
        ',''projectionsSis2'',projection2IndAutocorr' label{iLabel,1} ...
        ',''anglesSis2'',angle2IndAutocorr' label{iLabel,1} ');']);
    eval(['autocorr.indsis = struct(''projectionsSis1'',projection1SisAutocorr' label{iLabel,1} ...
        ',''anglesSis1'',angle1SisAutocorr' label{iLabel,1} ...
        ',''projectionsSis2'',projection2SisAutocorr' label{iLabel,1} ...
        ',''anglesSis2'',angle2SisAutocorr' label{iLabel,1} ');']);
    eval(['numSistersCat = iGlobal' label{iLabel,1} ';']);

    eval([label{iLabel,1} ' = struct('...
        '''numSisters'',numSistersCat,'...
        '''distribution'',distribution,'...
        '''meanStdMin25P50P75PMax'',meanStdMin25P50P75PMax,'...
        '''crosscorr'',crosscorr,'...
        '''autocorr'',autocorr);'])

end

inputParam = struct('samplingPeriod',samplingPeriod);
sisterMotionCoupling = struct('Inlier',Inlier,'Unaligned',Unaligned,...
    'Lagging',Lagging,'inputParam',inputParam);

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
    timeLapse = round(2*dataStruct(1).dataProperties.timeLapse)/2 * samplingPeriod;
    
    for iLabel = goodLabel

        %open figure and write title
        figure('Name',[fileName(1:end-4) ' - Motion coupling - ' label{iLabel,1}],'NumberTitle','off');

        %create subplot 1
        subplot(2,3,1);
        hold on;

        %plot projection crosscorrelation
        eval(['projectionCrosscorr = projectionCrosscorr' label{iLabel,1} ';']);
        eval(['autocorr1 = projection1Autocorr' label{iLabel,1} ';']);
        eval(['autocorr2 = projection2Autocorr' label{iLabel,1} ';']);
        plot((-maxLag:maxLag)*timeLapse,projectionCrosscorr(:,1),'k','marker','.');
        plot([-maxLag maxLag]*timeLapse,[0 0],'k--');
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
        subplot(2,3,2);
        hold on;
        
        %plot individual movie projection crosscorrelation
        eval(['tmpCorr = squeeze(projectionIndCrosscorr' label{iLabel,1} '(:,1,:));'])
        plot((-maxLag:maxLag)*timeLapse,tmpCorr,'marker','.');
        plot([-maxLag maxLag]*timeLapse,[0 0],'k--');
        
        %set axes limits
        minVal = min(tmpCorr(:));
        axis([-maxLag*timeLapse maxLag*timeLapse min(0,1.1*minVal) 1.1]);
        
        %write axes labels
        xlabel('Lag (s)');
        ylabel('Projection crosscorrelation - ind movies');
        
        %hold off figure
        hold off
        
        %create subplot 3
        subplot(2,3,3);
        hold on;
        
        %plot individual movie projection autocorrelations
        eval(['tmpCorr1 = squeeze(projection1IndAutocorr' label{iLabel,1} '(:,1,:));'])
        eval(['tmpCorr2 = squeeze(projection2IndAutocorr' label{iLabel,1} '(:,1,:));'])
        plot((-maxLag:0)*timeLapse,tmpCorr1(end:-1:1,:),'marker','.');
        plot((0:maxLag)*timeLapse,tmpCorr2,'marker','.');
        plot([-maxLag maxLag]*timeLapse,[0 0],'k--');
        
        %set axes limits
        minVal = min([tmpCorr1(:); tmpCorr2(:)]);
        axis([-maxLag*timeLapse maxLag*timeLapse min(0,1.1*minVal) 1.1]);
        
        %write axes labels
        xlabel('Lag (s)');
        ylabel('Projection autocorrelations - ind movies');
        
        %hold off figure
        hold off

        %create subplot 4
        subplot(2,3,4);
        hold on;

        %plot angle crosscorrelation
        eval(['angleCrosscorr = angleCrosscorr' label{iLabel,1} ';']);
        eval(['autocorr1 = angle1Autocorr' label{iLabel,1} ';']);
        eval(['autocorr2 = angle2Autocorr' label{iLabel,1} ';']);
        plot((-maxLag:maxLag)*timeLapse,angleCrosscorr(:,1),'k','marker','.');
        plot([-maxLag maxLag]*timeLapse,[0 0],'k--');
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
        
        %create subplot 5
        subplot(2,3,5);
        hold on;
        
        %plot individual movie angle crosscorrelation
        eval(['tmpCorr = squeeze(angleIndCrosscorr' label{iLabel,1} '(:,1,:));'])
        plot((-maxLag:maxLag)*timeLapse,tmpCorr,'marker','.');
        plot([-maxLag maxLag]*timeLapse,[0 0],'k--');
        
        %set axes limits
        minVal = min(tmpCorr(:));
        axis([-maxLag*timeLapse maxLag*timeLapse min(0,1.1*minVal) 1.1]);
        
        %write axes labels
        xlabel('Lag (s)');
        ylabel('Angle crosscorrelation - ind movies');
        
        %hold off figure
        hold off
        
        %create subplot 6
        subplot(2,3,6);
        hold on;
        
        %plot individual movie projection autocorrelations
        eval(['tmpCorr1 = squeeze(angle1IndAutocorr' label{iLabel,1} '(:,1,:));'])
        eval(['tmpCorr2 = squeeze(angle2IndAutocorr' label{iLabel,1} '(:,1,:));'])
        plot((-maxLag:0)*timeLapse,tmpCorr1(end:-1:1,:),'marker','.');
        plot((0:maxLag)*timeLapse,tmpCorr2,'marker','.');
        plot([-maxLag maxLag]*timeLapse,[0 0],'k--');
        
        %set axes limits
        minVal = min([tmpCorr1(:); tmpCorr2(:)]);
        axis([-maxLag*timeLapse maxLag*timeLapse min(0,1.1*minVal) 1.1]);
        
        %write axes labels
        xlabel('Lag (s)');
        ylabel('Angle autocorrelations - ind movies');
        
        %hold off figure
        hold off

    end
    
end

%% ~~~ the end ~~~ %%

